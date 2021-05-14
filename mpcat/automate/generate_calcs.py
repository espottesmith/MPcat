# coding: utf-8

import os
import subprocess
from typing import Optional, List, Dict, Union
import datetime
from pathlib import Path
import time
import random

from pymongo import UpdateOne

from monty.serialization import dumpfn

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph

from schrodinger.application.jaguar.autots_exceptions import UnsupportedReaction

from mpcat.aggregate.database import CatDB
from mpcat.apprehend.autots_input import AutoTSSet
from mpcat.utils.comparison import compositions_equal
from mpcat.utils.generate import mol_to_mol_graph


class AutoTSJob:
    """
    A helper class to prepare and execute AutoTS calculations.
    """

    def __init__(self,
                 reactants: List[Union[Molecule, MoleculeGraph]],
                 products: List[Union[Molecule, MoleculeGraph]],
                 path: Union[str, Path],
                 schrodinger_dir: Optional[str] = "SCHRODINGER",
                 job_name: Optional[str] = None,
                 num_cores: Optional[int] = 40,
                 host: Optional[str] = None,
                 save_scratch: Optional[bool] = False,
                 input_params: Optional[Dict] = None):

        """
        Args:
            reactants (list of Molecule objects): the reactants of the reaction
                to be examined
            products (list of Molecule objects): the products of the reaction
                to be examined
            path (str): The directory where this calculation will take place.
            schrodinger_dir (str): A path to the Schrodinger Suite of software.
                This is used to call AutoTS and other utilities. By default,
                this is "$SCHRODINGER", which should be an environment variable
                set at the time of installation.
            job_name (str): If provided (default None), this will set the name
                of the job in Schrodinger's jobcontrol system.
            num_cores (int): How many cores should the program be parallelized
                over (default 40). When multiple subjobs need to be run
                simultaneously, AutoTS will distribute these cores automatically
                between subjobs
            host (str): Which host should the calculation be run on? By default,
                this is "localhost", which should generally mean that the
                calculation is run on the current node without using a queueing
                system
            save_scratch (bool): If True (default False), save a *.zip file
                containing the contents of the calculation scratch directory
            input_params (dict): Keywords and associated values to be provided
                to AutoTSSet
        """

        self.reactants = list()
        self.products = list()

        for reactant in reactants:
            self.reactants.append(mol_to_mol_graph(reactant))

        for product in products:
            self.products.append(mol_to_mol_graph(product))

        if isinstance(path, Path):
            self.path = path
        else:
            self.path = Path(path)

        if schrodinger_dir == "SCHRODINGER":
            self.schrodinger_dir = Path(os.environ["SCHRODINGER"])
        else:
            if isinstance(schrodinger_dir, str):
                self.schrodinger_dir = Path(schrodinger_dir)
            else:
                self.schrodinger_dir = schrodinger_dir

        self.job_name = job_name
        self.num_cores = num_cores
        self.host = host
        self.save_scratch = save_scratch

        self.input_params = input_params

    def setup_calculation(self):
        """
        Prepare for calculation - prepare directory with all necessary input
            files.

        Args:
            None

        Returns:
            None
        """

        if not self.path.exists():
            self.path.mkdir()

        calc_set = AutoTSSet(self.reactants, self.products,
                             **self.input_params)

        calc_set.write(self.path / "autots.in", write_molecules=True,
                       jobname=self.job_name)

    def run(self, command_line_args: Optional[Dict] = None):
        """
        Execute the AutoTS job

        Args:
            command_line_args (dict): A dictionary of flag: value pairs to be
            provided to the autots command-line interface. Ex:
            {"WAIT": None,
             "subdir": None,
             "nsubjobs": 5}
             will be interpreted as "-WAIT -subdir -nsubjobs 5"

        Returns:
            None
        """

        cur_dir = Path.cwd()
        os.chdir(self.path.as_posix())

        command = [(self.schrodinger_dir / "autots").as_posix(),
                   "-PARALLEL", str(self.num_cores), "-use_one_node"]

        if self.host is not None:
            command.append("-HOST")
            command.append(self.host)

        if self.job_name is not None:
            command.append("-jobname")
            command.append(str(self.job_name))

        if self.save_scratch:
            command.append("-SAVE")

        if command_line_args is not None:
            for key, value in command_line_args.items():
                command.append("-" + str(key))

                if value is not None:
                    command.append(str(value))

        command.append("autots.in")

        process = subprocess.run(command)

        os.chdir(cur_dir.as_posix())

        if process.returncode != 0:
            raise RuntimeError("Job launch failed!")


def launch_mass_jobs(reactions: List[Dict[str, List[Union[Molecule, MoleculeGraph]]]],
                     base_dir: Union[str, Path],
                     schrodinger_dir: Optional[Union[str, Path]] = "$SCHRODINGER",
                     job_name_prefix: Optional[str] = None,
                     num_cores: Optional[int] = 40,
                     host: Optional[str] = "localhost",
                     save_scratch: Optional[bool] = False,
                     input_params: Optional[Dict] = None,
                     command_line_args: Optional[Dict] = None):

    """
    Create many AutoTSJobs from
        pymatgen.reaction_network.reaction_network.Reaction objects

    Args:
        reactions (list of dicts of Molecules): Dict, with keys "reactants" and
            "products", and values being lists of Molecules
        base_dir (str): Root directory where all calculation directories should
            be made
        schrodinger_dir (str): A path to the Schrodinger Suite of software.
            This is used to call AutoTS and other utilities. By default,
            this is "$SCHRODINGER", which should be an environment variable
            set at the time of installation.
        job_name_prefix (str): All jobs in this set of reactions will be given a
            unique name, but this prefix will be prepended to all calculations
            in this set.
        num_cores (int): How many cores should the program be parallelized
            over (default 40). When multiple subjobs need to be run
            simultaneously, AutoTS will distribute these cores automatically
            between subjobs
        host (str): Which host should the calculation be run on? By default,
            this is "localhost", which should generally mean that the
            calculation is run on the current node without using a queueing
            system
        save_scratch (bool): If True (default False), save a *.zip file
            containing the contents of the calculation scratch directory
        input_params (dict): Keywords and associated values to be provided
            to AutoTSSet
        command_line_args (dict): A dictionary of flag: value pairs to be
            provided to the autots command-line interface. Ex:
            {"WAIT": None,
             "subdir": None,
             "nsubjobs": 5}
             will be interpreted as "-WAIT -subdir -nsubjobs 5"


    Returns:
        None
    """

    if isinstance(base_dir, str):
        base_dir = Path(base_dir)

    jobs = list()

    jobs_made = 0

    # Make an AutoTSJob object for each reaction
    for rr, reaction in enumerate(reactions):

        # Verify that reaction is balanced in terms of charge and species
        charge_rct = 0
        charge_pro = 0

        reactants = list()
        products = list()

        for mol in reaction["reactants"]:
            if isinstance(mol, Molecule):
                charge_rct += mol.charge
                reactants.append(mol)
            else:
                charge_rct += mol.molecule.charge
                reactants.append(mol.molecule)

        for mol in reaction["products"]:
            if isinstance(mol, Molecule):
                charge_pro += mol.charge
                products.append(mol)
            else:
                charge_pro += mol.molecule.charge
                products.append(mol.molecule)

        charge_rct = sum([r.charge for r in reaction["reactants"]])
        charge_pro = sum([p.charge for p in reaction["products"]])

        if charge_rct != charge_pro:
            print("Reactants and products do not have balanced charge! Skipping reaction {} in reactions".format(rr))
            continue

        if not compositions_equal(reactants, products):
            print("Reactants and products are not balanced! Skipping reaction {} in reactions".format(rr))
            continue

        job_name = job_name_prefix + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + str(jobs_made)
        job = AutoTSJob(reaction["reactants"], reaction["products"],
                        (base_dir / job_name),
                        schrodinger_dir=schrodinger_dir,
                        job_name=job_name, num_cores=num_cores, host=host,
                        save_scratch=save_scratch, input_params=input_params)
        jobs.append(job)
        jobs_made += 1

    # Prepare all jobs
    for job in jobs:
        job.setup_calculation()

    return_point = Path.cwd()

    for job in jobs:
        try:
            job.run(command_line_args=command_line_args)
        except RuntimeError:
            print("Failed to run job {}".format(job.job_name))

        os.chdir(return_point.as_posix())


def launch_jobs_from_queue(database: CatDB,
                           base_dir: Union[str, Path],
                           num_launches: int = 1,
                           by_priority: bool = False,
                           randomize: bool = False,
                           query: Optional[Dict] = None,
                           schrodinger_dir: Optional[str] = "SCHRODINGER",
                           num_cores: Optional[int] = 40,
                           host: Optional[str] = None,
                           save_scratch: Optional[bool] = False,
                           command_line_args: Optional[Dict] = None):
    """
    Use a CatDB.queue_collection to automatically submit jobs to a job manager.

    Args:
        database (CatDB): Database to query and identify calculations to run
        base_dir (str): Base directory in which to launch the calculations
        num_launches (int): Number of jobs to launch
        by_priority (bool): If True (default False), then prioritize jobs with
            higher priority, and do not run jobs with no priority assigned
        query (dict): MongoDB query dictionary.
        schrodinger_dir (str): A path to the Schrodinger Suite of software.
            This is used to call AutoTS and other utilities. By default,
            this is "$SCHRODINGER", which should be an environment variable
            set at the time of installation.
        num_cores (int): How many cores should the program be parallelized
            over (default 40). When multiple subjobs need to be run
            simultaneously, AutoTS will distribute these cores automatically
            between subjobs
        host (str): Which host should the calculation be run on? By default,
            this is "localhost", which should generally mean that the
            calculation is run on the current node without using a queueing
            system
        save_scratch (bool): If True (default False), save a *.zip file
            containing the contents of the calculation scratch directory
        command_line_args (dict): A dictionary of flag: value pairs to be
            provided to the autots command-line interface. Ex:
            {"WAIT": None,
             "subdir": None,
             "nsubjobs": 5}
             will be interpreted as "-WAIT -subdir -nsubjobs 5"
    """

    if isinstance(base_dir, str):
        base_dir = Path(base_dir)

    # If query is None, just grab all possible calculations
    queue_collection = database.database[database.queue_collection]

    if query is None:
        initial_query = [e for e in queue_collection.find({"state": "READY"},
                                                          {"_id": 0,
                                                           "rxnid": 1,
                                                           "priority": 1})]
    else:
        if "state" not in query:
            query["state"] = "READY"
        initial_query = [e for e in queue_collection.find(query,
                                                          {"_id": 0,
                                                           "rxnid": 1,
                                                           "priority": 1})]

    if by_priority:
        initial_query = sorted([i for i in initial_query if i["priority"] is not None],
                               key=lambda x: x["priority"],
                               reverse=True)
    if randomize:
        random.shuffle(initial_query)

    if num_launches > len(initial_query):
        raise ValueError("num_launches too high! Only {} jobs available".format(len(initial_query)))

    rxn_ids = [r["rxnid"] for r in initial_query[0:num_launches]]
    to_calculate = list(queue_collection.find({"rxnid": {"$in": rxn_ids}}))

    for calc in to_calculate:
        time_now = datetime.datetime.now(datetime.timezone.utc)
        reactants = [MoleculeGraph.from_dict(r) for r in calc["reactants"]]
        products = [MoleculeGraph.from_dict(p) for p in calc["products"]]
        timestamp = time_now.strftime("%Y%m%d_%H%M%S_%f")
        name = "_".join(["launcher", str(calc["rxnid"]), timestamp])
        job = AutoTSJob(reactants, products,
                        base_dir / name,
                        schrodinger_dir=schrodinger_dir,
                        num_cores=num_cores,
                        host=host,
                        save_scratch=save_scratch,
                        input_params=calc["input"])

        try:
            job.setup_calculation()

            queue_collection.update_one({"rxnid": calc["rxnid"]}, {"$set": {"state": "SUBMITTED",
                                                                            "updated_on": time_now}})

            calc_dict = {"rxnid": calc.get("rxnid"), "name": calc.get("name"),
                         "charge": calc.get("charge"),
                         "spin_multiplicity": calc.get("spin_multiplicity"),
                         "nelectrons": calc.get("nelectrons"),
                         "priority": calc.get("priority"),
                         "input": calc.get("input"),
                         "reactants": calc.get("reactants"),
                         "products": calc.get("products"),
                         "reaction_graph": calc.get("reaction_graph"),
                         "molgraph_pro": calc.get("molgraph_pro"),
                         "molgraph_rct": calc.get("molgraph_rct"),
                         "tags": calc.get("tags"),
                         "additional_data": calc.get("additional_data")}

            dumpfn(calc_dict, (base_dir / name / "calc.json").as_posix(), indent=2)
            job.run(command_line_args=command_line_args)

        except UnsupportedReaction:
            queue_collection.update_one({"rxnid": calc["rxnid"]},
                                        {"$set": {"state": "UNSUPPORTED",
                                                  "updated_on": time_now}})
