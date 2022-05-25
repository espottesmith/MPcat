# coding: utf-8

import os
import subprocess
from typing import Optional, List, Dict, Union
import datetime
from pathlib import Path
import random

from monty.serialization import dumpfn

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph

from schrodinger.application.jaguar.autots_exceptions import UnsupportedReaction

from mpcat.aggregate.database import CatDB
from mpcat.apprehend.jaguar_input import (JagSet,
                                          OptSet,
                                          TSOptSet,
                                          FreqSet,
                                          ScanSet,
                                          IRCSet,
                                          ElectronTransferSet)
from mpcat.apprehend.autots_input import TSSet
from mpcat.apprehend.jaguar_output import JagOutput
from mpcat.apprehend.autots_output import TSOutput
from mpcat.utils.generate import mol_to_mol_graph
from mpcat.utils.types import JaguarJobType, job_type_mapping


#TODO: Fix documentation in this file
#TODO: Make a Job superclass; significant overlap between JaguarJob and AutoTSJob


class JaguarJob:
    """
    A helper class to prepare and execute Jaguar calculations.
    """

    def __init__(self,
                 molecule: Union[Molecule, MoleculeGraph],
                 job_type: Union[str, JaguarJobType],
                 path: Union[str, Path],
                 schrodinger_dir: Optional[Union[str, Path]] = "SCHRODINGER",
                 job_name: Optional[str] = None,
                 num_cores: Optional[int] = 40,
                 host: Optional[str] = None,
                 save_scratch: Optional[bool] = False,
                 input_params: Optional[Dict] = None):

        """

        :param molecule:
        :param job_type:
        :param path:
        :param schrodinger_dir:
        :param job_name:
        :param num_cores:
        :param host:
        :param save_scratch:
        :param input_params:
        """

        self.molecule = mol_to_mol_graph(molecule)

        if isinstance(job_type, str):
            if job_type.lower() not in job_type_mapping:
                raise ValueError("Job type {} unknown!".format(job_type))

            self.job_type = job_type_mapping[job_type.lower()]
        else:
            self.job_type = job_type

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

        job_name = self.job_name or "jaguar"

        if self.job_type == JaguarJobType.SP:
            calc_set = JagSet(self.molecule, name=job_name, **self.input_params)
        elif self.job_type == JaguarJobType.OPT:
            calc_set = OptSet(self.molecule, name=job_name, **self.input_params)
        elif self.job_type == JaguarJobType.TS:
            calc_set = TSOptSet(self.molecule, name=job_name, **self.input_params)
        elif self.job_type == JaguarJobType.FREQ:
            calc_set = FreqSet(self.molecule, name=job_name, **self.input_params)
        elif self.job_type == JaguarJobType.SCAN:
            calc_set = ScanSet(self.molecule, name=job_name, **self.input_params)
        elif self.job_type == JaguarJobType.IRC:
            calc_set = IRCSet(self.molecule, name=job_name, **self.input_params)
        elif self.job_type == JaguarJobType.ET:
            calc_set = ElectronTransferSet(self.molecule, name=job_name, **self.input_params)
        else:
            raise NotImplementedError("No calculation set available for given calculation type!")

        calc_set.write(self.path)

    def run(self, command_line_args: Optional[Dict] = None):
        """
        Execute the Jaguar job

        Args:
            command_line_args (dict): A dictionary of flag: value pairs to be
            provided to the autots command-line interface. Ex:
            {"WAIT": None,
            "max_threads": 20}
             will be interpreted as "-WAIT -max_threads 20"

        Returns:
            None
        """

        job_name = self.job_name or "jaguar"

        cur_dir = Path.cwd()
        os.chdir(self.path.as_posix())

        command = [(self.schrodinger_dir / "jaguar").as_posix(), "run",
                   "-PARALLEL", str(self.num_cores)]

        if self.job_name is not None:
            command.append("-jobname")
            command.append(str(self.job_name))

        if self.host is not None:
            command.append("-HOST")
            command.append(self.host)

        if self.save_scratch:
            command.append("-SAVE")

        if command_line_args is not None:
            for key, value in command_line_args.items():
                command.append("-" + str(key))

                if value is not None:
                    command.append(str(value))

        command.append("{}.in".format(job_name))

        process = subprocess.run(command)

        os.chdir(cur_dir.as_posix())

        if process.returncode != 0:
            raise RuntimeError("Job launch failed!")


class AutoTSJob:
    """
    A helper class to prepare and execute AutoTS calculations.
    """

    def __init__(self,
                 reactants: List[Union[Molecule, MoleculeGraph]],
                 products: List[Union[Molecule, MoleculeGraph]],
                 path: Union[str, Path],
                 schrodinger_dir: Optional[Union[str, Path]] = "SCHRODINGER",
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
                to TSSet
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

        calc_set = TSSet(self.reactants, self.products,
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
                   "-PARALLEL", str(self.num_cores)]

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


def launch_jaguar_from_queue(database: CatDB,
                             base_dir: Union[str, Path],
                             num_launches: int = 1,
                             by_priority: bool = False,
                             randomize: bool = False,
                             job_name: bool = False,
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
    queue_collection = database.database[database.jaguar_queue_collection]

    if query is None:
        query = {"state": "READY"}
    else:
        if "state" not in query:
            query["state"] = "READY"

    initial_query = [e for e in queue_collection.find(query,
                                                      {"_id": 0,
                                                       "calcid": 1,
                                                       "priority": 1})]

    if by_priority:
        initial_query = sorted([i for i in initial_query if i["priority"] is not None],
                               key=lambda x: x["priority"],
                               reverse=True)
    if randomize:
        random.shuffle(initial_query)

    if num_launches > len(initial_query):
        raise ValueError("num_launches too high! Only {} jobs available".format(len(initial_query)))

    calc_ids = [r["calcid"] for r in initial_query[0:num_launches]]
    to_calculate = list(queue_collection.find({"calcid": {"$in": calc_ids}}))

    for calc in to_calculate:
        time_now = datetime.datetime.now(datetime.timezone.utc)
        molecule = MoleculeGraph.from_dict(calc["molecule"])
        job_type = job_type_mapping[calc["job_type"].lower()]
        timestamp = time_now.strftime("%Y%m%d_%H%M%S_%f")
        name = "_".join(["launcher", str(calc["calcid"]), timestamp])
        job = JaguarJob(molecule,
                        job_type,
                        base_dir / name,
                        job_name=str(calc["calcid"]) if job_name else None,
                        schrodinger_dir=schrodinger_dir,
                        num_cores=num_cores,
                        host=host,
                        save_scratch=save_scratch,
                        input_params=calc["input"])

        try:
            job.setup_calculation()

            queue_collection.update_one({"calcid": calc["calcid"]}, {"$set": {"state": "SUBMITTED",
                                                                              "updated_on": time_now}})

            calc_dict = {"calcid": calc.get("calcid"),
                         "name": calc.get("name"),
                         "job_type": calc.get("job_type"),
                         "molecule": calc.get("molecule"),
                         "charge": calc.get("charge"),
                         "spin_multiplicity": calc.get("spin_multiplicity"),
                         "nelectrons": calc.get("nelectrons"),
                         "priority": calc.get("priority"),
                         "input": calc.get("input"),
                         "formula_alphabetical": calc.get("formula_alphabetical"),
                         "tags": calc.get("tags"),
                         "additional_data": calc.get("additional_data")}

            dumpfn(calc_dict, (base_dir / name / "calc.json").as_posix(), indent=2)
            job.run(command_line_args=command_line_args)

        except UnsupportedReaction:
            queue_collection.update_one({"rxnid": calc["rxnid"]},
                                        {"$set": {"state": "UNSUPPORTED",
                                                  "updated_on": time_now}})


def launch_autots_from_queue(database: CatDB,
                                  base_dir: Union[str, Path],
                                  num_launches: int = 1,
                                  by_priority: bool = False,
                                  randomize: bool = False,
                                  job_name: bool = False,
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
    queue_collection = database.database[database.autots_queue_collection]

    if query is None:
        query = {"state": "READY"}
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
                        job_name=str(calc["rxnid"]) if job_name else None,
                        schrodinger_dir=schrodinger_dir,
                        num_cores=num_cores,
                        host=host,
                        save_scratch=save_scratch,
                        input_params=calc["input"])

        try:
            job.setup_calculation()

            queue_collection.update_one({"rxnid": calc["rxnid"]}, {"$set": {"state": "SUBMITTED",
                                                                            "updated_on": time_now}})

            calc_dict = {"rxnid": calc.get("rxnid"),
                         "name": calc.get("name"),
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


def recover_job_autots(schrodinger_dir: Optional[Union[str, Path]] = "SCHRODINGER",
                       check_if_needed: bool = False,
                       database: Optional[CatDB] = None,
                       rxnid: Optional[int] = None,
                       directory: Optional[Path] = None,
                       recovery_name: Optional[str] = None,
                       output_name: Optional[str] = None):
    """
    Recover an unfinished AutoTS calculation.

    :param schrodinger_dir: A path to the Schrodinger Suite of software.
        This is used to call AutoTS and other utilities. By default,
        this is "$SCHRODINGER", which should be an environment variable
        set at the time of installation.
    :param check_if_needed: If True (default False), then parse output file to
        check if calculation has finished before attempting to recover the calculation.
    :param database: Connection to MongoDB database. Only needed if attempting to recover
        an unfinished calculation from DB. Default None.
    :param rxnid: Reaction ID of the calculation to be recovered. Only needed if attempting
        to recover an unfinished calculation from DB. Default None.
    :param directory: Path to the calculation directory of interest. Only needed if
        attempting to recover from directory without using DB. Default None.
    :param recovery_name: Name of recover file to use to restart calculation. If None (default),
        will look for file called "AutoTS.T9XnCsLi.recover".
    :param output_name: Name of output file to use to (optionally) check for calculation completion.
        If None (default), will look for file called "AutoTS.T9XnCsLi.out".
    :return:
    """

    if schrodinger_dir == "SCHRODINGER":
        schrodinger_dir = Path(os.environ["SCHRODINGER"])
    else:
        if isinstance(schrodinger_dir, str):
            schrodinger_dir = Path(schrodinger_dir)
        else:
            schrodinger_dir = schrodinger_dir

    # Make sure some necessary inputs are provided
    if (database is None or rxnid is None):
        if directory is None:
            raise ValueError("Either provide database and rxnid, or provide directory!")
        else:
            from_db = False
    else:
        from_db = True

    if from_db:
        calc = database.database[database.autots_data_collection].find_one({"rxnid": rxnid})
        path = Path(calc["path"])
    else:
        path = Path(directory)

    # Check if path exists
    if not path.is_dir():
        raise FileNotFoundError("Path does not exist!")

    files = [e.name for e in path.iterdir()]

    # Check if calculation has finished
    if check_if_needed:
        if output_name is not None:
            filestring = (path / output_name).as_posix()
        else:
            if "AutoTS.T9XnCsLi.out" in files:
                filestring = (path / "AutoTS.T9XnCsLi.out").as_posix()
            else:
                raise FileNotFoundError("Output file could not be found.")

        if not Path(filestring).exists():
            raise FileNotFoundError("Output file could not be found.")

        # If calculation has completed, return without doing anything
        out = TSOutput(filestring)
        if out.data["walltime"] is not None:
            print("Calculation completed! No recovery necessary")
            return

    if recovery_name is not None:
        filestring = (path / recovery_name).as_posix()
    else:
        if "AutoTS.T9XnCsLi.recover" in files:
            filestring = (path / "AutoTS.T9XnCsLi.recover").as_posix()
        else:
            raise FileNotFoundError("Recovery file could not be found.")

    if not Path(filestring).exists():
        raise FileNotFoundError("Recovery file could not be found.")

    cur_dir = Path.cwd()
    os.chdir(path.as_posix())

    command = [(schrodinger_dir / "jaguar").as_posix(), "run", filestring]

    process = subprocess.run(command)

    os.chdir(cur_dir.as_posix())

    if process.returncode != 0:
        raise RuntimeError("Job launch failed!")


def recover_job_jaguar(schrodinger_dir: Optional[Union[str, Path]] = "SCHRODINGER",
                       check_if_needed: bool = False,
                       database: Optional[CatDB] = None,
                       calcid: Optional[int] = None,
                       directory: Optional[Path] = None,
                       recovery_name: Optional[str] = None,
                       output_name: Optional[str] = None):
    """
    Recover an unfinished Jaguar calculation.

    :param schrodinger_dir: A path to the Schrodinger Suite of software.
        This is used to call AutoTS and other utilities. By default,
        this is "$SCHRODINGER", which should be an environment variable
        set at the time of installation.
    :param check_if_needed: If True (default False), then parse output file to
        check if calcultion has finished before attempting to recover the calculation.
    :param database: Connection to MongoDB database. Only needed if attempting to recover
        an unfinished calculation from DB. Default None.
    :param calcid: Calculation ID of the calculation to be recovered. Only needed if
        attempting to recover an unfinished calculation from DB. Default None.
    :param directory: Path to the calculation directory of interest. Only needed if
        attempting to recover from directory without using DB. Default None.
    :param recovery_name: Name of recover file to use to restart calculation. If None (default),
        will look for file called "AutoTS.T9XnCsLi.recover".
    :param output_name: Name of output file to use to (optionally) check for calculation completion.
        If None (default), will look for file called "AutoTS.T9XnCsLi.out".
    :return:
    """

    if schrodinger_dir == "SCHRODINGER":
        schrodinger_dir = Path(os.environ["SCHRODINGER"])
    else:
        if isinstance(schrodinger_dir, str):
            schrodinger_dir = Path(schrodinger_dir)
        else:
            schrodinger_dir = schrodinger_dir

    # Make sure some necessary inputs are provided
    if (database is None or calcid is None):
        if directory is None:
            raise ValueError("Either provide database and rxnid, or provide directory!")
        else:
            from_db = False
    else:
        from_db = True

    if from_db:
        calc = database.database[database.autots_data_collection].find_one({"calcid": calcid})
        path = Path(calc["path"])
    else:
        path = Path(directory)

    # Check if path exists
    if not path.is_dir():
        raise FileNotFoundError("Path does not exist!")

    files = [e.name for e in path.iterdir()]

    # Check if calculation has finished
    if check_if_needed:
        if output_name is not None:
            filestring = (path / output_name).as_posix()
        else:
            if "jaguar.out" in files:
                filestring = (path / "jaguar.out").as_posix()
            else:
                raise FileNotFoundError("Output file could not be found.")

        if not Path(filestring).exists():
            raise FileNotFoundError("Output file could not be found.")

        # If calculation has completed, return without doing anything
        out = JagOutput(filestring)
        if out.data["walltime"] is not None:
            print("Calculation completed! No recovery necessary")
            return

    if recovery_name is not None:
        filestring = (path / recovery_name).as_posix()
    else:
        if "jaguar.recover" in files:
            filestring = (path / "jaguar.recover").as_posix()
        else:
            raise FileNotFoundError("Recovery file could not be found.")

    if not Path(filestring).exists():
        raise FileNotFoundError("Recovery file could not be found.")

    cur_dir = Path.cwd()
    os.chdir(path.as_posix())

    command = [(schrodinger_dir / "jaguar").as_posix(), "run", filestring]

    process = subprocess.run(command)

    os.chdir(cur_dir.as_posix())

    if process.returncode != 0:
        raise RuntimeError("Job launch failed!")