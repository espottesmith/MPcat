# coding: utf-8

import os
import subprocess
import shutil
from typing import Optional, List, Dict, Type
import datetime

from pymatgen.core.structure import Molecule
# from pymatgen.reaction_network.reaction_network import (ReactionPath,
#                                                         ReactionNetwork)

from mpcat.apprehend.autots_input import AutoTSSet
from mpcat.utils.comparison import compositions_equal


class AutoTSJob:
    """
    A helper class to prepare and execute AutoTS calculations.
    """

    def __init__(self,
                 reactants: List[Molecule],
                 products: List[Molecule],
                 path: str,
                 schrodinger_dir: Optional[str] = "$SCHRODINGER",
                 job_name: Optional[str] = None,
                 num_cores: Optional[int] = 40,
                 host: Optional[str] = "localhost",
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

        self.reactants = reactants
        self.products = products
        self.path = path

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

        if not os.path.exists(self.path):
            os.mkdir(self.path)

        calc_set = AutoTSSet(self.reactants, self.products,
                             **self.input_params)

        calc_set.write(os.path.join(self.path, "autots.in"), write_molecules=True,
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

        os.chdir(self.path)

        command = [self.schrodinger_dir + "/autots",
                   "-jobname", self.job_name,
                   "-PARALLEL", str(self.num_cores),
                   "-HOST", self.host, "-use_one_node"]

        if self.save_scratch:
            command.append("-SAVE")

        if command_line_args is not None:
            for key, value in command_line_args.items():
                command.append("-" + str(key))

                if value is not None:
                    command.append(str(value))

        command.append("autots.in")

        process = subprocess.run(command,
                                 capture_output=True,
                                 text=True)
        if process.returncode != 0:
            raise RuntimeError("Job launch failed!")


def launch_mass_jobs(reactions: List[Dict[str, List[Molecule]]],
                     base_dir: str,
                     schrodinger_dir: Optional[str] = "$SCHRODINGER",
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

    jobs = list()

    jobs_made = 0

    # Make an AutoTSJob object for each reaction
    for rr, reaction in enumerate(reactions):

        # Verify that reaction is balanced in terms of charge and species

        charge_rct = sum([r.charge for r in reaction["reactants"]])
        charge_pro = sum([p.charge for p in reaction["products"]])

        if charge_rct != charge_pro:
            print("Reactants and products do not have balanced charge! Skipping reaction {} in reactions".format(rr))
            continue

        if not compositions_equal(reaction["reactants"], reaction["products"]):
            print("Reactants and products are not balanced! Skipping reaction {} in reactions".format(rr))
            continue

        job_name = job_name_prefix + datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + str(jobs_made)
        job = AutoTSJob(reaction["reactants"], reaction["products"],
                        os.path.join(base_dir, job_name),
                        schrodinger_dir=schrodinger_dir,
                        job_name=job_name, num_cores=num_cores, host=host,
                        save_scratch=save_scratch, input_params=input_params)
        jobs.append(job)
        jobs_made += 1

    # Prepare all jobs
    for job in jobs:
        job.setup_calculation()

    return_point = os.getcwd()

    for job in jobs:
        try:
            job.run(command_line_args=command_line_args)
        except RuntimeError:
            print("Failed to run job {}".format(job.job_name))

        os.chdir(return_point)


# def launch_reaction_path(reaction_network: ReactionNetwork,
#                          reaction_path: ReactionPath,
#                          base_dir: str,
#                          schrodinger_dir: Optional[str] = "$SCHRODINGER",
#                          job_name_prefix: Optional[str] = None,
#                          num_cores: Optional[int] = 40,
#                          host: Optional[str] = "localhost",
#                          save_scratch: Optional[bool] = False,
#                          input_params: Optional[Dict] = None,
#                          command_line_args: Optional[Dict] = None):
#     """
#     For all reactions along a reaction path, prepare and launch AutoTS
#         calculations for them.
#
#     Args:
#         reaction_network (ReactionNetwork): A network containing the reactions
#             in the ReactionPath
#         reaction_path (ReactionPath): The path to be studied
#         base_dir (str): Root directory where all calculation directories should
#             be made
#         schrodinger_dir (str): A path to the Schrodinger Suite of software.
#             This is used to call AutoTS and other utilities. By default,
#             this is "$SCHRODINGER", which should be an environment variable
#             set at the time of installation.
#         job_name_prefix (str): All jobs in this set of reactions will be given a
#             unique name, but this prefix will be prepended to all calculations
#             in this set.
#         num_cores (int): How many cores should the program be parallelized
#             over (default 40). When multiple subjobs need to be run
#             simultaneously, AutoTS will distribute these cores automatically
#             between subjobs
#         host (str): Which host should the calculation be run on? By default,
#             this is "localhost", which should generally mean that the
#             calculation is run on the current node without using a queueing
#             system
#         save_scratch (bool): If True (default False), save a *.zip file
#             containing the contents of the calculation scratch directory
#         input_params (dict): Keywords and associated values to be provided
#             to AutoTSSet
#         command_line_args (dict): A dictionary of flag: value pairs to be
#             provided to the autots command-line interface. Ex:
#             {"WAIT": None,
#              "subdir": None,
#              "nsubjobs": 5}
#              will be interpreted as "-WAIT -subdir -nsubjobs 5"
#
#     Returns:
#         None
#     """
#
#     reactions = list()
#
#     for node in reaction_path:
#         # We only care about reaction nodes
#         if "," in node:
#             sides = node.split(',')
#             rct_ids = [int(r.replace("PR_", "")) for r in sides[0].split("+")]
#             pro_ids = [int(p) for p in sides[1].split("+")]
#
#             rcts = [reaction_network.entries_list[e].molecule for e in rct_ids]
#             pros = [reaction_network.entries_list[e].molecule for e in pro_ids]
#             reactions.append({"reactants": rcts, "products": pros})
#
#     launch_mass_jobs(reactions, base_dir=base_dir,
#                      schrodinger_dir=schrodinger_dir,
#                      job_name_prefix=job_name_prefix,
#                      num_cores=num_cores,
#                      host=host, save_scratch=save_scratch,
#                      input_params=input_params,
#                      command_line_args=command_line_args)