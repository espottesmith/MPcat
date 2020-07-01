import os
import shutil
from typing import Optional, List, Dict

from pymatgen.core.structure import Molecule
from pymatgen.reaction_network.reaction_network import (Reaction,
                                                        IntramolSingleBondChangeReaction,
                                                        IntermolecularReaction,
                                                        CoordinationBondChangeReaction,
                                                        ConcertedReaction,
                                                        ReactionPath,
                                                        ReactionNetwork)

from mpcat.apprehend.autots_input import AutoTSInput, AutoTSSet


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
        Prepare for calculation - prepare directory with all necessary input files.
        """

        if not os.path.exists(self.path):
            os.mkdir(self.path)

        calc_set = AutoTSSet(self.reactants, self.products,
                             **self.input_params)

        calc_set.write(os.path.join(self.path, "autots.in"), write_molecules=True,
                       jobname=self.job_name)

    def run(self, command_line_args: Optional[Dict] = None):
        """

        :return:
        """
        #TODO
        pass