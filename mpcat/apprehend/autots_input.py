# coding: utf-8

import os
from typing import List, Dict, Optional
import copy

from monty.json import MSONable
from pymatgen.core.structure import Molecule

from schrodinger.application.jaguar.reactiq_input import (ReactiqInput)

from mpcat.adapt.schrodinger_adapter import (molecule_to_maestro_file,
                                             maestro_file_to_molecule)


class AutoTSInput(MSONable):
    """
    An object representing an AutoTS input, including the actual workflow input
    file, complete with all AutoTS-specific and Jaguar $gen keyword variables,
    as well as the reactants and products of the reaction.

    Args:
        reactants (list of Molecule objects)
        products (list of Molecule objects)
        autots_variables (dict): Dictionary of AutoTS-specific keywords
            variables and values. For instance,
            {"units": "ev", "use_template": False}
        gen_variables (dict): Dictionary of Jaguar $gen section keyword
            variables and values. For instance,
            {"irder": 1, "maxitg": 200, "basis": "def2-tzvpd", "babel": "xyz"}
    """

    def __init__(self, reactants: List[Molecule], products: List[Molecule],
                 autots_variables: Dict, gen_variables: Dict):

        self.reactants = reactants
        self.products = products
        self.autots_variables = autots_variables
        self.gen_variables = gen_variables

        self.autots_variables["charge"] = sum([m.charge for m in self.reactants])

        nelectrons = sum([m._nelectrons for m in self.reactants])
        self.autots_variables["multiplicity"] = 1 if nelectrons % 2 == 0 else 2

    def write(self, filename: str, write_molecules: Optional[bool] = True,
              jobname: Optional[str] = None):
        """
        Write an AutoTS input file.

        Args:
            filename (str): Path to AutoTS workflow input file.
            write_molecules (bool): If True (default), write all reactant and
                product molecules as *.mae files.
            jobname (str): If not None (default), this will be the name of the
                AutoTS workflow

        Returns:
            None
        """

        base_dir = os.path.dirname(filename)

        self.autots_variables["reactant"] = [os.path.join(base_dir, "rct_{}.mae".format(rr))
                                             for rr in range(len(self.reactants))]
        self.autots_variables["product"] = [os.path.join(base_dir, "pro_{}.mae".format(pp))
                                            for pp in range(len(self.products))]

        if write_molecules:
            for rr, reactant in enumerate(self.reactants):
                molecule_to_maestro_file(reactant, os.path.join(base_dir,
                                                                "rct_{}.mae".format(rr)))
            for pp, product in enumerate(self.products):
                molecule_to_maestro_file(product, os.path.join(base_dir,
                                                               "pro_{}.mae".format(pp)))

        input_file = ReactiqInput(keywords=self.autots_variables,
                                  jaguar_keywords=self.gen_variables,
                                  jobname=jobname)
        input_file.save(filename)

    @classmethod
    def from_file(cls, filename: str, read_molecules: bool = True):
        """
        Parse an AutoTS workflow input file and store its data into an
            AutoTSInput object.

        Args:
            filename (str): Path to an existing AutoTS input file.
            read_molecules (bool): If True, also process all reactant and
                product files

        Returns:
            autots_input: AutoTSInput object
        """

        input_file = ReactiqInput()
        input_file.read(filename)

        keywords = input_file._keywords.keys()
        autots_variables = dict()

        for keyword in keywords:
            autots_variables[keyword] = input_file[keyword]

        gen_variables = copy.deepcopy(input_file._jaguar_user_keys)

        if read_molecules:
            base_dir = os.path.dirname(filename)

            reactants = input_file.getValue("reactant")
            products = input_file.getValue("product")

            rct_mols = [maestro_file_to_molecule(os.path.join(base_dir, r))[0]
                        for r in reactants]
            pro_mols = [maestro_file_to_molecule(os.path.join(base_dir, p))[0]
                        for p in products]
        else:
            rct_mols = list()
            pro_mols = list()

        return cls(rct_mols, pro_mols, autots_variables, gen_variables)


class AutoTSSet(AutoTSInput):
    """
    Build an AutoTSInput given various input parameters.
    """

    def __init__(self,
                 reactants,
                 products,
                 basis_set="def2-tzvpd",
                 dft_rung=4,
                 pcm_dielectric=None,
                 max_scf_cycles=400,
                 geom_opt_max_cycles=250,
                 overwrite_inputs_autots=None,
                 overwrite_inputs_gen=None):
        """
        Args:
            reactants (list of Molecule objects):
            products (list of Molecule objects):
            basis_set (str):
            dft_rung (int):
            pcm_dielectric (float):
            max_scf_cycles (int):
            geom_opt_max_cycles (int):
            overwrite_inputs_autots (dict): Dictionary to overwrite default
                AutoTS-specific parameters
            overwrite_inputs_gen (dict): Dictionary to overwrite default
                Jaguar gen parameters
        """

        autots_variables = {"eliminate_multiple_frequencies": True,
                            "free_energy": True,
                            "require_irc_success": True,
                            "ts_vet_max_freq": -40.0,
                            "units": "ev",
                            "use_alternate": True}

        if dft_rung == 1:
            dftname = "b3lyp"
        elif dft_rung == 2:
            dftname = "cam-b3lyp"
        elif dft_rung == 3:
            dftname = "cam-b3lyp-d3"
        elif dft_rung == 4:
            dftname = "wb97x-d"
        else:
            raise ValueError("Invalid dft_rung provided!")

        gen_variables = {"dftname": dftname,
                         "basis": basis_set,
                         "babel": "xyz",
                         "ip472": 2,  # Output all steps of geometry optimization in *.mae
                         "ip172": 2,  # Print RESP file
                         "ip175": 2,  # Print XYZ files
                         "ifreq": 1,  # Frequency calculation
                         "irder": 1,  # IR vibrational modes calculated
                         "nmder": 2,  # Numerical second derivatives
                         "nogas": 2,  # Skip gas-phase optimization, if PCM is used
                         "maxitg": geom_opt_max_cycles,  # Maximum number of geometry optimization iterations
                         "intopt_switch": 0,  # Do not switch from internal to Cartesian coordinates
                         "optcoord_update": 0,  # Do not run checks to change coordinate system
                         "props_each_step": 1,  # Calculate properties at each optimization step
                         # "iaccg": 5  # Tight convergence criteria for optimization
                         "mulken": 1,  # Calculate Mulliken properties by atom
                         "maxit": max_scf_cycles,  # Maximum number of SCF iterations
                         "iacc": 2,  # Use "accurate" SCF convergence criteria
                         # "noauto": 3  # All calculations done on fine grid
                         "isymm": 0,  # Do not use symmetry
                         "espunit": 6  # Electrostatic potential in units of eV
                         }

        if pcm_dielectric is not None:
            gen_variables["isolv"] = 7
            gen_variables["epsout"] = pcm_dielectric
            gen_variables["pcm_model"] = "cosmo"

        if overwrite_inputs_autots is not None:
            for key, value in overwrite_inputs_autots.items():
                autots_variables[key] = value

        if overwrite_inputs_gen is not None:
            for key, value in overwrite_inputs_gen.items():
                gen_variables[key] = value

        super().__init__(reactants, products, autots_variables, gen_variables)