# coding: utf-8

from typing import List, Dict, Optional, Union
import copy
from pathlib import Path

from monty.json import MSONable
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph

from schrodinger.application.jaguar.reactiq_input import (ReactiqInput)

from mpcat.adapt.schrodinger_adapter import (mol_graph_to_maestro_file,
                                             maestro_file_to_molecule)
from mpcat.utils.generate import mol_to_mol_graph


class JaguarInput(MSONable):
    """
    An object representing a Jaguar input, including the molecule information
    and any calculation parameters.
    """
    pass


class AutoTSInput(MSONable):
    """
    An object representing an AutoTS input, including the actual workflow input
    file, complete with all AutoTS-specific and Jaguar $gen keyword variables,
    as well as the reactants and products of the reaction.

    Args:
        reactants (list of Molecule or MoleculeGraph objects)
        products (list of Molecule or MoleculeGraph objects)
        autots_variables (dict): Dictionary of AutoTS-specific keywords
            variables and values. For instance,
            {"units": "ev", "use_template": False}
        gen_variables (dict): Dictionary of Jaguar $gen section keyword
            variables and values. For instance,
            {"irder": 1, "maxitg": 200, "basis": "def2-tzvpd", "babel": "xyz"}
        spin_multiplicity (int): overall system spin multiplicity
    """

    def __init__(self, reactants: List[Union[Molecule, MoleculeGraph]], products: List[Molecule],
                 autots_variables: Dict, gen_variables: Dict,
                 spin_multiplicity: Optional[int] = None):

        self.reactants = list()
        for reactant in reactants:
            self.reactants.append(mol_to_mol_graph(reactant))

        self.products = list()
        for product in products:
            self.products.append(mol_to_mol_graph(product))

        self.autots_variables = autots_variables
        self.gen_variables = gen_variables

        rct_sum = int(sum([m.molecule.charge for m in self.reactants]))
        pro_sum = int(sum([m.molecule.charge for m in self.products]))

        if rct_sum != pro_sum:
            raise ValueError("Reactant and product charges do not match!")

        self.autots_variables["charge"] = rct_sum

        nelectrons = int(sum([m.molecule._nelectrons for m in self.reactants]))
        if spin_multiplicity is None:
            if len(self.reactants) == len(self.products) == 1:
                self.autots_variables["multiplicity"] = self.reactants[0].molecule.spin_multiplicity
            else:
                self.autots_variables["multiplicity"] = 1 if nelectrons % 2 == 0 else 2
        else:
            if (spin_multiplicity % 2) == (nelectrons % 2):
                raise ValueError("Invalid spin multiplicity: {} with {} electrons!".format(spin_multiplicity, nelectrons))
            else:
                self.autots_variables["multiplicity"] = spin_multiplicity


    def write(self, filename: Union[str, Path], write_molecules: Optional[bool] = True,
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

        if isinstance(filename, str):
            base_dir = Path(filename).parent
        else:
            base_dir = filename.parent

        self.autots_variables["reactant"] = [(base_dir / "rct_{}.mae".format(rr)).as_posix()
                                             for rr in range(len(self.reactants))]
        self.autots_variables["product"] = [(base_dir / "pro_{}.mae".format(pp)).as_posix()
                                            for pp in range(len(self.products))]

        if write_molecules:
            for rr, reactant in enumerate(self.reactants):
                mol_graph_to_maestro_file(reactant,
                                          base_dir / "rct_{}.mae".format(rr))
            for pp, product in enumerate(self.products):
                mol_graph_to_maestro_file(product,
                                          base_dir / "pro_{}.mae".format(pp))

        input_file = ReactiqInput(keywords=self.autots_variables,
                                  jaguar_keywords=self.gen_variables,
                                  jobname=jobname)
        input_file.save(filename)

    @classmethod
    def from_file(cls, filename: Union[str, Path], read_molecules: bool = True):
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

        if isinstance(filename, Path):
            fn = filename.as_posix()
        else:
            fn = filename

        input_file = ReactiqInput()
        input_file.read(fn)

        keywords = input_file._keywords.keys()
        autots_variables = dict()

        for keyword in keywords:
            autots_variables[keyword] = input_file[keyword]

        gen_variables = copy.deepcopy(input_file._jaguar_user_keys)

        if read_molecules:
            base_dir = Path(filename).parent

            reactants = input_file.getValue("reactant")
            products = input_file.getValue("product")

            rct_mols = [maestro_file_to_molecule(base_dir / r)[0]
                        for r in reactants]
            pro_mols = [maestro_file_to_molecule(base_dir / p)[0]
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
                 spin_multiplicity=None,
                 basis_set="def2-tzvpd",
                 dft_rung=4,
                 pcm_dielectric=None,
                 max_scf_cycles=400,
                 geom_opt_max_cycles=250,
                 overwrite_inputs_autots=None,
                 overwrite_inputs_gen=None):
        """
        Args:
            reactants (list of Molecule or MoleculeGraph objects):
            products (list of Molecule or MoleculeGraph objects):
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
                            "units": "ev"}

        if dft_rung == 1:
            dftname = "hfs"
        elif dft_rung == 2:
            dftname = "b97-d3"
        elif dft_rung == 3:
            dftname = "m06-l"
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

        super().__init__(reactants, products, autots_variables, gen_variables,
                         spin_multiplicity=spin_multiplicity)