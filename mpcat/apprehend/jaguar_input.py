# coding: utf-8

from typing import List, Dict, Optional, Union
import copy
from pathlib import Path
import os

import numpy as np

from monty.json import MSONable
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph

import schrodinger.infra.mm as mm
from schrodinger.application.jaguar.input import JaguarInput, read
from schrodinger.application.jaguar import validation

from mpcat.adapt.schrodinger_adapter import (molecule_to_schrodinger_struct,
                                             mol_graph_to_schrodinger_struct,
                                             schrodinger_struct_to_molecule)


def get_default_gen():
    """
    Returns reasonable default values for Jaguar DFT calculations.

    Args:
        None

    Returns;
        Dict
    """
    return {
        "babel": "xyz",
         "ip472": 2,  # Output all steps of geometry optimization in *.mae
         "ip175": 2,  # Print XYZ files
         "nogas": 2,  # Skip gas-phase optimization, if PCM is used
         "iacc": 2,  # Use "accurate" SCF convergence criteria
         "isymm": 0,  # Do not use symmetry
         # "iaccg": 5  # Tight convergence criteria for optimization
         # "noauto": 3  # All calculations done on fine grid
    }


def generate_gen(dft_rung: int,
                 basis_set: str,
                 pcm_settings: Optional[Dict] = None,
                 max_scf_cycles: int = 400,
                 max_geom_opt_cycles: Optional[int] = None,
                 overwrite_inputs_gen: Optional[Dict] = None):

    gen = get_default_gen()
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

    gen["dftname"] = dftname
    gen["basis"] = basis_set
    gen["maxit"] = max_scf_cycles

    if max_geom_opt_cycles is not None:
        gen["maxitg"] = max_geom_opt_cycles

    if pcm_settings is not None:
        gen["isolv"] = 7
        if pcm_settings["solvent"] == "other":
            gen["solvent"] = "other"
            gen["epsout"] = pcm_settings["dielectric"]
            gen["epsout_opt"] = pcm_settings["optical"]

            # Calculate the probe radius
            delta = 0.5
            rho = pcm_settings["density"]
            m = pcm_settings["molar_mass"]
            r = round(((3 * m / 6.023 * delta) / (4 * np.pi * rho) * 10) ** (1/3), 3)
            gen["radprb"] = r
        else:
            gen["solvent"] = pcm_settings["solvent"]
        gen["pcm_model"] = pcm_settings["model"]

    if overwrite_inputs_gen is not None:
        for key, value in overwrite_inputs_gen.items():
            gen[key] = value

    return gen

class JagInput(MSONable):
    """
    An object representing a Jaguar input, including the molecule information
    and any calculation parameters.

    This is a light wrapper around schrodinger.application.jaguar.input

    Args:
        molecule (Union[Molecule, MoleculeGraph]): molecule object that will
            be the subject of this calculation. Will be converted to a Schrodinger
            Structure object.
        gen_variables (Dict): Dictionary of Jaguar inputs
        name (str): Name for this Jaguar input. Default is "jaguar.in"
        jagin (Optional[JaguarInput]): Already-made JaguarInput. Used in
            other constructor methods like from_file
    """

    def __init__(self,
                 molecule: Union[Molecule, MoleculeGraph],
                 gen_variables: Dict,
                 name: str = "jaguar.in",
                 jagin: Optional[JaguarInput] = None
                 ):

        self.name = name

        # Construct Schrodinger Structure object
        self.mol = molecule
        if isinstance(molecule, Molecule):
            struct = molecule_to_schrodinger_struct(molecule)
        else:
            struct = mol_graph_to_schrodinger_struct(molecule)
        self.struct = struct

        if jagin:
            self.jagin = jagin
        else:
            self.jagin = JaguarInput(name=name, structure=self.struct)

        self.modify_gen({k: str(v) for k, v in gen_variables.items()})

    def modify_gen(self, new_gen: Dict):
        """
        Update the dictionary of gen key-value pairs for JaguarInput.

        Args:
            new_gen (Dict): Dictionary of Jaguar inputs

        Returns:
            None
        """
        # Validate gen keyword-value pairs
        for k, v in new_gen.items():
            validation.keyword_value_pair_is_valid(k, v)
        self.jagin.setValues(new_gen)

    def write(self, directory: Path):
        """
        Write Jaguar input file

        Args:
            directory (Path): Path to write the input file to

        Returns:
            None
        """

        # Move to target directory
        curdir = os.getcwd()
        os.chdir(directory.as_posix())

        self.jagin.save()

        # Move back
        os.chdir(curdir)

    @classmethod
    def from_file(cls, file: Path,
                  name: str = "jaguar.in",
                  overwrite_mol: Optional[Union[Molecule, MoleculeGraph]] = None,
                  overwrite_gen: Optional[Dict] = None):
        """
        Construct a JagInput from an existing Jaguar input file.

        Args:
            file (Path): Path to the existing Jaguar input file]
            name (str): Name for this Jaguar input. Default is "jaguar.in"
            overwrite_mol (Optional[Union[Molecule, MoleculeGraph]]): If not None
                (default), use this molecule instead of the molecule included in
                the existing Jaguar input file
            overwrite_gen (Optional[Dict]): If not None (default), provide
                additional/updated gen key-value pairs.
        """

        jagin = read(file.as_posix())

        if overwrite_mol is None:
            mol = schrodinger_struct_to_molecule(jagin.getStructure())
        else:
            mol = overwrite_mol

        if overwrite_gen is None:
            gen = dict()
        else:
            gen = overwrite_gen

        return cls(mol, gen, name=name, jagin=jagin)


class JagSet(JagInput):
    """
    JagInput (Jaguar input object) with some helpful defaults, aimed to improve
    ease of use.

    Args;
         molecule (Union[Molecule, MoleculeGraph]): molecule object that will
            be the subject of this calculation. Will be converted to a Schrodinger
            Structure object.
        name (str): Name for this Jaguar input. Default is "jaguar.in"
        dft_rung (int): Which type of DFT functional should be used. Values
            between 1-4 are accepted, with default functionals being chosen:
            1: "hfs"
            2: "b97-d3"
            3: "m06-l"
            4: "wb97x-d"
        basis_set (str): Basis set to be used. Default is "def2-svpd(-f)", a relatively
            small split-valence basis with polarization and diffuse functions
            included.
        pcm_settings (Dict): If not None (default), then the polarizable continuum model
            (PCM) will be used to construct an implicit solvation environment around
            the molecule of interest. pcm_settings requires at least two keys,
            solvent and model. If solvent is "other", then four additional keys are
            required, dielectric (solvent dielectric constant), optical (square of the
            solvent index of refraction), density (solvent density in g/cm^3) and
            molar_mass (solvent molar mass in g/mol).
        max_scf_cycles (int): Maximum SCF cycles to be allowed before the calculation
            fails. The detault is 400, but for simple calculations, a much lower number
            (perhaps 100) is reasonable.
        overwrite_inputs_gen (Dict): Dictionary of Jaguar inputs that should overwrite
            defaults

    """

    def __init__(self,
                 molecule: Union[Molecule, MoleculeGraph],
                 name: str = "jaguar.in",
                 dft_rung: int = 4,
                 basis_set: str = "def2-svpd(-f)",
                 pcm_settings: Optional[Dict] = None,
                 max_scf_cycles: int = 400,
                 overwrite_inputs_gen: Optional[Dict] = None):

        if isinstance(molecule, Molecule):
            mol = molecule
        else:
            mol = molecule.molecule

        gen = generate_gen(dft_rung, basis_set,
                           pcm_settings=pcm_settings,
                           max_scf_cycles=max_scf_cycles,
                           overwrite_inputs_gen=overwrite_inputs_gen)
        gen["molchg"] = mol.charge
        gen["multip"] = mol.spin_multiplicity

        super().__init__(molecule, gen, name=name)


class OptSet(JagSet):
    """
    JagSet object for use with geometry optimization calculations

    Args;
         molecule (Union[Molecule, MoleculeGraph]): molecule object that will
            be the subject of this calculation. Will be converted to a Schrodinger
            Structure object.
        name (str): Name for this Jaguar input. Default is "jaguar.in"
        dft_rung (int): Which type of DFT functional should be used. Values
            between 1-4 are accepted, with default functionals being chosen:
            1: "hfs"
            2: "b97-d3"
            3: "m06-l"
            4: "wb97x-d"
        basis_set (str): Basis set to be used. Default is "def2-svpd(-f)", a relatively
            small split-valence basis with polarization and diffuse functions
            included.
        pcm_settings (Dict): If not None (default), then the polarizable continuum model
            (PCM) will be used to construct an implicit solvation environment around
            the molecule of interest. pcm_settings requires at least two keys,
            solvent and model. If solvent is "other", then four additional keys are
            required, dielectric (solvent dielectric constant), optical (square of the
            solvent index of refraction), density (solvent density in g/cm^3) and
            molar_mass (solvent molar mass in g/mol).
        max_scf_cycles (int): Maximum SCF cycles to be allowed before the calculation
            fails. The detault is 400, but for simple calculations, a much lower number
            (perhaps 100) is reasonable.
        geom_opt_max_cycles (int): Maximum number of cycles allowed before a
            geometry optimization fails. The default is 250, but for simple calculations,
            a much lower number (perhaps 100) is reasonable.
        overwrite_inputs_gen (Dict): Dictionary of Jaguar inputs that should overwrite
            defaults
    """

    def __init__(self,
                 molecule: Union[Molecule, MoleculeGraph],
                 name: str = "jaguar.in",
                 dft_rung: int = 4,
                 basis_set: str = "def2-svpd(-f)",
                 pcm_settings: Optional[Dict] = None,
                 max_scf_cycles: int = 400,
                 geom_opt_max_cycles: int = 250,
                 overwrite_inputs_gen: Optional[Dict] = None):

        if overwrite_inputs_gen is None:
            gen = dict()
        else:
            gen = overwrite_inputs_gen

        gen["igeopt"] = "1"

        if "maxitg" not in gen:
            gen["maxitg"] = geom_opt_max_cycles
        if "check_min" not in gen:
            gen["check_min"] = "1"  # Calculate Hessian to ensure convergence to a minimum
        # if "geoconv_mode" not in gen:
        #     gen["geoconv_mode"] = "standard"  # By default, do not allow "flexible" approach
        if "check_min_eigcut" not in gen:
            gen["check_min_eigcut"] = "-15.0"  # Ignore imaginary frequencies larger than -15 cm^-1

        super().__init__(molecule, name=name, dft_rung=dft_rung, basis_set=basis_set,
                         pcm_settings=pcm_settings, max_scf_cycles=max_scf_cycles,
                         overwrite_inputs_gen=gen)


class TSOptSet(JagSet):
    """
    JagSet object for use with transition-state optimization calculations

    Args;
         molecule (Union[Molecule, MoleculeGraph]): molecule object that will
            be the subject of this calculation. Will be converted to a Schrodinger
            Structure object.
            NOTE: If the quadratic synchronous transit (QST) method is being used,
            this molecule should be the TS guess. Reactant and product structures
            should be provided separately (see below).
        name (str): Name for this Jaguar input. Default is "jaguar.in"
        dft_rung (int): Which type of DFT functional should be used. Values
            between 1-4 are accepted, with default functionals being chosen:
            1: "hfs"
            2: "b97-d3"
            3: "m06-l"
            4: "wb97x-d"
        basis_set (str): Basis set to be used. Default is "def2-svpd(-f)", a relatively
            small split-valence basis with polarization and diffuse functions
            included.
        pcm_settings (Dict): If not None (default), then the polarizable continuum model
            (PCM) will be used to construct an implicit solvation environment around
            the molecule of interest. pcm_settings requires at least two keys,
            solvent and model. If solvent is "other", then four additional keys are
            required, dielectric (solvent dielectric constant), optical (square of the
            solvent index of refraction), density (solvent density in g/cm^3) and
            molar_mass (solvent molar mass in g/mol).
        max_scf_cycles (int): Maximum SCF cycles to be allowed before the calculation
            fails. The detault is 400, but for simple calculations, a much lower number
            (perhaps 100) is reasonable.
        geom_opt_max_cycles (int): Maximum number of cycles allowed before a
            geometry optimization fails. The default is 250, but for simple calculations,
            a much lower number (perhaps 100) is reasonable.
        use_qst (bool): If True (default False), perform a quadratic synchronous transit
            calculation. If this is True, then reactant_molecule and product_molecule
            must be set (see below).
        reactant_molecule (Optional[Union[Molecule, MoleculeGraph]]): Reactant molecule,
            used only for quadratic synchronous transit (QST) calculations. Default None.
        product_molecule (Optional[Union[Molecule, MoleculeGraph]]): Product molecule,
            used only for quadratic synchronous transit (QST) calculuations. Default None.
        use_analytic_hessian (bool): If True (default), then at the beginning of the
            calculation, an analytic Hessian will be calculated.
        overwrite_inputs_gen (Dict): Dictionary of Jaguar inputs that should overwrite
            defaults
    """

    def __init__(self,
                 molecule: Union[Molecule, MoleculeGraph],
                 name: str = "jaguar.in",
                 dft_rung: int = 4,
                 basis_set: str = "def2-svpd(-f)",
                 pcm_settings: Optional[Dict] = None,
                 max_scf_cycles: int = 400,
                 geom_opt_max_cycles: int = 250,
                 use_qst: bool = False,
                 reactant_molecule: Optional[Union[Molecule, MoleculeGraph]] = None,
                 product_molecule: Optional[Union[Molecule, MoleculeGraph]] = None,
                 use_analytic_hessian: bool = True,
                 overwrite_inputs_gen: Optional[Dict] = None):

        if overwrite_inputs_gen is None:
            gen = dict()
        else:
            gen = overwrite_inputs_gen

        gen["igeopt"] = "2"

        if "maxitg" not in gen:
            gen["maxitg"] = geom_opt_max_cycles
        if "iqst" not in gen and use_qst:
            if reactant_molecule is None or product_molecule is None:
                raise ValueError("For QST calculation, reactant_molecule and product_molecule must be set!")
            gen["iqst"] = "1"  # Use quadratic synchronous transit TS search method
        if "no_mul_imag_freq" not in gen:
            gen["no_mul_imag_freq"] = "1"  # Perturb TS structure to remove multiple imaginary frequencies
        if "inhess" not in gen and use_analytic_hessian:
            gen["inhess"] = "4"  # Calculate an analytic quantum mechanical Hessian initially

        super().__init__(molecule, name=name, dft_rung=dft_rung, basis_set=basis_set,
                         pcm_settings=pcm_settings, max_scf_cycles=max_scf_cycles,
                         overwrite_inputs_gen=gen)

        # TODO: test this
        if use_qst:
            if isinstance(reactant_molecule, Molecule):
                rct = molecule_to_schrodinger_struct(reactant_molecule)
            else:
                rct = mol_graph_to_schrodinger_struct(reactant_molecule)

            if isinstance(product_molecule, Molecule):
                pro = molecule_to_schrodinger_struct(product_molecule)
            else:
                pro = mol_graph_to_schrodinger_struct(product_molecule)

            self.jagin.setStructure(rct, zmat=1)
            self.jagin.setStructure(pro, zmat=2)


class FreqSet(JagSet):
    """
    JagSet object for use with vibrational frequency calculations

    Args;
         molecule (Union[Molecule, MoleculeGraph]): molecule object that will
            be the subject of this calculation. Will be converted to a Schrodinger
            Structure object.
        name (str): Name for this Jaguar input. Default is "jaguar.in"
        dft_rung (int): Which type of DFT functional should be used. Values
            between 1-4 are accepted, with default functionals being chosen:
            1: "hfs"
            2: "b97-d3"
            3: "m06-l"
            4: "wb97x-d"
        basis_set (str): Basis set to be used. Default is "def2-svpd(-f)", a relatively
            small split-valence basis with polarization and diffuse functions
            included.
        pcm_settings (Dict): If not None (default), then the polarizable continuum model
            (PCM) will be used to construct an implicit solvation environment around
            the molecule of interest. pcm_settings requires at least two keys,
            solvent and model. If solvent is "other", then four additional keys are
            required, dielectric (solvent dielectric constant), optical (square of the
            solvent index of refraction), density (solvent density in g/cm^3) and
            molar_mass (solvent molar mass in g/mol).
        max_scf_cycles (int): Maximum SCF cycles to be allowed before the calculation
            fails. The detault is 400, but for simple calculations, a much lower number
            (perhaps 100) is reasonable.
        overwrite_inputs_gen (Dict): Dictionary of Jaguar inputs that should overwrite
            defaults
    """

    def __init__(self,
                 molecule: Union[Molecule, MoleculeGraph],
                 name: str = "jaguar.in",
                 dft_rung: int = 4,
                 basis_set: str = "def2-svpd(-f)",
                 pcm_settings: Optional[Dict] = None,
                 max_scf_cycles: int = 400,
                 overwrite_inputs_gen: Optional[Dict] = None):

        if overwrite_inputs_gen is None:
            gen = dict()
        else:
            gen = overwrite_inputs_gen

        gen["ifreq"] = 1

        if "iraman" not in gen:
            gen["iraman"] = "1"  # Calculate Raman spectrum

        super().__init__(molecule, name=name, dft_rung=dft_rung, basis_set=basis_set,
                         pcm_settings=pcm_settings, max_scf_cycles=max_scf_cycles,
                         overwrite_inputs_gen=gen)


class ScanSet(JagSet):
    """
    JagSet object for use with potential energy surface scan calculations

    Args;
         molecule (Union[Molecule, MoleculeGraph]): molecule object that will
            be the subject of this calculation. Will be converted to a Schrodinger
            Structure object.
        scan_variables (list of dicts): Coordinates to be scanned, along with scan
            parameters. Ex:
            scan_variables=[{"coordinate_type": "distance", "atoms": [1, 2],
                             "initial": 1.0, "final": 2.0, "num_steps": 100, "step": 0.01}]
            Valid coordinate_type values are "distance", "angle", and "torsion".
            As many as 5 scan_variables can be given (at least 1 must be given for the
            scan to be valid), and atom numbers are 0-based
        name (str): Name for this Jaguar input. Default is "jaguar.in"
        dft_rung (int): Which type of DFT functional should be used. Values
            between 1-4 are accepted, with default functionals being chosen:
            1: "hfs"
            2: "b97-d3"
            3: "m06-l"
            4: "wb97x-d"
        basis_set (str): Basis set to be used. Default is "def2-svpd(-f)", a relatively
            small split-valence basis with polarization and diffuse functions
            included.
        pcm_settings (Dict): If not None (default), then the polarizable continuum model
            (PCM) will be used to construct an implicit solvation environment around
            the molecule of interest. pcm_settings requires at least two keys,
            solvent and model. If solvent is "other", then four additional keys are
            required, dielectric (solvent dielectric constant), optical (square of the
            solvent index of refraction), density (solvent density in g/cm^3) and
            molar_mass (solvent molar mass in g/mol).
        max_scf_cycles (int): Maximum SCF cycles to be allowed before the calculation
            fails. The detault is 400, but for simple calculations, a much lower number
            (perhaps 100) is reasonable.
        geom_opt_max_cycles (int): Maximum number of cycles allowed before a
            geometry optimization fails. The default is 250, but for simple calculations,
            a much lower number (perhaps 100) is reasonable.
        overwrite_inputs_gen (Dict): Dictionary of Jaguar inputs that should overwrite
            defaults
    """

    def __init__(self,
                 molecule: Union[Molecule, MoleculeGraph],
                 scan_variables: List[Dict],
                 name: str = "jaguar.in",
                 dft_rung: int = 4,
                 basis_set: str = "def2-svpd(-f)",
                 pcm_settings: Optional[Dict] = None,
                 max_scf_cycles: int = 400,
                 geom_opt_max_cycles: int = 250,
                 overwrite_inputs_gen: Optional[Dict] = None):

        if len(scan_variables) == 0 or len(scan_variables) > 5:
            raise ValueError("Between 1 and 5 scan_variables must be provided!")

        if overwrite_inputs_gen is None:
            gen = dict()
        else:
            gen = overwrite_inputs_gen

        if "maxitg" not in gen:
            gen["maxitg"] = geom_opt_max_cycles

        super().__init__(molecule, name=name, dft_rung=dft_rung,
                         basis_set=basis_set, pcm_settings=pcm_settings,
                         max_scf_cycles=max_scf_cycles, overwrite_inputs_gen=gen)

        for var in scan_variables:
            if not all([x in var for x in ["coordinate_type", "atoms", "initial", "final", "num_steps", "step"]]):
                raise ValueError("Invalid scan_variable format! Must have keys 'coordinate_type',"
                                 " 'atoms', 'initial', 'final', 'num_steps', and 'step'")

            if var["coordinate_type"].lower() not in ["distance", "angle", "torsion"]:
                raise ValueError("Invalid coordinate type for scan_variables! Must be one "
                                 "of 'distance', 'angle', or 'torsion'!")

            if var["coordinate_type"].lower() == "distance":
                ct = mm.MMJAG_COORD_DISTANCE
            elif var["coordinate_type"].lower() == "angle":
                ct = mm.MMJAG_COORD_ANGLE
            else:
                ct = mm.MMJAG_COORD_TORSION

            self.jagin.setScanCoordinate(ct,
                                         [a + 1 for a in var["atoms"]],
                                         var["initial"],
                                         var["final"],
                                         var["num_steps"],
                                         var["step"])


class IRCSet(JagSet):
    """
    JagSet object for use with intrinsic reaction coordinate (IRC) calculations.

    Args;
         molecule (Union[Molecule, MoleculeGraph]): molecule object that will
            be the subject of this calculation. Will be converted to a Schrodinger
            Structure object.
        name (str): Name for this Jaguar input. Default is "jaguar.in"
        dft_rung (int): Which type of DFT functional should be used. Values
            between 1-4 are accepted, with default functionals being chosen:
            1: "hfs"
            2: "b97-d3"
            3: "m06-l"
            4: "wb97x-d"
        basis_set (str): Basis set to be used. Default is "def2-svpd(-f)", a relatively
            small split-valence basis with polarization and diffuse functions
            included.
        pcm_settings (Dict): If not None (default), then the polarizable continuum model
            (PCM) will be used to construct an implicit solvation environment around
            the molecule of interest. pcm_settings requires at least two keys,
            solvent and model. If solvent is "other", then four additional keys are
            required, dielectric (solvent dielectric constant), optical (square of the
            solvent index of refraction), density (solvent density in g/cm^3) and
            molar_mass (solvent molar mass in g/mol).
        max_scf_cycles (int): Maximum SCF cycles to be allowed before the calculation
            fails. The detault is 400, but for simple calculations, a much lower number
            (perhaps 100) is reasonable.
        endpoints_only (bool): If True (default False), then perform a "three-point"
            IRC, where the goal is to identify the endpoints, not to characterize
            the entire reaction pathway.
        overwrite_inputs_gen (Dict): Dictionary of Jaguar inputs that should overwrite
            defaults
    """

    def __init__(self,
                 molecule: Union[Molecule, MoleculeGraph],
                 name: str = "jaguar.in",
                 dft_rung: int = 4,
                 basis_set: str = "def2-svpd(-f)",
                 pcm_settings: Optional[Dict] = None,
                 max_scf_cycles: int = 400,
                 endpoints_only: bool = False,
                 overwrite_inputs_gen: Optional[Dict] = None):

        if overwrite_inputs_gen is None:
            gen = dict()
        else:
            gen = overwrite_inputs_gen

        gen["irc"] = "1"
        gen["irc_grad_check"] = "0"

        if "inhess" not in gen:
            gen["inhess"] = 4

        if endpoints_only:
            if "three_pt_irc" not in gen:
                gen["three_pt_irc"] = "1"  # Only calculate endpoints, and not entire reaction pathway
        else:
            if "ircmax" not in gen:
                gen["ircmax"] = "100"  # Allow for up to 10 points along IRC on each side of TS
            if "ircmxcyc" not in gen:
                gen["ircmxcyc"] = "300"  # Maximum number of steps for optimization of each point along IRC

        super().__init__(molecule, name=name, dft_rung=dft_rung, basis_set=basis_set, pcm_settings=pcm_settings,
                         max_scf_cycles=max_scf_cycles, overwrite_inputs_gen=gen)