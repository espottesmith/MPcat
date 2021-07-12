# coding: utf-8

from typing import List, Dict, Optional, Union
import copy
from pathlib import Path
import os

import numpy as np

from monty.json import MSONable
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph

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
    gen["matit"] = max_scf_cycles

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
            self.jagin = JaguarInput(name=name, stucture=self.struct)

        self.modify_gen(gen_variables)

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

        gen["igeopt"] = 1

        if "maxitg" not in gen:
            gen["maxitg"] = geom_opt_max_cycles
        if "check_min" not in gen:
            gen["check_min"] = 1  # Calculate Hessian to ensure convergence to a minimum
        # if "geoconv_mode" not in gen:
        #     gen["geoconv_mode"] = "standard"  # By default, do not allow "flexible" approach
        if "check_min_eigcut" not in gen:
            gen["check_min_eigcut"] = -15.0  # Ignore imaginary frequencies larger than -15 cm^-1

        super().__init__(molecule, name=name, dft_rung=dft_rung, basis_set=basis_set,
                         pcm_settings=pcm_settings, max_scf_cycles=max_scf_cycles,
                         overwrite_inputs_gen=gen)


class TSSet(JagSet):
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
                 use_qst: bool = False,
                 use_analytic_hessian: bool = True,
                 overwrite_inputs_gen: Optional[Dict] = None):

        if overwrite_inputs_gen is None:
            gen = dict()
        else:
            gen = overwrite_inputs_gen

        gen["igeopt"] = 2

        if "maxitg" not in gen:
            gen["maxitg"] = geom_opt_max_cycles
        if "iqst" not in gen and use_qst:
            gen["iqst"] = 1  # Use quadratic synchronous transit TS search method
        if "no_mul_imag_freq" not in gen:
            gen["no_mul_imag_freq"] = 1  # Perturb TS structure to remove multiple imaginary frequencies
        if "inhess" not in gen and use_analytic_hessian:
            gen["inhess"] = 4  # Calculate an analytic quantum mechanical Hessian initially

        super().__init__(molecule, name=name, dft_rung=dft_rung, basis_set=basis_set,
                         pcm_settings=pcm_settings, max_scf_cycles=max_scf_cycles,
                         overwrite_inputs_gen=gen)


class FreqSet(JagSet):
    pass


class ScanSet(JagSet):
    pass


class IRCSet(JagSet):
    pass