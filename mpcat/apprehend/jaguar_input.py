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


class JagSet(JagInput):
    pass


class OptSet(JagSet):
    pass


class TSSet(JagSet):
    pass


class FreqSet(JagSet):
    pass


class ScanSet(JagSet):
    pass


class IRCSet(JagSet)
    pass