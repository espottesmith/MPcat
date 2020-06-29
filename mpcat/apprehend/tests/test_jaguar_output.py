import os
import unittest

from monty.serialization import loadfn, dumpfn

from pymatgen import Molecule
from mpcat.apprehend.jaguar_output import JagOutput
from mpcat.adapt.schrodinger_adapter import (file_to_schrodinger_structure,
                                             molecule_to_maestro_file,
                                             maestro_file_to_molecule,
                                             schrodinger_struct_to_molecule,
                                             molecule_to_schrodinger_struct)

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files")


def prepare_files():
    pass


class TestJagOutput(unittest.TestCase):

    def test_success(self):
        pass

    def test_failure(self):
        pass
