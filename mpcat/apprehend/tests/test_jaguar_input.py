# coding: utf-8

from pathlib import Path

import unittest

try:
    from openbabel import openbabel as ob
except ImportError:
    ob = None

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from mpcat.apprehend.jaguar_input import (JagInput,
                                          JagSet,
                                          OptSet,
                                          TSOptSet,
                                          FreqSet,
                                          ScanSet,
                                          IRCSet)

module_dir = Path(__file__).resolve().parent
test_dir = Path(__file__).resolve().parent.parent.parent.parent / "test_files"
input_dir = test_dir / "jaguar_inputs"


class TestJagInput(unittest.TestCase):
    pass


class TestJagSet(unittest.TestCase):
    pass


class TestOptSet(unittest.TestCase):
    pass


class TestTSOptSet(unittest.TestCase):
    pass


class TestFreqSet(unittest.TestCase):
    pass


class TestScanSet(unittest.TestCase):
    pass


class TestIRCSet(unittest.TestCase):
    pass