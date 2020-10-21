# coding: utf-8

import os
import unittest
import copy

import numpy as np

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN, CovalentBondNN
from pymatgen.analysis.fragmenter import metal_edge_extender

from mpcat.utils.generate import mol_to_mol_graph


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files")


class GenerateTest(unittest.TestCase):

    def test_mol_to_mol_graph(self):
        mol = Molecule.from_file(os.path.join(test_dir, "molecules", "li2co3_1.xyz"))
        mg = MoleculeGraph.with_local_env_strategy(mol, OpenBabelNN())
        mg = metal_edge_extender(mg)

        self.assertEqual(mg, mol_to_mol_graph(mol))


if __name__ == "__main__":
    unittest.main()
