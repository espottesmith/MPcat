# coding: utf-8

import unittest
from pathlib import Path

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN, CovalentBondNN
from pymatgen.analysis.fragmenter import metal_edge_extender

from mpcat.utils.generate import mol_to_mol_graph


test_dir = Path(__file__).resolve().parent.parent.parent / "test_files"


class GenerateTest(unittest.TestCase):

    def test_mol_to_mol_graph(self):
        mol = Molecule.from_file((test_dir / "molecules" / "li2co3_1.xyz").as_posix())
        mg = MoleculeGraph.with_local_env_strategy(mol, OpenBabelNN())
        mg = metal_edge_extender(mg)

        self.assertEqual(mg, mol_to_mol_graph(mol))


if __name__ == "__main__":
    unittest.main()
