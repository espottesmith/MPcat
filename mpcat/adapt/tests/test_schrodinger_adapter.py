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

from mpcat.adapt.schrodinger_adapter import (file_to_schrodinger_structure,
                                             maestro_file_to_molecule,
                                             schrodinger_struct_to_molecule,
                                             schrodinger_struct_to_mol_graph,
                                             molecule_to_schrodinger_struct,
                                             mol_graph_to_schrodinger_struct)

module_dir = Path(__file__).resolve().parent
test_dir = Path(__file__).resolve().parent.parent.parent.parent / "test_files"
molecule_dir = test_dir / "molecules"


class TestSchrodingerAdapter(unittest.TestCase):

    def setUp(self) -> None:
        if ob:
            self.molecule = Molecule.from_file((molecule_dir / "ethane.mol").as_posix())
            self.molecule.to("sdf", (molecule_dir / "ethane.sdf").as_posix())

    @unittest.skipIf(not ob, "Openbabel not present. Skipping...")
    def test_file_to_schrodinger_structure(self):
        structures = file_to_schrodinger_structure(molecule_dir / "ethane.sdf")
        self.assertEqual(len(structures), 1)

        elements = [a.element for a in structures[0].molecule[1].atom]
        self.assertSequenceEqual(elements, ["C", "C", "H", "H", "H", "H", "H", "H"])

        positions = [[-0.7520, 0.0010, -0.1410],
                     [0.7520, -0.0010, 0.1410],
                     [-1.1580, 0.9910, 0.0700],
                     [-1.2400, -0.7370, 0.4960],
                     [-0.9240, -0.2490, -1.1880],
                     [1.1580, -0.9910, -0.0700],
                     [0.9240, 0.2490, 1.1880],
                     [1.2400, 0.7370, -0.4960]]

        for ii, atom in enumerate(structures[0].molecule[1].atom):
            self.assertListEqual(positions[ii],
                                 list(atom.xyz))

    @unittest.skipIf(not ob, "Openbabel not present. Skipping...")
    def test_maestro_file_to_molecule(self):
        molecules = maestro_file_to_molecule(molecule_dir / "ec.01.mae")
        self.assertEqual(len(molecules), 1)

        elements = [str(s) for s in molecules[0].species]
        self.assertSequenceEqual(elements, ["O", "C", "C", "O", "O", "C",
                                            "H", "H", "H", "H"])

        positions = [[0.292719, -1.186596, -0.240504],
                     [-0.708427, -0.361958, 0.069717],
                     [1.542315, -0.497739, -0.080547],
                     [-1.855122, -0.696497, 0.151131],
                     [-0.274057, 0.881704, 0.277870],
                     [1.116388, 0.963738, -0.070584],
                     [2.189227, -0.755905, -0.912950],
                     [1.987489, -0.814710, 0.861406],
                     [1.194818, 1.433330, -1.050016],
                     [1.626845, 1.559988, 0.679031]]

        for ii, line in enumerate(positions):
            self.assertListEqual(line, list(molecules[0].cart_coords[ii]))

    @unittest.skipIf(not ob, "Openbabel not present. Skipping...")
    def test_schrodinger_struct_to_molecule(self):
        struct = file_to_schrodinger_structure(molecule_dir / "ethane.sdf")[0]
        for atom in struct.atom:
            atom.formal_charge = 0
        struct.atom[1].formal_charge = -1
        mol = schrodinger_struct_to_molecule(struct)

        self.assertListEqual([a.element for a in struct.molecule[1].atom],
                             [str(s) for s in mol.species])

        self.assertEqual(mol.charge, struct.formal_charge)
        self.assertEqual(mol.charge, -1)

        for ii in range(len(mol)):
            self.assertListEqual(list(mol.cart_coords[ii]),
                                 list(struct.molecule[1].atom[ii + 1].xyz))

    @unittest.skipIf(not ob, "Openbabel not present. Skipping...")
    def test_schrodinger_struct_to_mol_graph(self):
        struct = file_to_schrodinger_structure(molecule_dir / "ethane.sdf")[0]
        for atom in struct.atom:
            atom.formal_charge = 0
        struct.atom[1].formal_charge = -1
        mg = schrodinger_struct_to_mol_graph(struct)

        self.assertListEqual([a.element for a in struct.molecule[1].atom],
                             [str(s) for s in mg.molecule.species])

        self.assertEqual(mg.molecule.charge, struct.formal_charge)
        self.assertEqual(mg.molecule.charge, -1)

        for ii in range(len(mg.molecule)):
            self.assertListEqual(list(mg.cart_coords[ii]),
                                 list(struct.molecule[1].atom[ii + 1].xyz))

        struct_bonds = set()
        for bond in struct.bond:
            struct_bonds.add(tuple(sorted([bond.atom1.index - 1, bond.atom2.index])))

        mg_bonds = set()
        for bond in mg.graph.edges():
            mg_bonds.add(tuple(sorted([bond[0], bond[1]])))

        self.assertSetEqual(struct_bonds, mg_bonds)

    @unittest.skipIf(not ob, "Openbabel not present. Skipping...")
    def test_molecule_to_schrodinger_struct(self):
        self.molecule.set_charge_and_spin(charge=-1)
        struct = molecule_to_schrodinger_struct(self.molecule)

        self.assertListEqual([a.element for a in struct.molecule[1].atom],
                             [str(s) for s in self.molecule.species])

        self.assertEqual(struct.formal_charge, -1)

        for ii in range(len(self.molecule)):
            self.assertListEqual(list(self.molecule.cart_coords[ii]),
                                 list(struct.molecule[1].atom[ii + 1].xyz))

    @unittest.skipIf(not ob, "Openbabel not present. Skipping...")
    def test_mol_graph_to_schrodinger_struct(self):
        self.molecule.set_charge_and_spin(charge=-1)
        mg = MoleculeGraph.with_local_env_strategy(self.molecule, OpenBabelNN())

        struct = mol_graph_to_schrodinger_struct(mg)

        self.assertListEqual([a.element for a in struct.molecule[1].atom],
                             [str(s) for s in mg.molecule.species])

        self.assertEqual(struct.formal_charge, -1)

        for ii in range(len(mg.molecule)):
            self.assertListEqual(list(mg.molecule.cart_coords[ii]),
                                 list(struct.molecule[1].atom[ii + 1].xyz))

        struct_bonds = set()
        for bond in struct.bond:
            struct_bonds.add(tuple(sorted([bond.atom1.index - 1, bond.atom2.index])))

        mg_bonds = set()
        for bond in mg.graph.edges():
            mg_bonds.add(tuple(sorted([bond[0], bond[1]])))

        self.assertSetEqual(struct_bonds, mg_bonds)


if __name__ == "__main__":
    unittest.main()
