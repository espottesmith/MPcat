# coding: utf-8

import unittest
import copy
from pathlib import Path

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN, CovalentBondNN
from pymatgen.analysis.fragmenter import metal_edge_extender

from mpcat.utils.reaction import (get_reaction_graphs,
                                  union_molgraph,
                                  get_atom_mappings,
                                  get_reaction_template)


test_dir = Path(__file__).resolve().parent.parent.parent / "test_files"


class ReactionTest(unittest.TestCase):

    def test_get_reaction_graphs(self):

        # Test that isomorphic graphs return themselves
        li2co3 = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "molecules" / "li2co3_1.xyz").as_posix()),
            OpenBabelNN())
        li2co3 = metal_edge_extender(li2co3)
        li2co3_copy = copy.deepcopy(li2co3)

        self.assertEqual(get_reaction_graphs(li2co3, li2co3_copy), [li2co3])

        # Test for case where graphs cannot be compared (different species, etc.)
        ethane = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "molecules" / "ethane.mol").as_posix()),
            OpenBabelNN())

        self.assertEqual(get_reaction_graphs(li2co3, ethane), list())

        # Test for case of a single bond-breaking
        break_rct = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "reactions" / "isomorphism" / "break" / "rct.xyz").as_posix()),
            OpenBabelNN())
        break_rct = metal_edge_extender(break_rct)
        break_pro = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "reactions" / "isomorphism" / "break" / "pro.xyz").as_posix()),
            OpenBabelNN())
        break_pro = metal_edge_extender(break_pro)

        break_rxn = get_reaction_graphs(break_rct, break_pro)
        self.assertEqual(len(break_rxn), 1)
        self.assertTrue(break_rxn[0].isomorphic_to(break_rct))
        # for edge in break_rxn[0].graph.edges(data=True):
        #     print(edge)
        self.assertTrue(break_rxn[0].graph[4][5][0].get("broken", False))
        self.assertFalse(break_rxn[0].graph[4][5][0].get("formed", True))

        # Test for case of a single bond-forming
        form_rct = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "reactions" / "isomorphism" / "form" / "rct.xyz").as_posix()),
            OpenBabelNN())
        form_rct = metal_edge_extender(form_rct)
        form_pro = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "reactions" / "isomorphism" / "form" / "pro.xyz").as_posix()),
            OpenBabelNN())
        form_pro = metal_edge_extender(form_pro)

        form_rxn = get_reaction_graphs(form_rct, form_pro)
        self.assertEqual(len(form_rxn), 1)
        self.assertTrue(form_rxn[0].isomorphic_to(form_pro))
        self.assertTrue(form_rxn[0].graph[0][3][0].get("formed", False))
        self.assertFalse(form_rxn[0].graph[0][3][0].get("broken", True))

        # Test for case of complex bond breaking and formation
        bf_rct = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "reactions" / "isomorphism" / "break_form" / "rct.xyz").as_posix()),
            OpenBabelNN())
        bf_rct = metal_edge_extender(bf_rct)
        bf_pro = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "reactions" / "isomorphism" / "break_form" / "pro.xyz").as_posix()),
            OpenBabelNN())
        bf_pro = metal_edge_extender(bf_pro)

        bf_rxn = get_reaction_graphs(bf_rct, bf_pro)
        self.assertEqual(len(bf_rxn), 0)

        bf_rxn = get_reaction_graphs(bf_rct, bf_pro, allowed_break=2, allowed_form=1,
                                     stop_at_one=True)
        self.assertEqual(len(bf_rxn), 1)

        bf_rxn = get_reaction_graphs(bf_rct, bf_pro, allowed_break=2, allowed_form=1,
                                     stop_at_one=False)
        self.assertEqual(len(bf_rxn), 2)

        self.assertTrue(bf_rxn[0].graph[0][8][0].get("broken", False))
        self.assertTrue(bf_rxn[0].graph[6][12][0].get("broken", False))
        self.assertTrue(bf_rxn[0].graph[8][12][0].get("formed", False))

        self.assertTrue(bf_rxn[1].graph[0][9][0].get("broken", False))
        self.assertTrue(bf_rxn[1].graph[6][12][0].get("broken", False))
        self.assertTrue(bf_rxn[1].graph[9][12][0].get("formed", False))

    def test_union_molgraph(self):

        with self.assertRaises(ValueError):
            _ = union_molgraph([])

        # Test "good", well-behaved case
        good_one = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "molecules" / "union" / "good" / "1.xyz").as_posix()),
            OpenBabelNN())
        good_one = metal_edge_extender(good_one)
        good_two = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "molecules" / "union" / "good" / "2.xyz").as_posix()),
            OpenBabelNN())
        good_two = metal_edge_extender(good_two)
        good_union = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "molecules" / "union" / "good" / "union.xyz").as_posix()),
            OpenBabelNN())
        good_union = metal_edge_extender(good_union)

        good = union_molgraph([good_one, good_two])
        self.assertTrue(good_union.isomorphic_to(good))

        # Test "bad" case where proximity might be an issue
        bad_one = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "molecules" / "union" / "bad" / "1.xyz").as_posix()),
            OpenBabelNN())
        bad_one = metal_edge_extender(bad_one)
        bad_two = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "molecules" / "union" / "bad" / "2.xyz").as_posix()),
            OpenBabelNN())
        bad_two = metal_edge_extender(bad_two)
        bad_union = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "molecules" / "union" / "bad" / "union.xyz")),
            OpenBabelNN())
        bad_union = metal_edge_extender(bad_union)

        bad = union_molgraph([bad_one, bad_two])
        self.assertTrue(bad_union.isomorphic_to(bad))

        with self.assertRaises(ValueError):
            _ = union_molgraph([bad_one, bad_two], validate_proximity=True)

    def test_get_atom_mappings(self):

        rct = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "reactions" / "da" / "rct.xyz").as_posix()),
            CovalentBondNN())
        pro = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "reactions" / "da" / "pro.xyz").as_posix()),
            CovalentBondNN())

        # Test with no valid mapping
        mapping = get_atom_mappings(rct, pro, allowed_form=1, allowed_break=0)
        self.assertEqual(len(mapping), 0)

        # Test with stop_at_one
        mapping = get_atom_mappings(rct, pro, allowed_form=2, allowed_break=0)
        self.assertEqual(len(mapping), 1)

        # Test without stop_at_one
        mapping = get_atom_mappings(rct, pro, allowed_form=2, allowed_break=0,
                                    stop_at_one=False)
        self.assertEqual(len(mapping), 4)

        # Test with give_best
        mapping = get_atom_mappings(rct, pro, allowed_form=2, allowed_break=0,
                                    stop_at_one=False, give_best=True)
        self.assertEqual(len(mapping), 1)
        self.assertDictEqual(mapping[0], {2: 0, 5: 1, 11: 2, 0: 3, 12: 4, 4: 5,
                                          10: 6, 13: 7, 1: 8, 3: 9, 14: 10,
                                          6: 11, 8: 12, 9: 13, 15: 14, 16: 15,
                                          17: 16, 7: 17, 18: 18})

    def test_get_reaction_template(self):
        rct = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "reactions" / "da" / "rct.xyz").as_posix()),
            CovalentBondNN())
        pro = MoleculeGraph.with_local_env_strategy(Molecule.from_file(
            (test_dir / "reactions" / "da" / "pro.xyz").as_posix()),
            CovalentBondNN())

        rg = get_reaction_graphs(rct, pro, allowed_form=2,
                                 allowed_break=0, stop_at_one=True)[0]

        template_0 = get_reaction_template(rg)
        self.assertSetEqual(set(template_0.nodes), {2, 5, 11, 12})
        self.assertEqual(len(template_0.edges), 3)
        template_1 = get_reaction_template(rg, order=1)
        self.assertSetEqual(set(template_1.nodes), {0, 2, 4, 5, 8, 9, 10,
                                                    11, 12, 13, 14, 15, 16})
        self.assertEqual(len(template_1.edges), 14)
        template_2 = get_reaction_template(rg, order=2)
        self.assertSetEqual(set(template_2.nodes), {0, 1, 2, 3, 4, 5, 6, 8, 9,
                                                    10, 11, 12, 13, 14, 15, 16,
                                                    17, 18})
        self.assertEqual(len(template_2.edges), 20)
        template_3 = get_reaction_template(rg, order=3)
        self.assertSetEqual(set(template_3.nodes), {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                                    10, 11, 12, 13, 14, 15, 16,
                                                    17, 18})
        self.assertEqual(len(template_3.edges), 21)


if __name__ == "__main__":
    unittest.main()
