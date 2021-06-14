# coding: utf-8
# Distributed under the terms of the MIT License.


import logging
import os
import unittest

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.util.testing import PymatgenTest

from mpcat.automate.atomate.inputs import QCTemplate, GSMIsomerInput

__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2021"
__version__ = "0.1"
__email__ = "espottesmith@gmail.com"

logger = logging.getLogger(__name__)

test_dir = os.path.join(os.path.dirname(__file__), "test_files")


class TestQCTemplate(PymatgenTest):

    def test_create(self):
        # Test without basis
        rem = {"xc_grid": 3,
               "max_scf_cycles": 200,
               "scf_algorithm": "diis",
               "thresh": 14}

        with self.assertRaises(ValueError):
            test = QCTemplate(rem)

        # Test without method
        rem["basis"] = "6-311++g(d,p)"

        with self.assertRaises(ValueError):
            test = QCTemplate(rem)

        # Test with job type other than force
        rem["method"] = "wb97xd"
        rem["job_type"] = "fsm"

        test = QCTemplate(rem)

        self.assertEqual(test.rem["job_type"], "force")

        # Test correct input
        rem["job_type"] = "force"
        smx = {"solvent": "thf"}

        test = QCTemplate(rem, smx=smx)
        self.assertDictEqual(test.rem, rem)
        self.assertDictEqual(test.smx, smx)

    def test_str(self):
        rem = {"basis": "6-31*",
               "method": "b3lyp",
               "job_type": "force"}
        pcm = {"theory": "cpcm"}
        solvent = {"dielectric": 80.4}

        test = QCTemplate(rem, pcm=pcm, solvent=solvent)

        self.assertEqual(str(test),
                         '$rem\n   basis = 6-31*\n   method = b3lyp\n   job_type = force\n$end\n\n$pcm\n   theory cpcm\n$end\n\n$solvent\n   dielectric 80.4\n$end\n\n$molecule\n')

    def test_from_file(self):
        from_file = QCTemplate.from_file(os.path.join(test_dir, "qin_good"))

        self.assertDictEqual(from_file.rem, {"job_type": "force",
                                             "basis": "def2-tzvppd",
                                             "max_scf_cycles": "200",
                                             "gen_scfman": "true",
                                             "xc_grid": "3",
                                             "scf_algorithm": "diis",
                                             "resp_charges": "true",
                                             "symmetry": "false",
                                             "sym_ignore": "true",
                                             "method": "wb97x-v",
                                             "solvent_method": "smd",
                                             "ideriv": "1",
                                             "thresh": "14"})
        self.assertDictEqual(from_file.smx, {"solvent": "other"})

    def test_from_qcinput(self):
        from_qcinp = QCTemplate.from_file(os.path.join(test_dir, "qchem_input.qin"))
        self.assertDictEqual(from_qcinp.rem, {'job_type': 'force',
                                              'basis': 'def2-tzvppd',
                                              'max_scf_cycles': '200',
                                              'gen_scfman': 'true',
                                              'scf_algorithm': 'diis',
                                              'method': 'wb97xd',
                                              'geom_opt_max_cycles': '200',
                                              'xc_grid': '3',
                                              'symmetry': 'false',
                                              'sym_ignore': 'true',
                                              'resp_charges': 'true'})


class TestGSMIsomerInput(PymatgenTest):

    def test_create(self):
        # Test with too many coordinates
        with self.assertRaises(ValueError):
            too_many_coords = GSMIsomerInput(bonds_formed=[(0, 1), (1, 2)],
                                             angles=[(4, 5, 6), (9, 10, 11)],
                                             torsions=[(3, 4, 7, 8)])

        # Test with non-integer indices
        with self.assertRaises(ValueError):
            non_integer = GSMIsomerInput(bonds_broken=[("pi", "three")])

        # Test with wrong number of indices
        with self.assertRaises(ValueError):
            wrong_indices = GSMIsomerInput(bonds_formed=[(1, 2, 3)])
        with self.assertRaises(ValueError):
            wrong_indices = GSMIsomerInput(bonds_broken=[(1, 2, 3)])
        with self.assertRaises(ValueError):
            wrong_indices = GSMIsomerInput(angles=[(1, 3)])
        with self.assertRaises(ValueError):
            wrong_indices = GSMIsomerInput(torsions=[(1, 2, 3)])
        with self.assertRaises(ValueError):
            wrong_indices = GSMIsomerInput(out_of_planes=[(1, 2, 3, 4, 5)])

        # Test with molecule - indices too high
        mol = Molecule.from_file(os.path.join(test_dir, "ethane.mol"))
        with self.assertRaises(ValueError):
            too_high = GSMIsomerInput(molecule=mol, bonds_broken=[(7, 9)])

        # Test good
        good = GSMIsomerInput(molecule=mol, bonds_broken=[(0, 1)])
        self.assertEqual(good.molecule, mol)
        self.assertEqual(good.bonds_broken, [(0, 1)])

    def test_verify_with_graphs(self):
        ethane = Molecule.from_file(os.path.join(test_dir, "ethane.mol"))

        # Test bad bond formed
        with self.assertRaises(ValueError):
            bad_bond = GSMIsomerInput(molecule=ethane, bonds_formed=[(0, 1)],
                                      use_graph=True)

        # Test bad bond broken
        with self.assertRaises(ValueError):
            bad_bond = GSMIsomerInput(molecule=ethane, bonds_broken=[(1, 2)],
                                      use_graph=True)

        # Test bad angle
        with self.assertRaises(ValueError):
            bad_angle = GSMIsomerInput(molecule=ethane,
                                       angles=[(0, 1, 2)],
                                       use_graph=True)

        # Test bad torsion
        with self.assertRaises(ValueError):
            bad_torsion = GSMIsomerInput(molecule=ethane,
                                         torsions=[(0, 1, 2, 3)],
                                         use_graph=True)

        # Test bad out of plane bend
        with self.assertRaises(ValueError):
            bad_out_of_plane = GSMIsomerInput(molecule=ethane,
                                              out_of_planes=[(0, 1, 2, 3)],
                                              use_graph=True)

        # Test good
        good = GSMIsomerInput(molecule=ethane,
                              bonds_formed=[(0, 7)],
                              angles=[(1, 0, 4)],
                              torsions=[(2, 0, 1, 6)],
                              use_graph=True)
        mg = MoleculeGraph.with_local_env_strategy(ethane, OpenBabelNN())
        self.assertEqual(mg, good.molecule_graph)

    def test_str(self):
        len_two = GSMIsomerInput(bonds_formed=[(0, 1)],
                                 bonds_broken=[(2, 3), (5, 8)])

        self.assertEqual(str(len_two), "ADD 1 2\nBREAK 3 4\nBREAK 6 9")

        len_three_four = GSMIsomerInput(angles=[(1, 2, 3)],
                                        torsions=[(4, 5, 8, 9)],
                                        out_of_planes=[(1, 11, 68, 3)])

        self.assertEqual(str(len_three_four),
                         "ANGLE 2 3 4\nTORSION 5 6 9 10\nOOP 2 12 69 4")

    def test_to_from_file(self):
        len_two = GSMIsomerInput(bonds_formed=[(0, 1)],
                                 bonds_broken=[(2, 3), (5, 8)])
        len_two.write_file("test")

        check = GSMIsomerInput.from_file("test")
        self.assertEqual(len_two.bonds_broken, check.bonds_broken)
        self.assertEqual(len_two.bonds_formed, check.bonds_formed)
        os.remove("test")


if __name__ == "__main__":
    unittest.main()
