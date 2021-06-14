# coding: utf-8
# Distributed under the terms of the MIT License.

import os
import unittest

from pymatgen.util.testing import PymatgenTest
from mpcat.automate.atomate.sets import GSMDictSet, QCTemplate

__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2021"
__version__ = "0.1"


class GSMDictSetTest(PymatgenTest):
    def test_init(self):
        test_DictSet = GSMDictSet(
            basis_set='6-31G*',
            scf_algorithm='diis')
        self.assertEqual(
            test_DictSet.rem, {
                'job_type': 'force',
                'gen_scfman': 'true',
                'basis': '6-31g*',
                'max_scf_cycles': 200,
                'method': 'wb97xv',
                'scf_algorithm': 'diis',
                'xc_grid': '3',
                'symmetry': 'false',
                'sym_ignore': 'true',
                'resp_charges': 'true'
            })
        self.assertEqual(test_DictSet.pcm, {})
        self.assertEqual(test_DictSet.solvent, {})
        self.assertEqual(test_DictSet.smx, {})

    def test_full_init(self):

        test_DictSet = GSMDictSet(
            basis_set='6-31g*',
            scf_algorithm='diis',
            dft_rung=1,
            pcm_dielectric=10.0,
            max_scf_cycles=35)
        self.assertEqual(
            test_DictSet.rem, {
                'job_type': 'force',
                'gen_scfman': 'true',
                'basis': '6-31g*',
                'max_scf_cycles': 35,
                'method': 'b3lyp',
                'scf_algorithm': 'diis',
                'xc_grid': '3',
                'solvent_method': 'pcm',
                'symmetry': 'false',
                'sym_ignore': 'true',
                'resp_charges': 'true'
            })
        self.assertEqual(
            test_DictSet.pcm, {
                'heavypoints': '194',
                'hpoints': '194',
                'radii': 'uff',
                'theory': 'cpcm',
                'vdwscale': '1.1'
            })
        self.assertEqual(test_DictSet.solvent, {'dielectric': 10.0})

        test_DictSet = GSMDictSet(
            basis_set='6-31g*',
            scf_algorithm='diis',
            dft_rung=1,
            smd_solvent='water',
            max_scf_cycles=35)
        self.assertEqual(
            test_DictSet.rem, {
                'job_type': 'force',
                'gen_scfman': 'true',
                'basis': '6-31g*',
                'max_scf_cycles': 35,
                'method': 'b3lyp',
                'scf_algorithm': 'diis',
                'xc_grid': '3',
                'solvent_method': 'smd',
                'ideriv': '1',
                'symmetry': 'false',
                'sym_ignore': 'true',
                'resp_charges': 'true'
            })
        self.assertEqual(test_DictSet.smx, {'solvent': 'water'})

    def test_overwrite_input(self):
        overwrite_inputs = {
            "rem": {
                'method': 'b3lyp',
                'basis': '6-31g*',
                'thresh': 10,
                "xc_grid": "000150000302"
            }
        }
        test_set = GSMDictSet(overwrite_inputs=overwrite_inputs)
        act_rem = {
            'job_type': 'force',
            'gen_scfman': 'true',
            'basis': '6-31g*',
            'max_scf_cycles': 200,
            'method': 'b3lyp',
            'scf_algorithm': 'diis',
            'xc_grid': '000150000302',
            'thresh': 10,
            'symmetry': 'false',
            'sym_ignore': 'true',
            'resp_charges': 'true'
        }
        self.assertDictEqual(act_rem, test_set.rem)

    def test_double_solvation(self):
        raised_error = False
        dict_set = None
        try:
            dict_set = GSMDictSet(basis_set='6-31g*',
                                  scf_algorithm='diis',
                                  dft_rung=1,
                                  pcm_dielectric=10.0,
                                  smd_solvent="water",
                                  max_scf_cycles=35)
        except ValueError:
            raised_error = True

        self.assertTrue(raised_error)
        self.assertEqual(dict_set, None)

    def test_pcm_write(self):
        dict_set = GSMDictSet(basis_set='6-31g*',
                              scf_algorithm='diis',
                              dft_rung=5,
                              pcm_dielectric=10.0,
                              max_scf_cycles=35)
        dict_set.write("mol.qin")
        test_dict = QCTemplate.from_file("mol.qin").as_dict()
        rem = {
            "job_type": "force",
            "basis": "6-31G*",
            "max_scf_cycles": '35',
            "method": "wb97mv",
            "gen_scfman": 'true',
            "scf_algorithm": "diis",
            "xc_grid": '3',
            "solvent_method": "pcm",
            'symmetry': 'false',
            'sym_ignore': 'true',
            'resp_charges': 'true'
        }
        pcm = {
            "heavypoints": "194",
            "hpoints": "194",
            "radii": "uff",
            "theory": "cpcm",
            "vdwscale": "1.1"
        }
        qc_input = QCTemplate(rem=rem, pcm=pcm, solvent={"dielectric": "10.0"})
        for k, v in qc_input.as_dict().items():
            self.assertEqual(v, test_dict[k])
        os.remove("mol.qin")

    def test_smd_write(self):
        dict_set = GSMDictSet(basis_set='6-31g*',
                              scf_algorithm='diis',
                              dft_rung=5,
                              smd_solvent="water",
                              max_scf_cycles=35)
        dict_set.write("mol.qin")
        test_dict = QCTemplate.from_file("mol.qin").as_dict()
        rem = {
            "job_type": "force",
            "basis": "6-31G*",
            "max_scf_cycles": '35',
            "method": "wb97mv",
            "gen_scfman": 'true',
            "scf_algorithm": "diis",
            "xc_grid": '3',
            "solvent_method": "smd",
            "ideriv": "1",
            'symmetry': 'false',
            'sym_ignore': 'true',
            'resp_charges': 'true'
        }
        qc_input = QCTemplate(rem=rem, smx={"solvent": "water"})
        for k, v in qc_input.as_dict().items():
            self.assertEqual(v, test_dict[k])
        os.remove("mol.qin")

    def test_custom_smd_write(self):
        dict_set = GSMDictSet(basis_set='6-31g*',
                              scf_algorithm='diis',
                              dft_rung=5,
                              smd_solvent="custom",
                              custom_smd="90.00,1.415,0.00,0.735,20.2,0.00,0.00",
                              max_scf_cycles=35)
        dict_set.write("mol.qin")
        test_dict = QCTemplate.from_file("mol.qin").as_dict()
        rem = {
            "job_type": "force",
            "basis": "6-31G*",
            "max_scf_cycles": '35',
            "method": "wb97mv",
            "gen_scfman": 'true',
            "scf_algorithm": "diis",
            "xc_grid": '3',
            "solvent_method": "smd",
            "ideriv": "1",
            'symmetry': 'false',
            'sym_ignore': 'true',
            'resp_charges': 'true'
        }
        qc_input = QCTemplate(rem=rem, smx={"solvent": "other"})
        for k, v in qc_input.as_dict().items():
            self.assertEqual(v, test_dict[k])
        os.remove("mol.qin")
        with open("solvent_data") as sd:
            lines = sd.readlines()
            self.assertEqual(lines[0], "90.00,1.415,0.00,0.735,20.2,0.00,0.00")
        os.remove("solvent_data")

if __name__ == '__main__':
    unittest.main()
