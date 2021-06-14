# coding: utf-8
# Distributed under the terms of the MIT License.

import os
import unittest

from monty.serialization import loadfn

from pymatgen.core.structure import Molecule
from pymatgen.util.testing import PymatgenTest

from mpcat.automate.atomate.outputs import (GSMOptimizedStringParser,
                                            GSMInternalCoordinateDataParser)


__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2020 The Materials Project"
__version__ = "0.1"
__email__ = "espottesmith@gmail.com"
__maintainer__ = "Evan Spotte-Smith"

gsm_dir = os.path.join(os.path.dirname(__file__), "test_files")


class TestGSMOutput(PymatgenTest):
    # Going to hold off here until I have more data points
    # After G2 run, should be possible
    pass


class TestGSMOptimizedStringParser(PymatgenTest):
    def test_compare_with_dict(self):
        self.maxDiff = None

        reference = loadfn(os.path.join(gsm_dir, "optimized_string.json"))
        compare = GSMOptimizedStringParser(os.path.join(gsm_dir, "opt_converged_000_000.xyz")).as_dict()

        self.assertEqual(compare["text"], reference["text"])
        self.assertSequenceEqual(compare["lines"], reference["lines"])
        for key in compare["data"]:
            if key == "molecules":
                for i, e in enumerate(compare["data"]["molecules"]):
                    self.assertEqual(Molecule.from_dict(e), reference["data"]["molecules"][i])
            else:
                try:
                    self.assertEqual(compare["data"][key], reference["data"][key])
                except ValueError:
                    self.assertSequenceEqual(compare["data"][key], reference["data"][key])


class TestGSMInternalCoordinateDataParser(PymatgenTest):
    def test_compare_with_dict(self):
        self.maxDiff = None

        reference = loadfn(os.path.join(gsm_dir, "internals.json"))
        compare = GSMInternalCoordinateDataParser(os.path.join(gsm_dir, "IC_data_0000.txt")).as_dict()

        self.assertEqual(compare["text"], reference["text"])

        self.assertEqual(compare["data"]["rct_node"], reference["data"]["rct_node"])
        self.assertEqual(compare["data"]["ts_node"], reference["data"]["ts_node"])
        self.assertEqual(compare["data"]["pro_node"], reference["data"]["pro_node"])

        for key in compare["data"]["rct"]:
            self.assertSequenceEqual(compare["data"]["rct"][key],
                                     reference["data"]["rct"][key])
        for key in compare["data"]["ts"]:
            self.assertSequenceEqual(compare["data"]["ts"][key],
                                     reference["data"]["ts"][key])
        for key in compare["data"]["pro"]:
            self.assertSequenceEqual(compare["data"]["pro"][key],
                                     reference["data"]["pro"][key])



if __name__ == "__main__":
    unittest.main()
