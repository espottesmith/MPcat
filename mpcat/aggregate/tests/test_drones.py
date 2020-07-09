# coding: utf-8

import os
import unittest

from mpcat.aggregate.drones import AutoTSCalcDrone, AutoTSBuilderDrone
from monty.serialization import loadfn, dumpfn


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files")


class TestAutoTSDronePMG(unittest.TestCase):

    def test_compare_to_reference(self):
        drone = AutoTSCalcDrone(path=os.path.join(test_dir, "decomp_ro1"))
        doc = drone.assimilate()
        dumpfn(doc, "test.json")
        doc = loadfn("test.json")

        reference = loadfn(os.path.join(test_dir, "decomp_ro1_doc.json"))

        self.assertEqual(doc["schema"], reference["schema"])
        self.assertEqual(doc["input"], reference["input"])
        self.assertEqual(doc["output"], reference["output"])
        self.assertEqual(doc["calcs"], reference["calcs"])
        self.assertEqual(doc["completed"], reference["completed"])
        self.assertEqual(doc["calculation_names"], reference["calculation_names"])

        os.remove("test.json")


if __name__ == "__main__":
    unittest.main()