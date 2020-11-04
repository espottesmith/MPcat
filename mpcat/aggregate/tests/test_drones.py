# coding: utf-8

import unittest
from pathlib import Path

from mpcat.aggregate.drones import AutoTSCalcDrone
from monty.serialization import loadfn, dumpfn


test_dir = Path(__file__).resolve().parent.parent.parent.parent / "test_files"


class TestAutoTSDronePMG(unittest.TestCase):

    def test_compare_to_reference(self):
        drone = AutoTSCalcDrone(path=test_dir / "decomp_ro1")
        doc = drone.assimilate()
        dumpfn(doc, "test.json")
        doc = loadfn("test.json")

        reference = loadfn((test_dir / "decomp_ro1_doc.json").as_posix())

        self.assertEqual(doc["schema"], reference["schema"])
        self.assertEqual(doc["input"], reference["input"])
        self.assertEqual(doc["output"], reference["output"])
        self.assertEqual(doc["calcs"], reference["calcs"])
        self.assertEqual(doc["completed"], reference["completed"])
        self.assertEqual(doc["calculation_names"], reference["calculation_names"])

        Path("test.json").unlink()


if __name__ == "__main__":
    unittest.main()
