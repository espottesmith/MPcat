import os
import unittest

from mpcat.aggregate.drones import AutoTSDrone
from monty.serialization import loadfn


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files")


class TestAutoTSDrone(unittest.TestCase):

    def test_compare_to_reference(self):
        drone = AutoTSDrone(path=os.path.join(test_dir, "decomp_ro1"))
        doc = drone.assimilate()

        reference = loadfn(os.path.join(test_dir, "decomp_ro1_doc.json"))

        self.assertEqual(doc, reference)


if __name__ == "__main__":
    unittest.main()