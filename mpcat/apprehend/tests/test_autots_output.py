import os
import unittest

from monty.serialization import loadfn, dumpfn

from mpcat.apprehend.autots_output import AutoTSOutput
from mpcat.utils.comparison import all_equal

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files")


def prepare_files():
    successful_file = AutoTSOutput(os.path.join(test_dir, "autots_success_partial",
                                                "AutoTS.T9XnCsLi.out"))
    failed_file = AutoTSOutput(os.path.join(test_dir, "autots_failure_partial",
                                            "AutoTS.T9XnCsLi.out"))

    dumpfn(successful_file, os.path.join(test_dir, "autots_success_partial", "success_ts.json"))
    dumpfn(failed_file, os.path.join(test_dir, "autots_failure_partial", "failure_ts.json"))


class TestJagOutput(unittest.TestCase):

    def setUp(self) -> None:
        self.successful_file = loadfn(os.path.join(test_dir, "autots_success_partial", "success_ts.json"))
        self.failed_file = loadfn(os.path.join(test_dir, "autots_failure_partial", "failure_ts.json"))

    def test_success(self):
        successful_file = AutoTSOutput(os.path.join(test_dir, "autots_success_partial",
                                                    "AutoTS.T9XnCsLi.out"))

        self.assertTrue(all_equal(successful_file.data, self.successful_file.data))

    def test_failure(self):
        failed_file = AutoTSOutput(os.path.join(test_dir, "autots_failure_partial",
                                                "AutoTS.T9XnCsLi.out"))

        self.assertTrue(all_equal(failed_file.data, self.failed_file.data))


if __name__ == "__main__":
    unittest.main()
