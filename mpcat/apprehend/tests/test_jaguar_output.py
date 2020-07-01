import os
import unittest

from monty.serialization import loadfn, dumpfn

from mpcat.apprehend.jaguar_output import JagOutput
from mpcat.utils.comparison import all_equal

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files")


def prepare_files():
    successful_file = JagOutput(os.path.join(test_dir, "autots_success_partial",
                                             "AutoTS.T9XnCsLi_opt_0.out"))
    failed_file = JagOutput(os.path.join(test_dir, "autots_failure_partial",
                                         "AutoTS.T9XnCsLi_opt_0.out"))

    dumpfn(successful_file, os.path.join(test_dir, "autots_success_partial", "success_jag.json"))
    dumpfn(failed_file, os.path.join(test_dir, "autots_failure_partial", "failure_jag.json"))


class TestJagOutput(unittest.TestCase):

    def setUp(self) -> None:
        self.successful_file = loadfn(os.path.join(test_dir, "autots_success_partial", "success_jag.json"))
        self.failed_file = loadfn(os.path.join(test_dir, "autots_failure_partial", "failure_jag.json"))

    def test_success(self):
        successful_file = JagOutput(os.path.join(test_dir, "autots_success_partial",
                                                 "AutoTS.T9XnCsLi_opt_0.out"))

        self.assertTrue(all_equal(successful_file.data, self.successful_file.data))

    def test_failure(self):
        failed_file = JagOutput(os.path.join(test_dir, "autots_failure_partial",
                                             "AutoTS.T9XnCsLi_opt_0.out"))

        self.assertTrue(all_equal(failed_file.data, self.failed_file.data))
