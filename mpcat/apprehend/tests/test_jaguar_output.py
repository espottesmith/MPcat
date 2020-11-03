# coding: utf-8

import unittest
from pathlib import Path

from monty.serialization import loadfn, dumpfn

from mpcat.apprehend.jaguar_output import JagOutput
from mpcat.utils.comparison import all_equal


test_dir = Path(__file__).resolve().parent.parent.parent / "test_files"


def prepare_files():
    successful_file = JagOutput(
        test_dir / "autots_success_partial" / "AutoTS.T9XnCsLi_opt_0.out",
        parse_molecules=False)
    failed_file = JagOutput(
        test_dir / "autots_failure_partial" / "AutoTS.T9XnCsLi_opt_0.out",
        parse_molecules=False)

    dumpfn(successful_file.as_dict(),
           test_dir / "autots_success_partial" / "success_jag.json")
    dumpfn(failed_file.as_dicct(),
           test_dir / "autots_failure_partial" / "failure_jag.json")


class TestJagOutput(unittest.TestCase):

    def setUp(self) -> None:
        self.successful_file = JagOutput.from_dict(loadfn(
            test_dir / "autots_success_partial" / "success_jag.json"))
        self.failed_file = JagOutput.from_dict(loadfn(
            test_dir / "autots_failure_partial" / "failure_jag.json"))

    def test_success(self):
        successful_file = JagOutput(
            test_dir / "autots_success_partial" / "AutoTS.T9XnCsLi_opt_0.out",
            parse_molecules=False)

        self.assertTrue(all_equal(successful_file.data, self.successful_file.data,
                                  exclude_keys=["vibrational_frequency_modes"]))

    def test_failure(self):
        failed_file = JagOutput(
            test_dir / "autots_failure_partial" / "AutoTS.T9XnCsLi_opt_0.out",
            parse_molecules=False)

        self.assertTrue(all_equal(failed_file.data, self.failed_file.data,
                                  exclude_keys=["vibrational_frequency_modes"]))


if __name__ == "__main__":
    unittest.main()
