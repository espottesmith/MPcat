import os
import unittest

from mpcat.apprehend.autots_input import AutoTSInput, AutoTSSet
from mpcat.adapt.schrodinger_adapter import maestro_file_to_molecule


test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files")


class TestAutoTSInput(unittest.TestCase):
    def setUp(self) -> None:
        self.reactant_1 = maestro_file_to_molecule(os.path.join(test_dir,
                                                                "autots_success_partial",
                                                                "rct1.mae"))

        self.reactant_2 = maestro_file_to_molecule(os.path.join(test_dir,
                                                                "autots_success_partial",
                                                                "rct2.mae"))

        self.product = maestro_file_to_molecule(os.path.join(test_dir,
                                                             "autots_success_partial",
                                                             "pro.mae"))

        self.autots_variables = {"eliminate_multiple_frequencies": True,
                                 "free_energy": True,
                                 "require_irc_success": True,
                                 "ts_vet_max_freq": -40.0,
                                 "units": "ev",
                                 "use_alternate": True}

        self.gen_variables = {"dftname": "wb97X-D",
                              "basis": "def2-tzvpd",
                              "babel": "xyz",
                              "ip472": 2,  # Output all steps of geometry optimization in *.mae
                              "ip172": 2,  # Print RESP file
                              "ip175": 2,  # Print XYZ files
                              "ifreq": 1,  # Frequency calculation
                              "irder": 1,  # IR vibrational modes calculated
                              "nmder": 2,  # Numerical second derivatives
                              "nogas": 2,  # Skip gas-phase optimization, if PCM is used
                              "maxitg": 250,  # Maximum number of geometry optimization iterations
                              "intopt_switch": 0,  # Do not switch from internal to Cartesian coordinates
                              "optcoord_update": 0,  # Do not run checks to change coordinate system
                              "props_each_step": 1,  # Calculate properties at each optimization step
                              "mulken": 1,  # Calculate Mulliken properties by atom
                              "maxit": 250,  # Maximum number of SCF iterations
                              "iacc": 2,  # Use "accurate" SCF convergence criteria
                              "isymm": 0,  # Do not use symmetry
                              "espunit": 6  # Electrostatic potential in units of eV
                              }

    def test_generate(self):
        autots_input = AutoTSInput([self.reactant_1, self.reactant_2],
                                   [self.product],
                                   self.autots_variables,
                                   self.gen_variables)

        self.assertEqual(self.reactant_1, autots_input.reactants[0])
        self.assertEqual(self.reactant_2, autots_input.reactants[1])
        self.assertEqual(self.product, autots_input.products[0])
        self.assertDictEqual(self.autots_variables, autots_input.autots_variables)
        self.assertDictEqual(self.gen_variables, autots_input.gen_variables)

    def test_to_from_file(self):
        writing = AutoTSInput([self.reactant_1, self.reactant_2],
                              [self.product],
                              self.autots_variables,
                              self.gen_variables)
        writing.write("test.in", write_molecules=True)

        self.assertEqual(self.reactant_1, maestro_file_to_molecule("rct_0.mae")[0])
        self.assertEqual(self.reactant_2, maestro_file_to_molecule("rct_1.mae")[0])
        self.assertEqual(self.product, maestro_file_to_molecule("pro_0.mae")[0])

        reading = AutoTSInput.from_file("test.in", read_molecules=True)
        self.assertEqual(reading.reactants[0], writing.reactants[0])
        self.assertEqual(reading.reactants[1], writing.reactants[1])
        self.assertEqual(reading.products[0], writing.products[0])

        for k, v in writing.autots_variables.items():
            self.assertEqual(v, reading.autots_variables[k])
        for k, v in writing.gen_variables.items():
            self.assertEqual(v, reading.gen_variables[k])


class TestAutoTSSet(unittest.TestCase):

    def test_defaults(self):
        pass

    def test_not_defaults(self):
        pass


if __name__ == "__main__":
    unittest.main()
