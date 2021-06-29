# coding: utf-8

import unittest
from pathlib import Path

from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN, metal_edge_extender

from mpcat.apprehend.autots_input import TSInput, TSSet
from mpcat.adapt.schrodinger_adapter import maestro_file_to_molecule


test_dir = Path(__file__).resolve().parent.parent.parent.parent / "test_files"


class TestTSInput(unittest.TestCase):
    def setUp(self) -> None:
        self.reactant_1 = maestro_file_to_molecule(
            test_dir / "autots_success_partial" / "rct1.mae")[0]

        self.reactant_2 = maestro_file_to_molecule(
            test_dir / "autots_success_partial" / "rct2.mae")[0]

        self.product = maestro_file_to_molecule(
            test_dir / "autots_success_partial" / "pro.mae")[0]

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
        autots_input = TSInput([self.reactant_1, self.reactant_2],
                                   [self.product],
                                   self.autots_variables,
                                   self.gen_variables)

        self.assertEqual(
            metal_edge_extender(
                MoleculeGraph.with_local_env_strategy(
                    self.reactant_1, OpenBabelNN()
                )
            ),
            autots_input.reactants[0])
        self.assertEqual(
            metal_edge_extender(
                MoleculeGraph.with_local_env_strategy(
                    self.reactant_2, OpenBabelNN()
                )
            ),
            autots_input.reactants[1])
        self.assertEqual(
            metal_edge_extender(
                MoleculeGraph.with_local_env_strategy(
                    self.product, OpenBabelNN()
                )
            ),
            autots_input.products[0])
        self.assertDictEqual(self.autots_variables, autots_input.autots_variables)
        self.assertDictEqual(self.gen_variables, autots_input.gen_variables)

    def test_to_from_file(self):
        writing = TSInput([self.reactant_1, self.reactant_2],
                              [self.product],
                              self.autots_variables,
                              self.gen_variables)
        writing.write("test.in", write_molecules=True)

        self.assertEqual(self.reactant_1, maestro_file_to_molecule("rct_0.mae")[0])
        self.assertEqual(self.reactant_2, maestro_file_to_molecule("rct_1.mae")[0])
        self.assertEqual(self.product, maestro_file_to_molecule("pro_0.mae")[0])

        reading = TSInput.from_file("test.in", read_molecules=True)
        self.assertEqual(reading.reactants[0], writing.reactants[0])
        self.assertEqual(reading.reactants[1], writing.reactants[1])
        self.assertEqual(reading.products[0], writing.products[0])

        for k, v in writing.autots_variables.items():
            self.assertEqual(v, reading.autots_variables[k])
        for k, v in writing.gen_variables.items():
            self.assertEqual(v, reading.gen_variables[k])


class TestTSSet(unittest.TestCase):

    def setUp(self) -> None:
        self.reactant_1 = maestro_file_to_molecule(
            test_dir / "autots_success_partial" / "rct1.mae")[0]

        self.reactant_2 = maestro_file_to_molecule(
            test_dir / "autots_success_partial" / "rct2.mae")[0]

        self.product = maestro_file_to_molecule(
            test_dir / "autots_success_partial" / "pro.mae")[0]

        self.autots_variables = {"eliminate_multiple_frequencies": True,
                                 "free_energy": True,
                                 "require_irc_success": True,
                                 "ts_vet_max_freq": -40.0,
                                 "units": "ev",
                                 "use_alternate": True}

        self.gen_variables = {"dftname": "wb97x-d",
                              "basis": "def2-tzvpd",
                              "babel": "xyz",
                              "ip472": 2,  # Output all steps of geometry optimization in *.mae
                              # "ip172": 2,  # Print RESP file
                              "ip175": 2,  # Print XYZ files
                              # "ifreq": 1,  # Frequency calculation
                              # "irder": 1,  # IR vibrational modes calculated
                              # "nmder": 2,  # Numerical second derivatives
                              "nogas": 2,  # Skip gas-phase optimization, if PCM is used
                              "maxitg": 250,  # Maximum number of geometry optimization iterations
                              # "intopt_switch": 0,  # Do not switch from internal to Cartesian coordinates
                              # "optcoord_update": 0,  # Do not run checks to change coordinate system
                              # "props_each_step": 1,  # Calculate properties at each optimization step
                              # "iaccg": 5  # Tight convergence criteria for optimization
                              # "mulken": 1,  # Calculate Mulliken properties by atom
                              "maxit": 400,  # Maximum number of SCF iterations
                              "iacc": 2,  # Use "accurate" SCF convergence criteria
                              # "noauto": 3  # All calculations done on fine grid
                              "isymm": 0,  # Do not use symmetry
                              # "espunit": 6  # Electrostatic potential in units of eV
                              }

    def test_defaults(self):
        default_set = TSSet([self.reactant_1, self.reactant_2],
                                [self.product])

        for key, value in self.gen_variables.items():
            self.assertEqual(value, default_set.gen_variables[key])
        for key, value in self.autots_variables.items():
            self.assertEqual(value, default_set.autots_variables[key])

    def test_not_defaults(self):
        non_default_set = TSSet([self.reactant_1, self.reactant_2],
                                [self.product],
                                basis_set="6-31+g(d)",
                                dft_rung=1,
                                pcm_settings={"solvent": "benzene",
                                              "model": "cosmo"},
                                max_scf_cycles=500,
                                geom_opt_max_cycles=500,
                                overwrite_inputs_autots={"ts_vet_max_freq": -10.0},
                                overwrite_inputs_gen={"noauto": 3})

        self.assertEqual(non_default_set.gen_variables["basis"], "6-31+g(d)")
        self.assertEqual(non_default_set.gen_variables["dftname"], "hfs")
        self.assertEqual(non_default_set.gen_variables["maxit"], 500)
        self.assertEqual(non_default_set.gen_variables["maxitg"], 500)
        self.assertEqual(non_default_set.gen_variables["noauto"], 3)
        self.assertEqual(non_default_set.gen_variables["isolv"], 7)
        self.assertEqual(non_default_set.gen_variables["solvent"], "benzene")
        self.assertEqual(non_default_set.gen_variables["pcm_model"], "cosmo")

        self.assertEqual(non_default_set.autots_variables["ts_vet_max_freq"], -10.0)


if __name__ == "__main__":
    unittest.main()
