# coding: utf-8

import unittest

import numpy as np

from pymatgen.core.structure import Molecule

from mpcat.utils.comparison import all_equal, compositions_equal


class ComparisonTest(unittest.TestCase):

    def test_all_equal(self):
        # Test with basic dict
        self.assertTrue(all_equal({"a": {"b": 1, "c": True}},
                                  {"a": {"b": 1, "c": True}}))
        self.assertFalse(all_equal({"a": {"b": 1, "c": True}},
                                    {"a": {"b": 2, "c": True}}))

        # Test with dict with list
        self.assertTrue(all_equal({"a": {"b": 1, "c": True},
                                   "d": [1, 2, 3, 4, 5]},
                                  {"a": {"b": 1, "c": True},
                                   "d": [1, 2, 3, 4, 5]}))
        self.assertFalse(all_equal({"a": {"b": 1, "c": True},
                                    "d": [1, 2, 3, 4, 5]},
                                   {"a": {"b": 1, "c": True},
                                    "d": [1, 2, 4, 5]}))
        self.assertFalse(all_equal({"a": {"b": 1, "c": True},
                                    "d": [1, 2, 3, 4, 5]},
                                   {"a": {"b": 1, "c": True},
                                    "d": [1, 2, 3, 4, 6]}))

        # Test with dict with np.array
        self.assertTrue(all_equal({"a": {"b": 1, "c": True},
                                   "d": [1, 2, 3, 4, 5],
                                   "e": np.array([["a", "b", "c"],
                                                  ["D", "E", "F"]])},
                                  {"a": {"b": 1, "c": True},
                                   "d": [1, 2, 3, 4, 5],
                                   "e": np.array([["a", "b", "c"],
                                                  ["D", "E", "F"]])}))
        self.assertTrue(all_equal({"a": {"b": 1, "c": True},
                                   "d": [1, 2, 3, 4, 5],
                                   "e": np.array([["a", "b", "c"],
                                                  ["D", "E", "F"]])},
                                  {"a": {"b": 1, "c": True},
                                   "d": [1, 2, 3, 4, 5],
                                   "e": [["a", "b", "c"],
                                         ["D", "E", "F"]]}))
        self.assertFalse(all_equal({"a": {"b": 1, "c": True},
                                    "d": [1, 2, 3, 4, 5],
                                    "e": np.array([["a", "b", "c"],
                                                   ["D", "E", "F"]])},
                                   {"a": {"b": 1, "c": True},
                                    "d": [1, 2, 3, 4, 5],
                                    "e": np.array([["a", "b", "c"],
                                                   ["f", "e", "d"]])}))
        self.assertFalse(all_equal({"a": {"b": 1, "c": True},
                                    "d": [1, 2, 3, 4, 5],
                                    "e": np.array([["a", "b", "c"],
                                                   ["D", "E", "F"]])},
                                   {"a": {"b": 1, "c": True},
                                    "d": [1, 2, 3, 4, 5],
                                    "e": [["a", "c"],
                                          ["f", "e"]]}))

        # Test different types
        self.assertFalse(all_equal({"a": {"b": 1, "c": True, "f": set()},
                                    "d": [1, 2, 3, 4, 5],
                                    "e": np.array([["a", "b", "c"],
                                                   ["D", "E", "F"]])},
                                   {"a": {"b": 1, "c": True, "f": "SET"},
                                    "d": [1, 2, 3, 4, 5],
                                    "e": [["a", "b", "c"],
                                          ["D", "E", "F"]]}))
        self.assertFalse(all_equal({"a": {"b": 1, "c": True, "f": set(),
                                          "g": None},
                                    "d": [1, 2, 3, 4, 5],
                                    "e": np.array([["a", "b", "c"],
                                                   ["D", "E", "F"]])},
                                   {"a": {"b": 1, "c": True, "f": set()},
                                    "d": [1, 2, 3, 4, 5],
                                    "e": [["a", "b", "c"],
                                          ["D", "E", "F"]]}))

    def test_compositions_equal(self):
        h = Molecule(["H"], [[0.0, 0.0, 0.0]])
        oh = Molecule(["O", "H"], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        h2o = Molecule(["O", "H", "H"], [[-0.013, -0.019, 0.0], [-0.299, 0.919, 0.0],
                                         [0.966, 0.024, 0.0]])
        h2 = Molecule(["H", "H"], [[0.0, 0.0, 0.0],
                                   [0.0, 0.0, 0.7]])
        o2 = Molecule(["O", "O"], [[0.0, 0.0, 0.0],
                                   [0.0, 0.0, 1.2]])

        self.assertTrue(compositions_equal([h, oh], [h2o]))
        self.assertFalse(compositions_equal([h2o], [h2, o2]))


if __name__ == "__main__":
    unittest.main()
