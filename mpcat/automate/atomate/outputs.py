# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

import logging

import numpy as np

from monty.io import zopen
from monty.json import MSONable, jsanitize

from pymatgen.core.structure import Molecule

from pymatgen.io.qchem.utils import (read_pattern,
                                     read_table_pattern)

# Classes for reading/manipulating/writing output files for use with pyGSM.

__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2021"
__version__ = "0.1"
__email__ = "espottesmith@gmail.com"
__credit__ = "Sam Blau"

logger = logging.getLogger(__name__)


class GSMOutput(MSONable):

    """
    An object that parses and contains data from pyGSM output files.

    TODO: Add comments

    Args:
        filename (str): The path to the pyGSM output file.
    """

    def __init__(self, filename):
        self.filename = filename
        self.data = dict()
        self.data["errors"] = list()
        self.data["warnings"] = list()
        self.text = ""
        with zopen(filename, "rt") as f:
            self.text = f.read()

        completion = read_pattern(self.text, {"key": r"Finished GSM!"},
                                  terminate_on_match=True).get("key")
        if completion is None:
            self.data["completion"] = False
        else:
            self.data["completion"] = True

        self._parse_input_values()
        self._parse_initial_energies()
        self._parse_node_opt_trajectory()

        if self.data["inputs"]["gsm_type"] == "SE_GSM":
            self._parse_driving_coordinates()
            self._parse_coordinate_trajectory()

        self._parse_opt_summary()

        if self.data["completion"]:
            self._parse_summary_info()

        self._parse_warnings()
        self._parse_errors()

    def _parse_input_values(self):
        header_pattern = r"#=+#\n#\|.+\[92m\s+Parsed GSM Keys : Values.+\[0m\s+\|#\n#=+#"
        table_pattern = r"(?P<key>[A-Za-z_]+)\s+(?P<value>[A-Za-z0-9\[\]_\.\-]+)\s*\n"
        footer_pattern = r"\-+"

        temp_inputs = read_table_pattern(self.text, header_pattern,
                                         table_pattern, footer_pattern)

        if temp_inputs is None or len(temp_inputs) == 0:
            self.data["inputs"] = dict()
        else:
            self.data["inputs"] = dict()
            temp_inputs = temp_inputs[0]

            for row in temp_inputs:
                key = row["key"]
                value = row["value"]

                if value == "True":  # Deal with bools
                    self.data["inputs"][key] = True
                elif value == "False":
                    self.data["inputs"][key] = False
                elif value.startswith("[") and value.endswith("]"):  # Deal with lists
                    val = value[1:-1].split(", ")
                    try:  # ints
                        val = [int(v) for v in val]
                        self.data["inputs"][key] = val
                    except ValueError:
                        self.data["inputs"][key] = val
                else:
                    # int
                    is_int = True
                    is_float = True
                    val = value
                    try:
                        val = int(value)
                    except ValueError:
                        is_int = False
                    if is_int:
                        self.data["inputs"][key] = val
                        continue
                    else:
                        try:
                            val = float(value)
                        except ValueError:
                            is_float = False

                    if is_float:
                        self.data["inputs"][key] = val
                        continue
                    else:
                        self.data["inputs"][key] = value

        if "charge" not in self.data["inputs"]:
            self.data["inputs"]["charge"] = 0

    def _parse_initial_energies(self):

        if self.data["inputs"]["gsm_type"] == "SE_GSM":
            temp_init_energy = read_pattern(self.text,
                                            {"single": r"\s*Initial energy is ([\-\.0-9]+)"},
                                            terminate_on_match=True).get("single")
            self.data["initial_energy"] = float(temp_init_energy[0][0])

        elif self.data["inputs"]["gsm_type"] == "DE_GSM":
            temp_init_energy = read_pattern(self.text, {
                "double": r"\s*Energy of the end points are ([\-\.0-9]+), ([\-\.0-9]+)",
                "double_relative": r"\s*relative E ([\-\.0-9]+), ([\-\.0-9]+)"},
                                            terminate_on_match=True)
            if temp_init_energy.get("double"):
                self.data["initial_energy_rct"] = float(temp_init_energy.get("double")[0][0])
                self.data["initial_energy_pro"] = float(temp_init_energy.get("double")[0][1])
                self.data["initial_relative_energy_rct"] = float(temp_init_energy.get("double_relative")[0][0])
                self.data["initial_relative_energy_pro"] = float(temp_init_energy.get("double_relative")[0][1])

    def _parse_driving_coordinates(self):
        temp_coords = read_pattern(self.text, {
            "key": r"driving coordinates \[((\['(ADD|BREAK|ANGLE|TORSION|OOP)', ([0-9]+,? ?)+\],? ?)+)\]"
        }, terminate_on_match=True).get("key")

        self.data["driving_coords"] = dict()
        self.data["driving_coords"]["add"] = list()
        self.data["driving_coords"]["break"] = list()
        self.data["driving_coords"]["angle"] = list()
        self.data["driving_coords"]["torsion"] = list()
        self.data["driving_coords"]["out_of_plane"] = list()

        coord_sets = temp_coords[0][0].split("], [")
        for coord_set in coord_sets:
            tokens = coord_set.strip("[]").split(", ")
            self.data["driving_coords"][tokens[0].strip("'").lower()].append(tuple([int(e) for e in tokens[1:]]))

    def _parse_coordinate_trajectory(self):
        #TODO: Hack pyGSM so that all of these return the indices associated with the coordinate
        #TODO: Also use consistent units?

        temp_coord_trajectory = read_pattern(self.text, {
            "add": r"\s*bond \(([0-9]+), ([0-9]+)\) target \(greater than\): ([\.0-9]+), current d: ([\.0-9]+) diff: [\-\.0-9]+",
            "break": r"s*bond \(([0-9]+), ([0-9]+)\) target \(greater than\): ([\.0-9]+), current d: ([\.0-9]+) diff: [\-\.0-9]+",
            "angle": r"s*anglev: ([\-\.0-9]+) align to ([\-\.0-9]+) diff(rad): [\-\.0-9]+",
            "torsion": r"s*current torv: ([\-\.0-9]+) align to ([\-\.0-9]+) diff(deg): [\-\.0-9]+"
        })

        # Initialize dictionary for driving coordinate trajectories
        self.data["driving_coord_trajectories"] = {"add": dict(),
                                                   "break": dict(),
                                                   "angle": dict(),
                                                   "torsion": dict()}

        self.data["driving_coord_goals"] = {"add": dict(),
                                            "break": dict(),
                                            "angle": dict(),
                                            "torsion": dict()}

        for coord_type, coords in self.data["driving_coords"].items():
            for coord in coords:
                self.data["driving_coord_trajectories"][coord_type][coord] = list()

        for add_coord in temp_coord_trajectory.get("add", list()):
            bond = (int(add_coord[0]), int(add_coord[1]))
            if bond in self.data["driving_coords"]["add"]:
                if bond not in self.data["driving_coord_goals"]["add"]:
                    self.data["driving_coord_goals"]["add"][bond] = float(add_coord[2])
                self.data["driving_coord_trajectories"]["add"][bond].append(float(add_coord[3]))

        for break_coord in temp_coord_trajectory.get("break", list()):
            bond = (int(break_coord[0]), int(break_coord[1]))
            if bond in self.data["driving_coords"]["break"]:
                if bond not in self.data["driving_coord_goals"]["break"]:
                    self.data["driving_coord_goals"]["break"][bond] = float(break_coord[2])
                self.data["driving_coord_trajectories"]["break"][bond].append(float(break_coord[3]))

        for e, ang_coord in enumerate(temp_coord_trajectory.get("angle", list())):
            #TODO: Fix this hack once the indices are printed for angles
            angle_ind = e % len(self.data["driving_coords"]["angle"])
            angle = self.data["driving_coords"]["angle"][angle_ind]
            if angle not in self.data["driving_coord_goals"]["angle"]:
                self.data["driving_coord_goals"]["angle"][angle] = float(ang_coord[1])
            self.data["driving_coord_trajectories"]["angle"][angle].append(float(ang_coord[0]))

        for e, tors_coord in enumerate(temp_coord_trajectory.get("torsion", list())):
            #TODO: Fix this hack once the indices are printed for angles
            tors_ind = e % len(self.data["driving_coords"]["torsion"])
            tors = self.data["driving_coords"]["torsion"][tors_ind]
            if tors not in self.data["driving_coord_goals"]["torsion"]:
                self.data["driving_coord_goals"]["torsion"][tors] = float(tors_coord[1])
            self.data["driving_coord_trajectories"]["torsion"][tors].append(float(tors_coord[0]))

    def _parse_node_opt_trajectory(self):
        self.data["opt_trajectory_energies"] = dict()
        self.data["opt_trajectory_gradrms"] = dict()

        header_pattern = r"\s*converged\n\s*opt-summary [0-9]+"
        body_pattern = r"\s*Node: ([0-9]+) Opt step: [0-9]+ E: ([\.\-0-9]+) predE: [\.\-0-9]+ ratio: [\.\-0-9]+ gradrms: ([\.0-9]+) ss: [\-\.0-9]+ DMAX: ([\.0-9]+)"
        footer_pattern = r""
        temp_opt_trajectories = read_table_pattern(self.text,
                                                   header_pattern=header_pattern,
                                                   row_pattern=body_pattern,
                                                   footer_pattern=footer_pattern)

        for table in temp_opt_trajectories:
            energies = list()
            grads = list()
            node_num = int(table[0][0])
            if node_num not in self.data["opt_trajectory_energies"]:
                self.data["opt_trajectory_energies"][node_num] = list()
            if node_num not in self.data["opt_trajectory_gradrms"]:
                self.data["opt_trajectory_gradrms"][node_num] = list()
            for row in table:
                energies.append(float(row[1]))
                grads.append(float(row[2]))

            self.data["opt_trajectory_energies"][node_num].append(energies)
            self.data["opt_trajectory_gradrms"][node_num].append(grads)

    def _parse_opt_summary(self):
        temp_opt_summary = read_pattern(self.text, {
            "v_profile": r"\s*V_profile:((\s+[\.\-0-9]+)+)",
            "v_profile_re": r"s*V_profile \(after reparam\):((\s+[\.\-0-9]+)+)",
            "all_uphill": r"\s*all uphill\?\s+(True|False)",
            "dissociative": r"\s*dissociative\?\s+(True|False)",
            "min_nodes": r"\s*min nodes\s+\[(([0-9]+,? ?)+)\]",
            "max_nodes": r"\s*max nodes\s+\[(([0-9]+,? ?)+)\]",
            "emax_nmax": r"\s*emax and nmax in find peaks ([\.\-0-9]+),([0-9]+)"
        })

        self.data["energy_profiles"] = list()
        self.data["reparameterized_energy_profiles"] = list()
        self.data["path_uphill"] = list()
        self.data["path_dissociative"] = list()
        self.data["path_min_nodes"] = list()
        self.data["path_max_nodes"] = list()
        self.data["max_nodes"] = list()
        self.data["max_energies"] = list()

        for v_profile in temp_opt_summary.get("v_profile", []):
            this_profile = list()
            for entry in v_profile[0].split(" "):
                try:
                    this_profile.append(float(entry))
                except ValueError:
                    continue
            self.data["energy_profiles"].append(this_profile)

        for v_profile in temp_opt_summary.get("v_profile_re", []):
            this_profile = list()
            for entry in v_profile[0].split(" "):
                try:
                    this_profile.append(float(entry))
                except ValueError:
                    continue
            self.data["reparameterized_energy_profiles"].append(this_profile)

        for uphill in temp_opt_summary.get("all_uphill", []):
            if uphill[0] == "True":
                self.data["path_uphill"].append(True)
            elif uphill[0] == "False":
                self.data["path_uphill"].append(False)
            else:
                # Don't know how this would happen
                self.data["path_uphill"].append(None)

        for dissoc in temp_opt_summary.get("dissociative", []):
            if dissoc[0] == "True":
                self.data["path_dissociative"].append(True)
            elif dissoc[0] == "False":
                self.data["path_dissociative"].append(False)
            else:
                # Don't know how this would happen
                self.data["path_dissociative"].append(None)

        for min_node_set in temp_opt_summary.get("min_nodes", []):
            nodes = [int(n) for n in min_node_set[0].split(", ")]
            self.data["path_min_nodes"].append(nodes)

        for max_node_set in temp_opt_summary.get("max_nodes", []):
            nodes = [int(n) for n in max_node_set[0].split(", ")]
            self.data["path_max_nodes"].append(nodes)

        if len(self.data["path_min_nodes"]) == 0 and len(self.data["path_max_nodes"]) != 0:
            for i in range(len(self.data["path_max_nodes"])):
                self.data["path_min_nodes"].append([])

        if len(self.data["path_max_nodes"]) == 0 and len(self.data["path_min_nodes"]) != 0:
            for i in range(len(self.data["path_min_nodes"])):
                self.data["path_max_nodes"].append([])

        for maxes in temp_opt_summary.get("emax_nmax", []):
            self.data["max_nodes"].append(int(maxes[1]))
            self.data["max_energies"].append(float(maxes[0]))

        self.data["final_energy_profile"] = self.data["energy_profiles"][-1]
        self.data["final_path_uphill"] = self.data["path_uphill"][-1]
        self.data["final_path_dissociative"] = self.data["path_dissociative"][-1]
        self.data["final_min_nodes"] = self.data["path_min_nodes"][-1]
        self.data["final_max_nodes"] = self.data["path_max_nodes"][-1]
        self.data["final_max_node"] = self.data["max_nodes"][-1]
        self.data["final_max_energy"] = self.data["max_energies"][-1]

    def _parse_summary_info(self):
        temp_summary = read_pattern(self.text, {
            "ts_energy": r"\s*TS energy: ([\.\-0-9]+)",
            "absolute_ts_energy": r"\s*absolute energy TS node ([\.\-0-9]+)",
            "nodes": r"\s*min reactant node: ([0-9]+) min product node ([0-9]+) TS node is ([0-9]+)",
            "delta_e": r"\s*Delta E is ([\.\-0-9]+)"
        })

        temp_ts_energy = temp_summary.get("ts_energy")
        if temp_ts_energy is None:
            self.data["ts_energy"] = None
        else:
            self.data["ts_energy"] = float(temp_ts_energy[0][0])

        temp_absolute_ts_energy = temp_summary.get("absolute_ts_energy")
        if temp_absolute_ts_energy is None:
            self.data["absolute_ts_energy"] = None
        else:
            self.data["absolute_ts_energy"] = float(temp_absolute_ts_energy[0][0])

        temp_nodes = temp_summary.get("nodes")
        if temp_nodes is None:
            self.data["min_rct_node"] = None
            self.data["min_pro_node"] = None
            self.data["ts_node"] = None
        else:
            self.data["min_rct_node"] = int(temp_nodes[0][0])
            self.data["min_pro_node"] = int(temp_nodes[0][1])
            self.data["ts_node"] = int(temp_nodes[0][2])

        temp_delta_e = temp_summary.get("delta_e")
        if temp_delta_e is None:
            self.data["delta_e"] = None
        else:
            self.data["delta_e"] = float(temp_delta_e[0][0])

    def _parse_warnings(self):
        temp_warnings = read_pattern(self.text, {
            "charged_molecule": r"\s*Warning: charge is not implemented for all level of theories\. Make sure this is correct for your package\.",
            "out_of_plane_ordering": r"\s*Warning: OutOfPlane atoms are the same, ordering is different",
            "last_node_not_fully_optimized": r"\s*Warning last node still not optimized fully",
            "possible_memory_leak": r"Warning: more than 100 B-matrices stored, memory leaks likely",
            "svd_failure_perturb": r"SVD fails, perturbing coordinates and trying again",
            "no_networkx": r"NetworkX cannot be imported \(topology tools won't work\)\.  Most functionality should still work though\.",
            "too_many_atoms": r"Warning: Large number of atoms [0-9]+, topology building may take a long time",
            "topology_broken_mol": r"Warning: Topology building will not work with broken molecules in nonorthogonal cells\.",
            "failed_constraint": r"Warning: Failed to apply Constraint",
            "redundant_constraint": r"Constraint [0-9]+ is almost redundant; after projection norm is [\.\-0-9]+ of original",
            "cartesian_different_weights": r"Warning: Cartesian[XYZ] same atoms, different weights \([0-9\.\-]+ [0-9\.\-]+\)",
            "translation_different_weights": r"Warning: Translation[XYZ] same atoms, different weights",
            "rotator_different_reference_positions": r"Warning: Rotator same atoms, different reference positions"
        })

        if temp_warnings.get("charged_molecule"):
            self.data["warnings"].append("charged_molecule")
        if temp_warnings.get("out_of_plane_ordering"):
            self.data["warnings"].append("out_of_plane_ordering")
        if temp_warnings.get("last_node_not_fully_optimized"):
            self.data["warnings"].append("last_node_not_fully_optimized")
        if temp_warnings.get("possible_memory_leak"):
            self.data["warnings"].append("possible_memory_leak")
        if temp_warnings.get("svd_failure_perturb"):
            self.data["warnings"].append("svd_failure_perturb")
        if temp_warnings.get("no_networkx"):
            self.data["warnings"].append("no_networkx")
        if temp_warnings.get("too_many_atoms"):
            self.data["warnings"].append("too_many_atoms")
        if temp_warnings.get("topology_broken_mol"):
            self.data["warnings"].append("topology_broken_mol")
        if temp_warnings.get("failed_constraint"):
            self.data["warnings"].append("failed_constraint")
        if temp_warnings.get("redundant_constraint"):
            self.data["warnings"].append("redundant_constraint")
        if temp_warnings.get("cartesian_different_weights"):
            self.data["warnings"].append("cartesian_different_weights")
        if temp_warnings.get("translation_different_weights"):
            self.data["warnings"].append("translation_different_weights")
        if temp_warnings.get("rotator_different_reference_positions"):
            self.data["warnings"].append("rotator_different_reference_positions")

    def _parse_errors(self):
        # Q-Chem somehow failed or didn't produce a gradient file
        # Requires analysis of the Q-Chem output file
        if read_pattern(self.text, {
            "key": r"FileNotFoundError: \[Errno 2\] No such file or directory: '.*GRAD'"
        }, terminate_on_match=True).get("key") == [[]]:
            self.data["errors"].append("qchem_failure_nograd")

        # Not sure why this error arises
        # It has something to do with the constraints formed
        # Might have been fixed by recent work; keeping it for now
        if read_pattern(self.text, {
            "key": r"if \(C\[:\]==0\.\)\.all\(\):\nAttributeError: 'bool' object has no attribute 'all'"
        }, terminate_on_match=True).get("key") == [[]]:
            self.data["errors"].append("coord_constraint_bool")

        # This is a bug that has been identified by the lead developer
        # Might have been fixed in the most recent release
        if read_pattern(self.text, {
            "key": r"NameError: name 'fp' is not defined"
        }, terminate_on_match=True).get("key") == [[]]:
            self.data["errors"].append("product_geom_not_fixed_bug")

    def as_dict(self):
        d = dict()
        d["data"] = self.data
        d["text"] = self.text
        d["filename"] = self.filename
        return jsanitize(d, strict=True)


class GSMOptimizedStringParser(MSONable):

    """
    An object that parses and contains data from pyGSM optimized geometry files.

    NOTE: This class is probably even less elegant than the usual Output classes
        (QCOutput, etc.) because the formatting of the optimized string files
        does not lend itself well to table-based or key-based parsing.

    Args:
        filename (str): The path to the pyGSM optimized geometry file.
    """

    def __init__(self, filename):
        self.filename = filename
        self.data = dict()
        self.text = ""
        with zopen(self.filename, "rt") as f:
            self.text = f.read()

        self.lines = self.text.split("\n")
        self.data["num_atoms"] = int(self.lines[2])

        self._parse_xyzs()
        self._parse_energies()
        self._parse_forces()
        self._parse_max_steps()

    def _parse_xyzs(self):
        species = list()
        geoms = list()

        num_atoms = self.data["num_atoms"]

        line_num = 2
        while not self.lines[line_num].startswith("[GEOCONV]") and line_num < len(self.lines):
            if self.lines[line_num] == str(num_atoms):
                # start new geometry
                geom = list()
                line_num += 2
                for i in range(num_atoms):
                    line = self.lines[line_num + i].split(" ")
                    if len(species) < num_atoms:
                        species.append(line[0])

                    elements = [e for e in line[1:] if e != ""]
                    elements = [float(e) for e in elements]

                    geom.append(elements)

                geoms.append(geom)

                line_num += num_atoms

        self.data["species"] = species
        self.data["geometries"] = list()
        self.data["molecules"] = list()
        for geom in geoms:
            self.data["geometries"].append(np.array(geom))
            self.data["molecules"].append(Molecule(species, geom))

        self.data["num_nodes"] = len(self.data["geometries"])

    def _parse_energies(self):
        for l, line in enumerate(self.lines):
            if line.startswith("energy"):
                line_num = l
                break

        self.data["energies"] = list()
        for i in range(self.data["num_nodes"]):
            self.data["energies"].append(float(self.lines[line_num + i + 1]))

    def _parse_forces(self):
        for l, line in enumerate(self.lines):
            if line.startswith("max-force"):
                line_num = l
                break

        self.data["forces"] = list()
        for i in range(self.data["num_nodes"]):
            self.data["forces"].append(float(self.lines[line_num + i + 1]))

    def _parse_max_steps(self):
        for l, line in enumerate(self.lines):
            if line.startswith("max-step"):
                line_num = l
                break

        self.data["max_steps"] = list()
        for i in range(self.data["num_nodes"]):
            self.data["max_steps"].append(float(self.lines[line_num + i + 1]))

    def as_dict(self):
        d = dict()
        d["data"] = self.data
        d["lines"] = self.lines
        d["text"] = self.text
        d["filename"] = self.filename
        return jsanitize(d, strict=True)


class GSMInternalCoordinateDataParser(MSONable):
    """
    An object that parses and contains data from pyGSM optimized geometry files.

    NOTE: This class is probably even less elegant than the usual Output classes
        (QCOutput, etc.) because the formatting of the optimized string files
        does not lend itself well to table-based or key-based parsing.

    Args:
        filename (str): The path to the pyGSM optimized geometry file.
    """

    def __init__(self, filename):
        self.filename = filename
        self.data = dict()
        self.text = ""
        with zopen(self.filename, "rt") as f:
            self.text = f.read()

        self.data["rct"] = {"translation": list(),
                            "rotation": list(),
                            "distance": list(),
                            "angle": list(),
                            "torsion": list(),
                            "out_of_plane": list()}

        self.data["ts"] = {"translation": list(),
                           "rotation": list(),
                           "distance": list(),
                           "angle": list(),
                           "torsion": list(),
                           "out_of_plane": list()}

        self.data["pro"] = {"translation": list(),
                            "rotation": list(),
                            "distance": list(),
                            "angle": list(),
                            "torsion": list(),
                            "out_of_plane": list()}

        # Parse for node numbers
        temp_node_numbers = read_pattern(self.text, {
            "key": r"Internals\s+minnodeR: ([0-9]+)\s+TSnode: ([0-9]+)\s+minnodeP: ([0-9]+)"
        }, terminate_on_match=True).get("key")
        if temp_node_numbers is None:
            raise ValueError("Could not parse first line; problem with file?")
        else:
            self.data["rct_node"] = int(temp_node_numbers[0][0])
            self.data["ts_node"] = int(temp_node_numbers[0][1])
            self.data["pro_node"] = int(temp_node_numbers[0][2])

        # Parse translations
        temp_translation = read_pattern(self.text, {
            "key": r"Translation-([XYZ]) ([0-9]+)\-([0-9]+)\s+([\.\-0-9]+)\s+([\.\-0-9]+)\s+([\.\-0-9]+)"
        })
        for trans in temp_translation.get("key", []):
            direction = trans[0]
            atom_1 = int(trans[1])
            atom_2 = int(trans[2])

            self.data["rct"]["translation"].append((direction, atom_1, atom_2, float(trans[3])))
            self.data["ts"]["translation"].append((direction, atom_1, atom_2, float(trans[4])))
            self.data["pro"]["translation"].append((direction, atom_1, atom_2, float(trans[5])))

        # Parse rotations
        temp_rotations = read_pattern(self.text, {
            "key": r"Rotation-([ABC]) ([0-9]+)\-([0-9]+)\s+([\.\-0-9]+)\s+([\.\-0-9]+)\s+([\.\-0-9]+)"
        })
        for rot in temp_rotations.get("key", []):
            direction = rot[0]
            atom_1 = int(rot[1])
            atom_2 = int(rot[2])

            self.data["rct"]["rotation"].append((direction, atom_1, atom_2, float(rot[3])))
            self.data["ts"]["rotation"].append((direction, atom_1, atom_2, float(rot[4])))
            self.data["pro"]["rotation"].append((direction, atom_1, atom_2, float(rot[5])))

        # Parse bonds
        temp_distance = read_pattern(self.text, {
            "key": r"Distance ([0-9]+)\-([0-9]+)\s+([\.\-0-9]+)\s+([\.\-0-9]+)\s+([\.\-0-9]+)"
        })
        for dist in temp_distance.get("key", []):
            atom_1 = int(dist[0])
            atom_2 = int(dist[1])

            self.data["rct"]["distance"].append((atom_1, atom_2, float(dist[2])))
            self.data["ts"]["distance"].append((atom_1, atom_2, float(dist[3])))
            self.data["pro"]["distance"].append((atom_1, atom_2, float(dist[4])))

        # Parse angles
        temp_angle = read_pattern(self.text, {
            "key": r"Angle ([0-9]+)\-([0-9]+)\-([0-9]+)\s+([\.\-0-9]+)\s+([\.\-0-9]+)\s+([\.\-0-9]+)"
        })
        for ang in temp_angle.get("key", []):
            atom_1 = int(ang[0])
            atom_2 = int(ang[1])
            atom_3 = int(ang[2])

            self.data["rct"]["angle"].append((atom_1, atom_2, atom_3, float(ang[3])))
            self.data["ts"]["angle"].append((atom_1, atom_2, atom_3, float(ang[4])))
            self.data["pro"]["angle"].append((atom_1, atom_2, atom_3, float(ang[5])))

        # Parse torsions
        temp_dihedral = read_pattern(self.text, {
            "key": r"Dihedral ([0-9]+)\-([0-9]+)\-([0-9]+)\-([0-9]+)\s+([\.\-0-9]+)\s+([\.\-0-9]+)\s+([\.\-0-9]+)"
        })
        for tors in temp_dihedral.get("key", []):
            atom_1 = int(tors[0])
            atom_2 = int(tors[1])
            atom_3 = int(tors[2])
            atom_4 = int(tors[3])

            self.data["rct"]["torsion"].append((atom_1, atom_2, atom_3, atom_4, float(tors[4])))
            self.data["ts"]["torsion"].append((atom_1, atom_2, atom_3, atom_4, float(tors[5])))
            self.data["pro"]["torsion"].append((atom_1, atom_2, atom_3, atom_4, float(tors[6])))

        # Parse out-of-plane angles
        temp_oop = read_pattern(self.text, {
            "key": r"Out\-of\-Plane ([0-9]+)\-([0-9]+)\s+([\.\-0-9]+)\s+([\.\-0-9]+)\s+([\.\-0-9]+)"
        })
        for oop in temp_oop.get("key", []):
            atom_1 = int(oop[0])
            atom_2 = int(oop[1])
            atom_3 = int(oop[2])
            atom_4 = int(oop[3])

            self.data["rct"]["torsion"].append((atom_1, atom_2, atom_3, atom_4, float(oop[4])))
            self.data["ts"]["torsion"].append((atom_1, atom_2, atom_3, atom_4, float(oop[5])))
            self.data["pro"]["torsion"].append((atom_1, atom_2, atom_3, atom_4, float(oop[6])))

    def as_dict(self) -> dict:
        d = dict()
        d["data"] = self.data
        d["text"] = self.text
        d["filename"] = self.filename
        return jsanitize(d, strict=True)