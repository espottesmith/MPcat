# coding: utf-8

import os
import datetime
import json
import glob
import traceback

from monty.io import zopen
from monty.json import jsanitize

from pymatgen.apps.borg.hive import AbstractDrone
from pymatgen.io.babel import BabelMolAdaptor
from mpcat.automate.atomate.inputs import (QCTemplate, GSMIsomerInput,
                                    parse_multi_xyz)
from mpcat.automate.atomate.outputs import (GSMOutput,
                                            GSMOptimizedStringParser,
                                            GSMInternalCoordinateDataParser)
from pymatgen.symmetry.analyzer import PointGroupAnalyzer

from atomate.utils.utils import get_logger
from atomate import __version__ as atomate_version

__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2020, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Alpha"
__date__ = "02/18/20"
__credits__ = "Sam Blau, Brandon Wood, Shyam Dwaraknath, Xiaohui Qu, Kiran Mathew, Shyue Ping Ong, Anubhav Jain"

logger = get_logger(__name__)

#TODO: Parse charge properly
# Also as an int (rather than float)


class GSMDrone(AbstractDrone):
    """
    A drone to parse calculations from pyGSM and insert an organized, searchable entry into the database.
    """

    # note: the version is inserted into the task doc
    __version__ = atomate_version

    # Schema def of important keys and sub-keys; used in validation
    schema = {
        "root": {
            "dir_name", "input", "output", "calc", "smiles", "formula_pretty",
            "formula_anonymous", "chemsys", "pointgroup", "formula_alphabetical"
        },
        "input": {"initial_reactants", "initial_products", "mode", "num_nodes",
                  "reactants_fixed", "products_fixed"},
        "output": {"optimized_node_molecules", "ts_node", "ts_molecule",
                   "ts_energy", "absolute_ts_energy"}
    }

    def __init__(self, additional_fields=None):
        """
        Initialize a GSM drone to parse pyGSM calculations
        Args:
            additional_fields (dict): dictionary of additional fields to add to output document
        """

        self.additional_fields = additional_fields or dict()

    def assimilate(self, path, molecule_file="input.xyz", template_file="qin",
                   output_file="gsm.out", isomers_file=None):
        """
        Parses qchem input and output files and insert the result into the db.

        Args:
            path (str): Path to the directory containing output file.
            molecule_file (str): Name of the input molecule geometry file.
                Default is "input.xyz".
            template_file (str): Name of the input QChem template file.
                Default is "qin".
            output_file (str): Name of the pyGSM output file. Default is
                "gsm.out".
            isomers_file (str): For single-ended calculations, this is the
                name of the isomers file that defines what coordinates to
                vary. Default is None; however, note that this should be
                provided for single-ended calculations.

        Returns:
            d (dict): a task dictionary
        """
        logger.info("Getting task doc for base dir :{}".format(path))

        all_files = os.listdir(path)
        all_scratch_files = list()
        if "scratch" in all_files:
            all_scratch_files = os.listdir(os.path.join(path, "scratch"))

        # Populate important input and output files
        mol_file = None
        temp_file = None
        out_file = None
        iso_file = None
        ic_file = None
        opt_file = None
        for file in all_files:
            if isomers_file is not None:
                if file == isomers_file and iso_file is None:
                    iso_file = os.path.join(path, file)
                    continue

            if file == molecule_file and mol_file is None:
                mol_file = os.path.join(path, file)
            elif file == template_file and temp_file is None:
                temp_file = os.path.join(path, file)
            elif file == output_file and out_file is None:
                out_file = os.path.join(path, file)
            elif "IC_data" in file and ic_file is None:
                ic_file = os.path.join(path, file)
            elif "opt_converged" in file and opt_file is None:
                opt_file = os.path.join(path, file)

        # Only check scratch directory if we're missing files
        any_needed_none = any([e is None for e in [mol_file, temp_file,
                                                   out_file, ic_file,
                                                   opt_file]])
        iso_needed_none = isomers_file is not None and iso_file is None

        if any_needed_none or iso_needed_none:
            for file in all_scratch_files:
                if isomers_file is not None:
                    if file == isomers_file and iso_file is None:
                        iso_file = os.path.join(path, "scratch", file)
                        continue

                if file == molecule_file and mol_file is None:
                    mol_file = os.path.join(path, "scratch", file)
                elif file == template_file and temp_file is None:
                    temp_file = os.path.join(path, "scratch", file)
                elif file == output_file and out_file is None:
                    out_file = os.path.join(path, "scratch", file)
                elif "IC_data" in file and ic_file is None:
                    ic_file = os.path.join(path, "scratch", file)
                elif "opt_converged" in file and opt_file is None:
                    opt_file = os.path.join(path, "scratch", file)

        any_needed_none = any([e is None for e in [mol_file, temp_file,
                                                   out_file, opt_file]])
        iso_needed_none = isomers_file is not None and iso_file is None

        if not(any_needed_none or iso_needed_none):
            d = self.generate_doc(path=path, molecule_file=mol_file,
                                  template_file=temp_file, output_file=out_file,
                                  isomers_file=iso_file, internal_coordinate_file=ic_file,
                                  optimized_geom_file=opt_file)
            self.post_process(path, d)
        else:
            raise ValueError("Either input or output not found!")
        self.validate_doc(d)
        return jsanitize(d, strict=True, allow_bson=True)

    def generate_doc(self, path, molecule_file, template_file, output_file,
                     isomers_file, internal_coordinate_file,
                     optimized_geom_file):

        try:
            fullpath = os.path.abspath(path)

            d = jsanitize(self.additional_fields, strict=True)

            d["schema"] = {
                "code": "atomate",
                "version": GSMDrone.__version__
            }

            d["dir_name"] = fullpath

            # TODO: Consider error handlers
            # Include an "orig" section to the doc

            # Parse all relevant files
            initial_mol = parse_multi_xyz(molecule_file)
            for mol in initial_mol:
                mol.set_charge_and_spin(d.get("true_charge", 0))

            temp_file = QCTemplate.from_file(template_file)
            if isomers_file is not None:
                iso_file = GSMIsomerInput.from_file(isomers_file)
            out_file = GSMOutput(output_file)
            if internal_coordinate_file is not None:
                ic_file = GSMInternalCoordinateDataParser(internal_coordinate_file)
            opt_file = GSMOptimizedStringParser(optimized_geom_file)

            d["warnings"] = dict()

            # INPUTS
            d["input"] = dict()
            d["input"]["initial_reactants"] = None
            d["input"]["initial_products"] = None

            if len(initial_mol) == 1:
                d["input"]["initial_reactants"] = initial_mol[0]
            elif len(initial_mol) == 2:
                d["input"]["initial_reactants"] = initial_mol[0]
                d["input"]["initial_products"] = initial_mol[1]

            d["input"]["mode"] = out_file.data["inputs"]["gsm_type"]

            num_nodes = out_file.data["inputs"].get("num_nodes")
            if num_nodes is None:
                if "SE" in d["input"]["mode"]:
                    d["input"]["num_nodes"] = 30
                else:
                    d["input"]["num_nodes"] = 9
            else:
                d["input"]["num_nodes"] = int(num_nodes)

            d["input"]["reactants_fixed"] = out_file.data["inputs"].get("reactant_geom_fixed", False)
            d["input"]["products_fixed"] = out_file.data["inputs"].get("product_geom_fixed", False)

            d["input"]["template"] = {"rem": temp_file.rem,
                                      "pcm": temp_file.pcm,
                                      "solvent": temp_file.solvent,
                                      "smx": temp_file.smx}

            if "SE" in d["input"]["mode"]:
                if isomers_file is None:
                    raise ValueError("No isomers file provided for single-ended calculation.")
                else:
                    d["input"]["isomers"] = {"bonds_formed": iso_file.bonds_formed,
                                             "bonds_broken": iso_file.bonds_broken,
                                             "angles": iso_file.angles,
                                             "torsions": iso_file.torsions,
                                             "out_of_planes": iso_file.out_of_planes}

            d["input"]["parameters"] = out_file.data["inputs"]

            # OUTPUTS
            d["output"] = dict()

            d["output"]["completion"] = out_file.data["completion"]

            if "SE" in d["input"]["mode"]:
                d["output"]["initial_energy"] = out_file.data.get("initial_energy", None)
                d["driving_coord_trajectories"] = out_file.data.get("driving_coord_trajectories", None)
            else:
                d["output"]["initial_energy_rct"] = out_file.data.get("initial_energy_rct", None)
                d["output"]["initial_energy_pro"] = out_file.data.get("initial_energy_pro", None)

            d["output"]["energy_profile"] = out_file.data.get("final_energy_profile", None)
            d["output"]["path_uphill"] = out_file.data.get("final_energy_profile", None)
            d["output"]["path_dissociative"] = out_file.data.get("final_path_dissociative", None)
            d["output"]["minima_nodes"] = out_file.data.get("final_min_nodes", None)
            d["output"]["maxima_nodes"] = out_file.data.get("final_max_nodes", None)
            d["output"]["minima_nodes"] = out_file.data.get("final_min_nodes", None)
            d["output"]["maximum_node"] = out_file.data.get("final_max_node", None)
            d["output"]["maximum_energy"] = out_file.data.get("final_max_energy", None)

            if d["output"]["completion"]:
                d["output"]["reactant_node"] = out_file.data["min_rct_node"]
                d["output"]["product_node"] = out_file.data["min_pro_node"]
                d["output"]["ts_node"] = out_file.data["ts_node"]
                d["output"]["absolute_ts_energy"] = out_file.data["absolute_ts_energy"]
                d["output"]["ts_energy"] = out_file.data["ts_energy"]
                d["output"]["delta_e"] = out_file.data["delta_e"]
            else:
                d["output"]["reactant_node"] = None
                d["output"]["product_node"] = None
                d["output"]["ts_node"] = None
                d["output"]["ts_energy"] = None
                d["output"]["absolute_ts_energy"] = None
                d["output"]["delta_e"] = None

            if d["output"]["completion"]:
                if internal_coordinate_file is not None:
                    d["output"]["internal_coords"] = ic_file.data
                else:
                    d["output"]["internal_coords"] = None
                d["output"]["species"] = opt_file.data["species"]
                d["output"]["optimized_node_geometries"] = opt_file.data["geometries"]
                d["output"]["optimized_node_molecules"] = opt_file.data["molecules"]
                for mol in d["output"]["optimized_node_molecules"]:
                    mol.set_charge_and_spin(d.get("true_charge", 0))
                d["output"]["optimized_node_energies"] = opt_file.data["energies"]
                d["output"]["optimized_node_forces"] = opt_file.data["forces"]
                if d["output"]["ts_node"] is not None:
                    d["output"]["ts_molecule"] = d["output"]["optimized_node_molecules"][d["output"]["ts_node"]]
                    d["output"]["ts_molecule"].set_charge_and_spin(d.get("true_charge", 0))
                else:
                    d["output"]["ts_molecule"] = None
                if d["output"]["reactant_node"] is not None:
                    d["output"]["reactant_molecule"] = d["output"]["optimized_node_molecules"][d["output"]["reactant_node"]]
                    d["output"]["reactant_molecule"].set_charge_and_spin(d.get("true_charge", 0))
                else:
                    d["output"]["reactant_molecule"] = None
                if d["output"]["product_node"] is not None:
                    d["output"]["product_molecule"] = d["output"]["optimized_node_molecules"][d["output"]["product_node"]]
                    d["output"]["product_molecule"].set_charge_and_spin(d.get("true_charge", 0))
                else:
                    d["output"]["product_molecule"] = None
            else:
                d["output"]["internal_coords"] = None
                d["output"]["species"] = None
                d["output"]["optimized_node_geometries"] = None
                d["output"]["optimized_node_molecules"] = None
                d["output"]["optimized_node_energies"] = None
                d["output"]["optimized_node_forces"] = None
                d["output"]["ts_molecule"] = None
                d["output"]["reactant_molecule"] = None
                d["output"]["product_molecule"] = None

            d["calc"] = out_file.data

            d["warnings"] = out_file.data["warnings"]
            d["errors"] = out_file.data["errors"]

            # if d_calc_final["completion"]:
            #     total_cputime = 0.0
            #     total_walltime = 0.0
            #     for calc in d["calcs_reversed"]:
            #         if "walltime" in calc and "cputime" in calc:
            #             if calc["walltime"] is not None:
            #                 total_walltime += calc["walltime"]
            #             if calc["cputime"] is not None:
            #                 total_cputime += calc["cputime"]
            #     d["walltime"] = total_walltime
            #     d["cputime"] = total_cputime
            # else:
            #     d["walltime"] = None
            #     d["cputime"] = None

            comp = d["input"]["initial_reactants"].composition
            d["formula_pretty"] = comp.reduced_formula
            d["formula_anonymous"] = comp.anonymized_formula
            d["formula_alphabetical"] = comp.alphabetical_formula

            elements = list()
            for component in d["formula_alphabetical"].split(" "):
                elements.append("".join([i for i in component if not i.isdigit()]))
            d["chemsys"] = "-".join(sorted(set(elements)))

            if d["output"]["ts_molecule"] is not None:
                try:
                    d["pointgroup_ts"] = PointGroupAnalyzer(d["output"]["ts_molecule"]).sch_symbol
                except ValueError:
                    d["pointgroup_ts"] = "PGA_error"
            else:
                d["pointgroup_ts"] = None

            if d["output"]["reactant_molecule"] is not None:
                try:
                    d["pointgroup_reactant"] = PointGroupAnalyzer(d["output"]["reactant_molecule"]).sch_symbol
                except ValueError:
                    d["pointgroup_reactant"] = "PGA_error"
            else:
                d["pointgroup_reactant"] = None

            if d["output"]["product_molecule"] is not None:
                try:
                    d["pointgroup_product"] = PointGroupAnalyzer(d["output"]["product_molecule"]).sch_symbol
                except ValueError:
                    d["pointgroup_product"] = "PGA_error"
            else:
                d["pointgroup_product"] = None

            if d["output"]["ts_molecule"] is not None:
                bb = BabelMolAdaptor(d["output"]["ts_molecule"])
                pbmol = bb.pybel_mol
                smiles = pbmol.write(str("smi")).split()[0]
                d["smiles"] = smiles
            else:
                d["smiles"] = None

            d["state"] = "successful" if d["output"]["completion"] else "unsuccessful"

            d["last_updated"] = datetime.datetime.utcnow()
            return d

        except Exception:
            logger.error(traceback.format_exc())
            logger.error("Error in " + os.path.abspath(path) + ".\n" +
                         traceback.format_exc())
            raise

    @staticmethod
    def post_process(path, d):
        """
        Post-processing for various files other than the QChem input and output files.
        """
        logger.info("Post-processing dir:{}".format(path))
        fullpath = os.path.abspath(path)
        filenames = glob.glob(os.path.join(fullpath, "custodian.json*"))
        if len(filenames) >= 1:
            with zopen(filenames[0], "rt") as f:
                d["custodian"] = json.load(f)
        filenames = glob.glob(os.path.join(fullpath, "solvent_data*"))
        if len(filenames) >= 1:
            with zopen(filenames[0], "rt") as f:
                d["custom_smd"] = f.readlines()[0].strip()

    def validate_doc(self, d):
        """
        Sanity check, aka make sure all the important keys are set. Note that a failure
        to pass validation is unfortunately unlikely to be noticed by a user.
        """
        for k, v in self.schema.items():
            diff = v.difference(set(d.get(k, d).keys()))
            if diff:
                logger.warning("The keys {0} in {1} not set".format(diff, k))

    @staticmethod
    def get_valid_paths(self, path):
        return [path]
