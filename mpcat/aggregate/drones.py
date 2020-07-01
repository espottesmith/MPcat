# coding: utf-8

import os
import datetime

from typing import Optional

from monty.json import jsanitize

from pymatgen.apps.borg.hive import AbstractDrone

from mpcat.adapt.schrodinger_adapter import maestro_file_to_molecule
from mpcat.apprehend.autots_input import AutoTSInput
from mpcat.apprehend.autots_output import AutoTSOutput
from mpcat.apprehend.jaguar_output import JagOutput


version = "0.0.1"


class AutoTSDrone(AbstractDrone):
    """
    A drone to parse AutoTS calculations and convert them into a task document.
    """

    schema = {
        "root": {
            "path", "input", "output", "calcs", "walltime"
        },
        "input": {"reactants", "products", "autots_variables", "gen_variables"},
    }

    def __init__(self, path: str):
        """
        Args:
            path (str): Path to root directory where the AutoTS was conducted.
        """
        self.path = path

        contents = os.listdir(self.path)

        self.directories = [c for c in contents if os.path.isdir(os.path.join(path, c))]
        self.files = [c for c in contents if os.path.isfile(os.path.join(path, c))]

        self.input_files = [f for f in self.files if f.endswith(".in") or f.endswith(".in.gz")]
        self.output_files = [f for f in self.files if f.endswith(".out") or f.endswith(".out.gz")]
        self.mae_files = [f for f in self.files if f.endswith(".mae") or f.endswith(".mae.gz")]

    def assimilate(self, workflow_input: Optional[str] = "autots.in",
                     workflow_output: Optional[str] = None):
        """
        Generate the task doc and perform any additional steps on it to prepare
        for insertion into a DB.

        Args:
            workflow_input (str): AutoTS input file. Default is "autots.in"
            workflow_output (str): AutoTS output file. Default is None, meaning
                that the AutoTSDrone will infer the output file.

        Returns:
            task_doc (dict): The compiled results from the calculation.
        """

        d = self.generate_doc(workflow_input, workflow_output)
        self.validate_doc(d)
        return jsanitize(d, strict=True, allow_bson=True)

    def generate_doc(self, workflow_input: Optional[str] = "autots.in",
                     workflow_output: Optional[str] = None):
        """
        Generate a dictionary from the inputs and outputs of the various
        Jaguar calculations involved in the AutoTS calculation.

        Args:
            workflow_input (str): AutoTS input file. Default is "autots.in"
            workflow_output (str): AutoTS output file. Default is None, meaning
                that the AutoTSDrone will infer the output file.

        Returns:
            d (dict): The compiled results from the calculation.
        """

        if workflow_input not in self.files:
            raise ValueError("Invalid input file given for workflow_input; "
                             "input file is not in path!")

        d = dict()
        d["schema"] = {"code": "mpcat", "version": version}
        d["path"] = self.path

        autots_input = AutoTSInput.from_file(workflow_input, read_molecules=True)

        d["input"] = {"reactants": autots_input.reactants,
                      "products": autots_input.products,
                      "autots_variables": autots_input.autots_variables,
                      "gen_variables": autots_input.gen_variables}

        if workflow_output is not None:
            if workflow_output in self.files:
                chosen_output = workflow_output
            else:
                raise ValueError("Invalid output file given for workflow_output; "
                                 "output file is not in path!")
        else:
            if len(self.output_files) == 0:
                raise ValueError("No output files found!")
            if len(self.output_files) == 1:
                chosen_output = self.output_files[0]
            else:
                probable_outputs = [o for o in self.output_files if "AutoTS" in o]
                if len(probable_outputs) > 1:
                    raise ValueError("Unclear what the correct AutoTS output "
                                     "file is! Please provide a value for "
                                     "workflow_output")
                else:
                    chosen_output = probable_outputs[0]

        autots_output = AutoTSOutput(chosen_output)

        d["calculation_names"] = autots_output.data["calculations"]

        d["warnings"] = autots_output.data["warnings"]
        d["errors"] = autots_output.data["errors"]
        d["completed"] = autots_output.data["complete"]
        d["walltime"] = autots_output.data["walltime"]
        if d["completed"]:
            d["output"] = {"temperature": autots_output.data["output"]["temperature_dependence"],
                           "free_energy": autots_output.data["output"]["gibbs_energies"],
                           "separated_energies": autots_output.data["output"]["separated_energies"],
                           "complex_energies": autots_output.data["output"]["complex_energies"],
                           "optimized_structure_energies": autots_output.data["output"]["struct_energies"]}
            d["diagnostic"] = autots_output.data["diagnostic"]
        else:
            d["output"] = None
            #TODO: could actually parse diagnostic for failed calcs; just requires tweaking of regex
            d["diagnostic"] = None

        d["calcs"] = list()
        for calculation in d["calculation_names"]:
            found = False
            for directory in self.directories:
                directory_files = os.listdir(os.path.join(self.path,
                                                          directory))
                if (calculation + ".out") in directory_files:
                    found = True
                    jag_out = JagOutput(os.path.join(self.path, directory,
                                                     calculation + ".out"),
                                        allow_failure=True,
                                        parse_molecules=True)
                    d["calcs"].append(jag_out.data)

            if not found:
                raise RuntimeWarning("Expected calculation {} missing from path!".format(calculation))

        full_paths = [f for f in self.files if f.endswith("full_path.mae")]
        if len(full_paths) > 0:
            d["output"]["path_molecules"] = maestro_file_to_molecule(os.path.join(self.path,
                                                                                  full_paths[0]))

        rct_complexes = [f for f in self.files if f.endswith("reactant_complex.mae")]
        if len(rct_complexes) > 0:
            d["output"]["reactant_complex"] = maestro_file_to_molecule(os.path.join(self.path,
                                                                                    rct_complexes[0]))[0]

        pro_complexes = [f for f in self.files if f.endswith("product_complex.mae")]
        if len(rct_complexes) > 0:
            d["output"]["product_complex"] = maestro_file_to_molecule(os.path.join(self.path,
                                                                                   pro_complexes[0]))[0]

        d["last_updated"] = datetime.datetime.utcnow()
        return d

    def validate_doc(self, d):
        """
        Sanity check, aka make sure all the important keys are set. Note that a failure
        to pass validation is unfortunately unlikely to be noticed by a user.
        """
        for k, v in self.schema.items():
            diff = v.difference(set(d.get(k, d).keys()))
            if diff:
                raise RuntimeWarning("The keys {0} in {1} not set".format(diff, k))