# coding: utf-8

import json
import os

from fireworks import FiretaskBase, FWAction, explicit_serialize
from fireworks.utilities.fw_serializers import DATETIME_HANDLER

from atomate.common.firetasks.glue_tasks import get_calc_loc
from atomate.utils.utils import env_chk
from atomate.utils.utils import get_logger

from mpcat.automate.atomate.database import GSMCalcDb
from mpcat.automate.atomate.drones import GSMDrone

__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2021"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Alpha"
__date__ = "4/24/20"
__credits__ = "Sam Blau"

logger = get_logger(__name__)


@explicit_serialize
class GSMToDb(FiretaskBase):
    """
    Enter a pyGSM run into the database. Uses current directory unless you
    specify calc_dir or calc_loc.

    Optional params:
        calc_dir (str): path to dir (on current filesystem) that contains QChem
            input and output files. Default: use current working directory.
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
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
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        fw_spec_field (str): if set, will update the task doc with the contents
            of this key in the fw_spec.
        additional_fields (dict): dict of additional fields to add
    """
    optional_params = [
        "calc_dir", "calc_loc", "molecule_file", "template_file", "output_file",
        "isomers_file", "db_file", "fw_spec_field", "additional_fields"
    ]

    def run_task(self, fw_spec):
        # get the directory that contains the QChem dir to parse
        calc_dir = os.getcwd()
        if "calc_dir" in self:
            calc_dir = self["calc_dir"]
        elif self.get("calc_loc"):
            calc_dir = get_calc_loc(self["calc_loc"],
                                    fw_spec["calc_locs"])["path"]

        molecule_file = self.get("molecule_file", "input.xyz")
        template_file = self.get("template_file", "qin")
        output_file = self.get("output_file", "gsm.out")
        isomers_file = self.get("isomers_file", "isomers.txt")

        # parse the QChem directory
        logger.info("PARSING DIRECTORY: {}".format(calc_dir))

        additional_fields = self.get("additional_fields", [])

        drone = GSMDrone(additional_fields=additional_fields)

        # assimilate (i.e., parse)
        task_doc = drone.assimilate(
            path=calc_dir, molecule_file=molecule_file,
            template_file=template_file, output_file=output_file,
            isomers_file=isomers_file)

        if "tags" in fw_spec:
            task_doc.update({"tags": fw_spec["tags"]})

        # Check for additional keys to set based on the fw_spec
        if self.get("fw_spec_field"):
            task_doc.update({self.get("fw_spec_field"): fw_spec.get(self.get("fw_spec_field"))})

        # Update fw_spec with final/optimized structure
        update_spec = dict()
        if task_doc.get("output").get("ts_molecule"):
            update_spec["prev_calc_molecule"] = task_doc["output"]["ts_molecule"]
            # update_spec["prev_calc_mulliken"] = task_doc["output"]["mulliken"]
            # if "RESP" in task_doc["output"]:
            #     update_spec["prev_calc_resp"] = task_doc["output"]["RESP"]
            # elif "ESP" in task_doc["output"]:
            #     update_spec["prev_calc_esp"] = task_doc["output"]["ESP"]

        # get the database connection
        db_file = env_chk(self.get("db_file"), fw_spec)

        # db insertion or taskdoc dump
        if not db_file:
            with open(os.path.join(calc_dir, "task.json"), "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
        else:
            mmdb = GSMCalcDb.from_db_file(db_file, admin=True)
            t_id = mmdb.insert(task_doc)
            logger.info("Finished parsing with task_id: {}".format(t_id))

        return FWAction(
            stored_data={"task_id": task_doc.get("task_id", None)},
            update_spec=update_spec)
