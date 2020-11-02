# coding: utf-8

import os
from datetime import datetime

from typing import Optional, List, Dict
from pathlib import Path
import hashlib

from monty.json import jsanitize
from monty.serialization import loadfn

from pymatgen.apps.borg.hive import AbstractDrone
# from maggma.core.drone import Document, Drone, RecordIdentifier

from mpcat.adapt.schrodinger_adapter import maestro_file_to_molecule
from mpcat.apprehend.autots_input import AutoTSInput
from mpcat.apprehend.autots_output import AutoTSOutput
from mpcat.apprehend.jaguar_output import JagOutput
from mpcat.aggregate.database import CatDB


version = "0.0.1"


# Credit to the Materials Project and the Maggma development team for compute_state_hash
def compute_state_hash(documents: List[Path]) -> str:
    """
    compute the hash of the state of some set of documents

    This functions is taken from maggma.core.drone.RecordIdentifier. All credit
    belongs to the Materials Project.

    Args:
        documents (list of Paths): Documents to be digested and conveted into
            a hash

    Returns:
        hash of the list of documents passed in
    """
    digest = hashlib.md5()
    for doc in documents:
        digest.update(doc.name.encode())
        with open(doc.as_posix(), "rb") as file:
            buf = file.read()
            digest.update(buf)
    return str(digest.hexdigest())


class AutoTSCalcDrone(AbstractDrone):
    """
    A drone to parse AutoTS calculations and convert them into a task document.
    """

    schema = {
        "root": {
            "path", "input", "output", "calcs", "walltime", "hash"
        },
        "input": {"reactants", "products", "autots_variables", "gen_variables"},
    }

    def __init__(self, path: Path):
        """
        Args:
            path (str): Path to root directory where the AutoTS was conducted.
        """
        self.path = path

        self.documents = self.get_documents_calc_dir(path)
        self.file_names = [d.name for d in self.documents]

    @staticmethod
    def get_documents_calc_dir(calc_dir: Path) -> List[Path]:
        """
        Identify all important documents within an AutoTS calculation directory.

        Args:
            calc_dir: Path object that indicates a path to an individual AutoTS
                calculation

        Returns:
            List of Documents
        """

        root_contents = [f for f in calc_dir.iterdir()]

        subdirs = [f for f in root_contents if (calc_dir / f).is_dir() and "AutoTS" in f.name]
        files = [f for f in root_contents if (calc_dir / f).is_file()]

        allowed_suffixes = ["mae", "out", "in", "mae.gz", "out.gz", "in.gz", ".json"]

        files_paths = list()
        for file in files:
            f = file.as_posix()
            if "AutoTS" in f and any([f.endswith(x) for x in allowed_suffixes]):
                files_paths.append(calc_dir / file)
            elif ("pro" in f or "rct" in f) and f.split(".", maxsplit=2)[-1] in ["mae",
                                                                                 "mae.gz"]:
                files_paths.append(calc_dir / file)
            elif f.endswith(".in") or f.endswith(".in.gz"):
                files_paths.append(calc_dir / file)

        for subdir in subdirs:
            sub_files = [f.as_posix() for f in (calc_dir / subdir).iterdir() if (calc_dir / subdir / f).is_file()]
            sub_paths = [subdir / f for f in sub_files]
            for ff, file in enumerate(sub_files):
                if any([f.endswith(x) for x in allowed_suffixes]):
                    files_paths.append(sub_paths[ff])

        return files_paths

    def assimilate(self):
        """
        Generate the task doc and perform any additional steps on it to prepare
        for insertion into a DB.

        Args:
            None

        Returns:
            task_doc (dict): The compiled results from the calculation.
        """

        d = self.generate_doc()
        self.validate_doc(d)
        return jsanitize(d, strict=True, allow_bson=True)

    def generate_doc(self):
        """
        Generate a dictionary from the inputs and outputs of the various
        Jaguar calculations involved in the AutoTS calculation.

        Args:
            None

        Returns:
            d (dict): The compiled results from the calculation.
        """

        if "autots.in" not in self.file_names:
            raise ValueError("Input file is not in path!")

        d = dict()
        d["schema"] = {"code": "mpcat", "version": version}
        d["path"] = self.path.as_posix()
        d["hash"] = compute_state_hash(self.documents)

        calc_data = dict()
        for document in self.documents:
            if document.name == "calc.json":
                calc_data = loadfn(document.as_posix())
                break

        d["rxnid"] = calc_data.get("rxnid")
        d["tags"] = calc_data.get("tags")
        d["name"] = calc_data.get("name")
        d["charge"] = calc_data.get("charge")
        d["spin_multiplicity"] = calc_data.get("spin_multiplicity")
        d["nelectrons"] = calc_data.get("nelectrons")

        autots_input = AutoTSInput.from_file((self.path / "autots.in").as_posix(),
                                             read_molecules=True)

        d["input"] = {"reactants": autots_input.reactants,
                      "products": autots_input.products,
                      "autots_variables": autots_input.autots_variables,
                      "gen_variables": autots_input.gen_variables}

        output_document = None
        for document in self.documents:
            if document.name == "AutoTS.T9XnCsLi.out":
                output_document = document
                break
        if output_document is None:
            raise ValueError("Output file is not in path!")

        autots_output = AutoTSOutput(output_document.as_posix())

        d["calculation_names"] = autots_output.data.get("calculations")
        if d["calculation_names"] is None:
            d["calculation_names"] = list()

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
            d["output"] = dict()
            #TODO: could actually parse diagnostic for failed calcs; just requires tweaking of regex
            d["diagnostic"] = dict()

        d["calcs"] = list()
        for calculation in d.get("calculation_names", list()):
            found = False
            for document in self.documents:
                if (calculation + ".out") in document.name:
                    found = True
                    jag_out = JagOutput(document.as_posix())
                    d["calcs"].append(jag_out.data)

            if not found:
                print("Expected calculation {} missing from path {}!".format(calculation, self.path.as_posix()))

        full_paths = [d for d in self.documents if "full_path.mae" in d.name]
        if len(full_paths) > 0:
            d["output"]["path_molecules"] = maestro_file_to_molecule(full_paths[0].as_posix())

        rct_complexes = [d for d in self.documents if "reactant_complex.mae" in d.name]
        if len(rct_complexes) > 0:
            d["output"]["reactant_complex"] = maestro_file_to_molecule(rct_complexes[0].as_posix())[0]

        pro_complexes = [d for d in self.documents if "product_complex.mae" in d.name]
        if len(pro_complexes) > 0:
            d["output"]["product_complex"] = maestro_file_to_molecule(pro_complexes[0].as_posix())[0]

        d["last_updated"] = datetime.now()
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

    @staticmethod
    def get_valid_paths(path):
        return [path]


class AutoTSBuilderDrone:
    """
    A drone to parse through a collection of many AutoTS calculations and enter
    them into a database.
    """

    def __init__(self, db: CatDB, path: Path):
        """

        Args:
            db (CatDB): Database connection for storing calculations.
            path (Path): Path to the root directory where calculations are
                stored.
        """

        self.db = db
        self.path = path

    def find_valid_directories(self):
        """
        Examine all subdirectories in the main path, and determine which of them
        are expected to be valid calculation directories.

        Args:
             None

        Returns:
            valid (list of Paths): List of directories to be parsed
        """

        directories = [e for e in self.path.iterdir() if e.is_dir()]

        valid = list()

        for directory in directories:
            file_names = [e.name for e in directory.iterdir() if e.is_file()]
            in_found = False
            out_found = False

            for file_name in file_names:
                if "autots.in" in file_name:
                    in_found = True
                if "AutoTS.T9XnCsLi.out" in file_name:
                    out_found = True
                if in_found and out_found:
                    break

            if in_found and out_found:
                valid.append(directory)

        return valid

    def read(self):
        """
        Determine the state hashes for all of the directories in the path

        Args:
            None

        Returns:
             mapping (dict): A dictionary where keys are Paths to specific
                calculation directories and values are task docs, generated
                using AutoTSCalcDrone
        """

        valid = self.find_valid_directories()

        mapping = dict()
        for directory in valid:
            documents = AutoTSCalcDrone.get_documents_calc_dir(directory)
            mapping[directory] = compute_state_hash(documents)

        return mapping

    def find_records_to_update(self, records: Dict):
        """
        Determine which records need to be updated, based on a hashing function.

        Args:
            records (Dict): A dictionary where keys are Paths and values are
                hashes

        Return:
            to_update (Dict): A dictionary of the same format as records.
        """

        paths = [e for e in records.keys()]
        paths_names = [e.as_posix() for e in paths]

        cursor = self.db.database[self.db.data_collection].find(
            {"path": {"$in": paths_names}},
            {"_id": 0, "path": 1, "hash": 1, "last_updated": 1}
        )

        db_record_log = {doc["path"]: doc["hash"] for doc in cursor}

        to_update_list = [hash_value != db_record_log.get(path, None)
                          for path, hash_value in records.items()]

        return [path for path, to_update in zip(paths, to_update_list) if to_update]

    def update_targets(self, items: List[Path]):
        """
        Use AutoTSCalcDrone to update select entries in the database

        Args:
            items (list of Paths): Paths to be parsed/updated

        Return:
            None
        """

        docs = list()

        for path in items:
            drone = AutoTSCalcDrone(path)
            try:
                doc = drone.assimilate()
                docs.append(doc)
            except ValueError:
                print("Cannot parse {}".format(path.as_posix()))
                continue

        if len(docs) > 0:
            self.db.update_data_doc(docs, key="path")

    def build(self):
        """
        Perform the build sequence - find the relevant directories, determine
            if they need to be updated, and then perform the update.

        Args:
            None

        Return:
            None
        """

        mapping = self.read()
        to_update = self.find_records_to_update(mapping)

        self.update_targets(to_update)
