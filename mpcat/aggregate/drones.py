# coding: utf-8

from datetime import datetime

from typing import List, Dict, Optional, Union
from pathlib import Path
import hashlib
from json import loads

from pymongo import ReturnDocument

from monty.json import jsanitize
from monty.serialization import loadfn

from pymatgen.apps.borg.hive import AbstractDrone

from schrodinger.application.jaguar.textparser import JaguarParseError

from mpcat.adapt.schrodinger_adapter import maestro_file_to_molecule
from mpcat.apprehend.jaguar_input import JagInput
from mpcat.apprehend.jaguar_output import JagOutput, JaguarOutputParseError
from mpcat.apprehend.autots_input import TSInput
from mpcat.apprehend.autots_output import TSOutput
from mpcat.aggregate.database import CatDB
from mpcat.utils.types import JaguarJobType, job_type_mapping


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


class JaguarCalcDrone(AbstractDrone):

    schema = {
        "root": {
            "path", "input", "output", "walltime", "hash"
        },
        "input": {"molecule", "gen_variables"},
    }

    schema_traj = {
        "root": {
            "path",
            "input",
            "energy_trajectory",
            "gradient_trajectory",
            "molecule_trajectory",
            "origin",
        },
        "input": {"molecule", "gen_variables"},
    }

    def __init__(self,
                 path: Path,
                 prefix: str = "jaguar",
                 job_type: Optional[Union[str, JaguarJobType]] = None):
        self.path = path
        self.prefix = prefix

        if isinstance(job_type, JaguarJobType):
            self.job_type = job_type
        elif job_type is None:
            self.job_type = None
        else:
            self.job_type = job_type_mapping[job_type]

        self.documents = self.get_documents_calc_dir(path)
        self.file_names = [d.name for d in self.documents]

    @staticmethod
    def get_documents_calc_dir(calc_dir: Path) -> List[Path]:
        """
        Identify all important documents within a Jaguar calculation directory.

        Args:
            calc_dir: Path object that indicates a path to an individual AutoTS
                calculation

        Returns:
            List of Documents
        """
        files = [e for e in calc_dir.iterdir() if e.is_file()]

        allowed_suffixes = ["mae", "out", "in", "xyz", "mae.gz", "out.gz", "in.gz"]

        files_paths = list()
        for file in files:
            f = file.as_posix()
            if prefix in f and any([f.endswith(x) for x in allowed_suffixes]):
                files_paths.append(file)
            elif "calc.json" in f:
                files_paths.append(file)

        return files_paths

    def assimilate(self, parse_molecules=True):
        """
        Generate the task doc and perform any additional steps on it to prepare
        for insertion into a DB.

        Args:
            None

        Returns:
            task_doc (dict): The compiled results from the calculation.
        """

        d = self.generate_doc(parse_molecules=parse_molecules)
        self.validate_doc(d, self.schema)
        return jsanitize(d, strict=True, allow_bson=True)

    def assimilate_trajectory(self):
        """
        Generate the trajectory doc and perform any additional steps on it to prepare
        for insertion into a DB.

        Args:
            None

        Returns:
            traj_doc (dict): The compiled trajectory data from a geometry optimization-type
                calculation.
        """

        d = self.generate_trajectory_doc()
        self.validate_doc(d, self.schema_traj)
        return jsanitize(d, strict=True, allow_bson=True)

    def generate_doc(self, parse_molecules=True):
        """
        Generate a dictionary from the inputs and outputs of the various
        Jaguar calculations involved in the Jaguar calculation.

        Args:
            None

        Returns:
            d (dict): The compiled results from the calculation.
        """

        if "jaguar.in" not in self.file_names:
            print(self.file_names)
            raise ValueError("Input file is not in path!")

        d = dict()
        d["schema"] = {"code": "mpcat", "version": version}
        d["path"] = self.path.as_posix()
        d["hash"] = compute_state_hash(self.documents)

        calc_data = dict()
        for document in self.documents:
            if document.name == "calc.json":
                calc_data = loads(open(document.as_posix()).read())
                break

        d["calcid"] = calc_data.get("calcid")
        d["tags"] = calc_data.get("tags")
        d["additional_data"] = calc_data.get("additional_data")
        d["name"] = calc_data.get("name")
        d["charge"] = calc_data.get("charge")
        d["spin_multiplicity"] = calc_data.get("spin_multiplicity")
        d["nelectrons"] = calc_data.get("nelectrons")
        d["job_type"] = calc_data.get("job_type")
        d["formula_alphabetical"] = calc_data.get("formula_alphabetical")

        jaguar_input = JagInput.from_file(self.path / "jaguar.in")

        output_document = None
        for document in self.documents:
            if document.name == "jaguar.out":
                output_document = document
                break
        if output_document is None:
            raise ValueError("Output file is not in path!")

        jaguar_output = JagOutput(output_document.as_posix(),
                                  parse_molecules=parse_molecules)

        d["job_name"] = jaguar_output.data.get("job_name")
        d["full_filename"] = jaguar_output.data.get("full_filename")
        d["job_id"] = jaguar_output.data.get("job_id")

        d["errors"] = dict()
        d["errors"]["parsing"] = jaguar_output.data.get("parsing_error")
        d["errors"]["fatal"] = jaguar_output.data.get("fatal_error")

        d["success"] = jaguar_output.data["success"]
        d["walltime"] = jaguar_output.data["walltime"]

        d["input"] = jaguar_output.data["input"]
        d["input"]["molecule"] = jaguar_input.mol
        d["input"]["gen_variables"] = jaguar_input.gen_variables

        if d.get("formula_alphabetical") is None:
            d["formula_alphabetical"] = d["input"]["molecule"].composition.alphabetical_formula

        d["output"] = jaguar_output.data["output"]

        d["last_updated"] = datetime.now()
        return d

    def generate_trajectory_doc(self, origin=None):
        doc = self.generate_doc(parse_molecules=True)

        doc["energy_trajectory"] = doc["output"]["energy_trajectory"]
        doc["gradient_trajectory"] = doc["output"]["gradient_trajectory"]
        doc["molecule_trajectory"] = doc["output"]["molecule_trajectory"]
        del doc["output"]

        doc["origin"] = dict()
        if doc.get("calcid") is not None:
            doc["origin"]["calcid"] = doc["calcid"]
            del doc["calcid"]
        elif origin is not None:
            doc["origin"] = origin

        return doc

    def validate_doc(self, d: Dict, schema: Dict):
        """
        Sanity check, aka make sure all the important keys are set. Note that a failure
        to pass validation is unfortunately unlikely to be noticed by a user.
        """

        for k, v in schema.items():
            diff = v.difference(set(d.get(k, d).keys()))
            if diff:
                raise RuntimeWarning("The keys {0} in {1} not set".format(diff, k))

    @staticmethod
    def get_valid_paths(path):
        return [path]


class JaguarBuilderDrone:
    """
    A drone to parse through a collection of many Jaguar calculations and enter
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
            this_valid = False

            if "jaguar.in" in file_names and "jaguar.out" in file_names:
                this_valid = True
            elif "calc.json" in file_names:
                calc_data = loadfn((directory / "calc.json").as_posix())
                in_file = "{}.in".format(calc_data["calcid"])
                out_file = "{}.out".format(calc_data["calcid"])
                if in_file in file_names and out_file in file_names:
                    this_valid = True

            if this_valid:
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
                using JaguarCalcDrone
        """

        valid = self.find_valid_directories()

        mapping = dict()
        for directory in valid:
            documents = JaguarCalcDrone.get_documents_calc_dir(directory)
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

        cursor = self.db.database[self.db.jaguar_data_collection].find(
            {"path": {"$in": paths_names}},
            {"_id": 0, "path": 1, "hash": 1, "last_updated": 1}
        )

        db_record_log = {doc["path"]: doc["hash"] for doc in cursor}

        to_update_list = [hash_value != db_record_log.get(path, None)
                          for path, hash_value in records.items()]

        return [path for path, to_update in zip(paths, to_update_list) if to_update]

    def update_targets(self, items: List[Path]):
        """
        Use JaguarCalcDrone to update select entries in the database

        Args:
            items (list of Paths): Paths to be parsed/updated

        Return:
            None
        """

        docs = list()

        for path in items:
            drone = JaguarCalcDrone(path)
            try:
                doc = drone.assimilate(parse_molecules=False)
                docs.append(doc)
            except:
                print("Cannot parse {}".format(path.as_posix()))
                continue

        if len(docs) > 0:
            self.db.update_jaguar_data_docs(docs, key="path")

    def build(self, parse_molecules=True):
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


class JaguarTrajectoryDrone:
    """
    A drone to parse through a collection of many Jaguar calculations and enter
    any trajectory information (from a PES scan, geometry optimization,
    TS optimization, or IRC calculation) into a database.
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
            this_valid = False

            if "jaguar.in" in file_names and "jaguar.out" in file_names:
                this_valid = True
            elif "calc.json" in file_names:
                calc_data = loadfn((directory / "calc.json").as_posix())
                in_file = "{}.in".format(calc_data["calcid"])
                out_file = "{}.out".format(calc_data["calcid"])

                job_type = calc_data.get("job_type", "sp")
                if in_file in file_names and out_file in file_names:
                    if job_type in ["opt", "ts", "scan", "irc"]:
                        this_valid = True

            if this_valid:
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
                using JaguarCalcDrone
        """

        valid = self.find_valid_directories()

        mapping = dict()
        for directory in valid:
            documents = JaguarCalcDrone.get_documents_calc_dir(directory)
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

        cursor = self.db.database[self.db.trajectory_collection].find(
            {"path": {"$in": paths_names}},
            {"_id": 0, "path": 1, "hash": 1, "last_updated": 1}
        )

        db_record_log = {doc["path"]: doc["hash"] for doc in cursor}

        to_update_list = [hash_value != db_record_log.get(path, None)
                          for path, hash_value in records.items()]

        return [path for path, to_update in zip(paths, to_update_list) if to_update]

    def update_targets(self, items: List[Path], parse_molecules=True):
        """
        Use JaguarCalcDrone to update select entries in the database

        Args:
            items (list of Paths): Paths to be parsed/updated

        Return:
            None
        """

        docs = list()

        for path in items:
            drone = JaguarCalcDrone(path)
            try:
                doc = drone.assimilate_trajectory()
                docs.append(doc)
            except Exception as e:
                print("Cannot parse {}: {}".format(path.as_posix(), e))
                continue

        if len(docs) > 0:
            self.db.update_trajectory_docs(docs, key="path")

    def build(self, parse_molecules=True):
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

        self.update_targets(to_update, parse_molecules=parse_molecules)


class AutoTSCalcDrone(AbstractDrone):
    """
    A drone to parse AutoTS calculations and convert them into task documents.
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

        allowed_suffixes = ["mae", "out", "in", "mae.gz", "out.gz", "in.gz"]

        files_paths = list()
        for file in files:
            f = file.as_posix()
            if "AutoTS" in f and any([f.endswith(x) for x in allowed_suffixes]):
                files_paths.append(calc_dir / file)
            elif ("pro" in f or "rct" in f) and (f.endswith("mae") or f.endswith("mae.gz")):
                files_paths.append(calc_dir / file)
            elif f.endswith(".in") or f.endswith(".in.gz"):
                files_paths.append(calc_dir / file)
            elif "calc.json" in f:
                files_paths.append(calc_dir / file)

        for subdir in subdirs:
            sub_files = [f.as_posix() for f in (calc_dir / subdir).iterdir() if (calc_dir / subdir / f).is_file()]
            sub_paths = [subdir / f for f in sub_files]
            for ff, file in enumerate(sub_files):
                if any([file.endswith(x) for x in allowed_suffixes]):
                    files_paths.append(sub_paths[ff])

            sub_subdirs = [f.as_posix() for f in (calc_dir / subdir).iterdir() if (calc_dir / subdir / f).is_dir()]
            for sub_subdir in sub_subdirs:
                subsub_files = [f.as_posix() for f in (calc_dir / subdir / sub_subdir).iterdir() if (calc_dir / subdir / sub_subdir / f).is_file()]
                subsub_paths = [subdir/ sub_subdir / f for f in subsub_files]
                for ff, file in enumerate(subsub_files):
                    if any([file.endswith(x) for x in allowed_suffixes]):
                        files_paths.append(subsub_paths[ff])

        return files_paths

    def assimilate(self, parse_molecules=True):
        """
        Generate the task doc and perform any additional steps on it to prepare
        for insertion into a DB.

        Args:
            None

        Returns:
            task_doc (dict): The compiled results from the calculation.
        """

        d = self.generate_doc(parse_molecules=parse_molecules)
        self.validate_doc(d)
        return jsanitize(d, strict=True, allow_bson=True)

    def generate_doc(self, parse_molecules=True):
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
                calc_data = loads(open(document.as_posix()).read())
                break

        d["rxnid"] = calc_data.get("rxnid")
        d["tags"] = calc_data.get("tags")
        d["additional_data"] = calc_data.get("additional_data")
        d["name"] = calc_data.get("name")
        d["charge"] = calc_data.get("charge")
        d["spin_multiplicity"] = calc_data.get("spin_multiplicity")
        d["nelectrons"] = calc_data.get("nelectrons")

        autots_input = TSInput.from_file(self.path / "autots.in",
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

        autots_output = TSOutput(output_document.as_posix())

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
                    try:
                        jag_out = JagOutput(document.as_posix(),
                                            parse_molecules=parse_molecules)
                        d["calcs"].append(jag_out.data)
                    except (JaguarParseError, JaguarOutputParseError):
                        print("Error parsing " + calculation + " in path " + self.path.as_posix())
                        break

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

    def validate_doc(self, d: Dict):
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

    def __init__(self, db: CatDB, path: Path, parse_molecules: bool = True):
        """

        Args:
            db (CatDB): Database connection for storing calculations.
            path (Path): Path to the root directory where calculations are
                stored.
            parse_molecules (bool): Should all molecules (for instance, along
                all optimization trajectories) be stored?
        """

        self.db = db
        self.path = path
        self.parse_molecules = parse_molecules

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

        cursor = self.db.database[self.db.autots_data_collection].find(
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
                doc = drone.assimilate(parse_molecules=self.parse_molecules)
                docs.append(doc)
            except:
                print("Cannot parse {}".format(path.as_posix()))
                continue

        if len(docs) > 0:
            self.db.update_autots_data_docs(docs, key="path")

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

