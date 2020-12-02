# coding: utf-8

import datetime
from typing import Optional, Dict, List, Union
from pathlib import Path

from pymongo import MongoClient, ReturnDocument, ReplaceOne

from monty.json import jsanitize
from monty.serialization import loadfn

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph

from mpcat.utils.generate import mol_to_mol_graph
from mpcat.utils.reaction import (get_reaction_graphs,
                                  union_molgraph)


class CatDB:
    """
    A database to store the results of Jaguar and AutoTS calculations

    This class takes heavily from atomate.utils.database.CalcDb.
    Credit is due to Kiran Mathew and Anubhav Jain.
    """

    def __init__(self, host: str, port: int, database_name: str, user: str,
                 password: str, data_collection: Optional[str] = "mpcat_data",
                 queue_collection: Optional[str] = "mpcat_queue",
                 **kwargs):
        """

        Args:
            host (str): The URL of the host of the database
            port (int): The port number for the database
            database_name (str): The name of the database
            user (str): The name of the user for this connection to the database
            password (str): The password for user
            data_collection (str): The name of the collection in which to store
                tasks. Default is "mpcat_tasks"
            queue_collection (str): The name of the collection in which to store
                information about calculations to be run. Default is "mpcat_queue"
            kwargs: additional keyword arguments
        """

        self.host = host
        self.port = int(port)
        self.database_name = database_name
        self.user = user
        self.password = password
        self.data_collection = data_collection
        self.queue_collection = queue_collection

        try:
            self.client = MongoClient(host=self.host, port=self.port,
                                      username=self.user,
                                      password=self.password,
                                      **kwargs)

            self.database = self.client[self.database_name]
        except:
            raise RuntimeError("Cannot connect to DB!")

        try:
            self.database.authenticate(self.user, self.password,
                                       source=kwargs.get("authsource", None))
        except:
            raise RuntimeError("Cannot authenticate user!")

        if self.database["counter"].find({"_id": "datumid"}).count() == 0:
            self.database["counter"].insert_one({"_id": "datumid", "c": 0})
        if self.database["counter"].find({"_id": "rxnid"}).count() == 0:
            self.database["counter"].insert_one({"_id": "rxnid", "c": 0})

    def insert_data_doc(self, doc: Dict):
        """
        Insert a data doc into the database.

        Args:
            doc (dict): The task document to be inserted

        Returns:
            None
        """

        previous = self.database[self.data_collection].find_one({"path": doc["path"]},
                                                                ["path", "datumid"])

        doc["last_updated"] = datetime.datetime.now()

        if previous is None:
            if not doc.get("datumid", None):
                doc["datumid"] = self.database["counter"].find_one_and_update({"_id": "datumid"},
                                                                              {"$inc": {"c": 1}},
                                                                              return_document=ReturnDocument.AFTER)["c"]
        else:
            doc["datumid"] = previous["datumid"]

        doc = jsanitize(doc, allow_bson=True)
        self.database[self.data_collection].update_one({"path": doc["path"]},
                                                       {"$set": doc}, upsert=True)

    def update_data_doc(self, docs: List[Dict], key: Optional[str] = "path"):
        """
        Update a number of data docs at once.

        Args:
            docs (list of dicts): Task docs to be updated.
            key (str): Database key to query

        Returns:
            None
        """

        if not isinstance(docs, list):
            docs = [docs]

        requests = list()

        for doc in docs:
            doc = jsanitize(doc, allow_bson=True)

            requests.append(ReplaceOne({key: doc[key]}, doc, upsert=True))

        if len(requests) > 0:
            self.database[self.data_collection].bulk_write(requests,
                                                           ordered=False)

    def insert_calculation(self, reactants: Union[List[Molecule], List[MoleculeGraph]],
                           products: Union[List[Molecule], List[MoleculeGraph]],
                           name: Optional[str] = None,
                           input_params: Optional[Dict] = None,
                           tags: Optional[Dict] = None,
                           priority: Optional[int] = None,
                           include_reaction_graph: Optional[bool] = False):
        """
        Add a reaction to the "queue" (self.queue_collection collection).

        Args:
            reactants (list of Molecule objects): The reactants of the reaction.
                Can be separated molecules or a reaction complex
            products (list of Molecule objects): The products of the reaction.
                Can be separated molecules or a reaction complex
            name (str, or None): Name of the reaction. No
            input_params (Dict, or None): Dictionary with all input parameters
                for this calculation. These keywords and the associated values
                will be provided to AutoTSSet (or, eventually, JaguarSet).
            tags (Dict, or None): Dictionary with some calculation metadata
                Ex: {"class": "production", "time": 3}
            priority (int, or None): To identify jobs that should or should
                not be run, calculations can be prioritized. The higher the
                priority, the more important the calculation. If the priority is
                None (default), then the job will not be selected unless
                chosen specifically by ID or other query. If the number is
                negative (< 0), the calculation will never be selected.
            include_reaction_graph (bool): Should a reaction graph be generated
                from the reactant and product MoleculeGraphs? This might be
                skipped because it can be costly to perform subgraph isomorphisms
                and identify the appropriate reaction graph.

        Returns:
            None
        """

        entry = {"state": "READY"}

        if len(reactants) == 0 or len(products) == 0:
            raise ValueError("reactants and products must be non-empty lists!")

        entry["reactants"] = list()
        entry["products"] = list()

        for reactant in reactants:
            entry["reactants"].append(mol_to_mol_graph(reactant))

        for product in products:
            entry["products"].append(mol_to_mol_graph(product))

        rct_charge = sum([m.molecule.charge for m in entry["reactants"]])
        pro_charge = sum([m.molecule.charge for m in entry["products"]])

        if rct_charge != pro_charge:
            raise ValueError("Reactants and products do not have the same charge!")

        entry["charge"] = rct_charge

        rct_nelectrons = sum([m.molecule._nelectrons for m in entry["reactants"]])
        pro_nelectrons = sum([m.molecule._nelectrons for m in entry["products"]])

        if rct_nelectrons != pro_nelectrons:
            raise ValueError("Reactants and products do not have the same number of electrons!")

        entry["nelectrons"] = int(rct_nelectrons)
        entry["spin_multiplicity"] = 1 if entry["nelectrons"] % 2 == 0 else 2

        if name is None:
            rct_names = [m.molecule.composition.alphabetical_formula + "_" + str(m.molecule.charge)
                         for m in entry["reactants"]]
            pro_names = [m.molecule.composition.alphabetical_formula + "_" + str(m.molecule.charge)
                         for m in entry["products"]]
            entry["name"] = " + ".join(rct_names) + " -> " + " + ".join(pro_names)
        else:
            entry["name"] = name

        entry["input"] = input_params
        entry["tags"] = tags
        entry["priority"] = priority

        union_rct = union_molgraph(entry["reactants"])
        union_pro = union_molgraph(entry["products"])

        entry["molgraph_rct"] = union_rct.as_dict()
        entry["molgraph_pro"] = union_pro.as_dict()

        if include_reaction_graph:
            reaction_graph = get_reaction_graphs(union_rct, union_pro,
                                                 allowed_form=2, allowed_break=2,
                                                 stop_at_one=True)
            if len(reaction_graph) == 0:
                raise RuntimeError("No valid reaction could be found between "
                                   "reactants and products!")

            entry["reaction_graph"] = reaction_graph[0].as_dict()
        else:
            entry["reaction_graph"] = None

        entry["reactants"] = [r.as_dict() for r in entry["reactants"]]
        entry["products"] = [p.as_dict() for p in entry["products"]]

        entry["rxnid"] = self.database["counter"].find_one_and_update({"_id": "rxnid"},
                                                                      {"$inc": {"c": 1}},
                                                                      return_document=ReturnDocument.AFTER)["c"]
        entry["created_on"] = datetime.datetime.now(datetime.timezone.utc)
        entry["updated_on"] = datetime.datetime.now(datetime.timezone.utc)

        doc = jsanitize(entry, allow_bson=True)
        self.database[self.queue_collection].update_one({"rxnid": doc["rxnid"]},
                                                        {"$set": doc}, upsert=True)

    def sync_queue_data(self):
        # TODO: Implement this
        pass

    @classmethod
    def from_db_file(cls, db_file: Union[str, Path],
                     admin: Optional[bool]=True):
        """
        Create MMDB from database file. File requires host, port, database,
        collection, and optionally admin_user/readonly_user and
        admin_password/readonly_password

        This method is taken directly from atomate. Credit is due to Kiran
        Mathew and Anubhav Jain.

        Args:
            db_file (str): path to the file containing the credentials
            admin (bool): whether to use the admin user

        Returns:
            MMDb object
        """

        if isinstance(db_file, Path):
            creds = loadfn(db_file.as_posix())
        else:
            creds = loadfn(db_file)

        if admin and "admin_user" not in creds and "readonly_user" in creds:
            raise ValueError("Trying to use admin credentials, "
                             "but no admin credentials are defined. "
                             "Use admin=False if only read_only "
                             "credentials are available.")

        if admin:
            user = creds.get("admin_user")
            password = creds.get("admin_password")
        else:
            user = creds.get("readonly_user")
            password = creds.get("readonly_password")

        kwargs = creds.get("mongoclient_kwargs", {})  # any other MongoClient kwargs can go here ...

        if "authsource" in creds:
            kwargs["authsource"] = creds["authsource"]
        else:
            kwargs["authsource"] = creds["database"]

        return cls(creds["host"], int(creds.get("port", 27017)),
                   creds["database"],
                   user, password,
                   data_collection=creds.get("data_collection", "mpcat_data"),
                   queue_collection=creds.get("queue_collection", "mpcat_queue"),
                   **kwargs)
