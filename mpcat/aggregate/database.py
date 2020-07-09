# coding: utf-8

import datetime
from typing import Optional, Dict, List

from pymongo import MongoClient, ReturnDocument, ReplaceOne, UpdateOne

from monty.json import jsanitize
from monty.serialization import loadfn


class CatDB:
    """
    A database to store the results of Jaguar and AutoTS calculations

    This class takes heavily from atomate.utils.database.CalcDb.
    Credit is due to Kiran Mathew and Anubhav Jain.
    """

    def __init__(self, host: str, port: int, database_name: str, user: str,
                 password: str, collection: Optional[str] = "mpcat"):
        """

        Args:
            host (str): The URL of the host of the database
            port (int): The port number for the database
            database_name (str): The name of the database
            user (str): The name of the user for this connection to the database
            password (str): The password for user
            collection (str): The name of the collection in which to store
                tasks. Default is "mpcat"
        """

        self.host = host
        self.port = port
        self.database_name = database_name
        self.user = user
        self.password = password
        self.collection = collection

        try:
            self.client = MongoClient(host=self.host, port=self.port,
                                      username=self.user,
                                      password=self.password)

            self.database = self.client[self.database_name]
        except:
            raise RuntimeError("Cannot connect to DB!")

        try:
            self.database.authenticate(self.user, self.password)
        except:
            raise RuntimeError("Cannot authenticate user!")

        if self.database["counter"].find({"_id": "taskid"}).count() == 0:
            self.database["counter"].insert_one({"_id": "taskid", "c": 0})

    def insert(self, doc: Dict):
        """
        Insert a task doc into the database.

        Args:
            doc (dict): The task document to be inserted

        Returns:
            None
        """

        previous = self.database[self.collection].find_one({"path": doc["path"]},
                                                           ["path", "task_id"])

        doc["last_updated"] = datetime.datetime.now()

        if previous is None:
            if not doc.get("task_id", None):
                doc["task_id"] = self.database["counter"].find_one_and_update({"_id": "taskid"},
                                                                              {"$inc": {"c": 1}},
                                                                              return_document=ReturnDocument.AFTER)["c"]
        else:
            doc["task_id"] = previous["task_id"]

        doc = jsanitize(doc, allow_bson=True)
        self.database[self.collection].update_one({"path": doc["path"]},
                                                  {"$set": doc}, upsert=True)

    def update(self, docs: List[Dict], key: Optional[str] = "path"):
        """
        Update a number of docs at once.

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
            self.database[self.collection].bulk_write(requests, ordered=False)

    @classmethod
    def from_db_file(cls, db_file, admin=True):
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

        return cls(creds["host"], int(creds.get("port", 27017)),
                   creds["database"],
                   user, password,
                   creds.get("collection", "mpcat"))

