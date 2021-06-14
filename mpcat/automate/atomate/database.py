# coding: utf-8


# This module defines the database classes.

from atomate.utils.database import CalcDb
from atomate.utils.utils import get_logger

__author__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__version__ = "0.1"
__credits__ = "Kiran Mathew, Anubhav Jain, Brandon Wood, Samuel Blau"
__status__ = "Alpha"
__date__ = "02/17/20"

logger = get_logger(__name__)


class GSMCalcDb(CalcDb):
    """
    Class to help manage database insertions of GSM drones
    """

    def __init__(self,
                 host="localhost",
                 port=27017,
                 database="gsm",
                 collection="gsm_tasks",
                 user=None,
                 password=None,
                 **kwargs):
        super(GSMCalcDb, self).__init__(host, port, database, collection,
                                        user, password, **kwargs)

    def build_indexes(self, indexes=None, background=True):
        """
        Build the indexes.

        Args:
            indexes (list): list of single field indexes to be built.
            background (bool): Run in the background or not.

        TODO: make sure that the index building is sensible and check for
            existing indexes.
        """
        _indices = indexes or [
            "formula_pretty", "dir_name", "last_updated"
        ]
        self.collection.create_index(
            "task_id", unique=True, background=background)
        # build single field indexes
        for ii in _indices:
            self.collection.create_index(ii, background=background)
        # TODO: build compound indexes

    def reset(self):
        self.collection.delete_many({})
        self.db.counter.delete_one({"_id": "taskid"})
        self.db.counter.insert_one({"_id": "taskid", "c": 0})
        self.build_indexes()
