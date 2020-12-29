#!/usr/bin/env python

import datetime

from mpcat.utils.config import load_from_config


def main():
    config, db = load_from_config()

    completed = [e["rxnid"] for e in db.database[db.data_collection].find({"completed": True},
                                                                          {"_id": 0, "rxnid": 1})]
    failed = [e["rxnid"] for e in db.database[db.data_collection].find({"completed": False},
                                                                       {"_id": 0, "rxnid": 1})]

    result_comp = db.database[db.queue_collection].update_many({"rxnid": {"$in": completed},
                                                                "state": {"$in": ["READY", "SUBMITTED"]}},
                                                               {"$set": {"state": "COMPLETED",
                                                                         "updated_on": datetime.datetime.now(
                                                                             datetime.timezone.utc)
                                                                         }})

    result_fail = db.database[db.queue_collection].update_many({"rxnid": {"$in": failed},
                                                                "state": {"$in": ["READY", "SUBMITTED"]}},
                                                               {"$set": {"state": "FAILED",
                                                                         "updated_on": datetime.datetime.now(
                                                                             datetime.timezone.utc)
                                                                         }})

    print("UPDATED: {} SUCCESSFUL, {} FAILED".format(result_comp.modified_count,
                                                     result_fail.modified_count))


if __name__ == "__main__":
    main()
