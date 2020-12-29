#!/usr/bin/env python

from pathlib import Path
import datetime

from mpcat.utils.config import load_from_config


def main():

    config, db = load_from_config()

    base_dir = Path(config["base_dir"])

    to_restart = list()

    for calc_dir in [x for x in base_dir.iterdir() if x.is_dir()]:
        files = [e.name for e in calc_dir.iterdir() if e.is_file()]

        if "autots.in" in files and "AutoTS.T9XnCsLi.out" not in files:
            rxn_id = None
            if calc_dir.name.startswith("launcher"):
                rxn_id = int(calc_dir.name.split("_")[1])
                to_restart.append(rxn_id)
            if rxn_id is None:
                print("Could not determine rxnid for", calc_dir.name)

    result = db.database[db.queue_collection].update_many({"rxnid": {"$in": to_restart},
                                                           "state": "SUBMITTED"},
                                                          {"$set": {"state": "READY",
                                                                    "updated_on": datetime.datetime.now(
                                                                        datetime.timezone.utc
                                                                    )}})
    print("UPDATED {} ENTRIES".format(result.modified_count))


if __name__ == "__main__":
    main()
