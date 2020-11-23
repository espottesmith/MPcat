#!/usr/bin/env python

import os
from pathlib import Path

from monty.serialization import loadfn

from mpcat.aggregate.database import CatDB


def main():
    mpcat_config = os.getenv("MPCAT_CONFIG")
    if mpcat_config:
        config_file = Path(mpcat_config) / "config.json"
    else:
        config_file = Path("config.json")

    config = loadfn(config_file.as_posix())

    if mpcat_config:
        db = CatDB.from_db_file(Path(mpcat_config) / config["db_file"])
    else:
        db = CatDB.from_db_file(config["db_file"])

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

    result = db.database[db.queue_collection].update_many({"rxnid": {"$in": to_restart}},
                                                          {"$set": {"state": "READY"}})
    print("Updated {} entries".format(result.modified_count))


if __name__ == "__main__":
    main()
