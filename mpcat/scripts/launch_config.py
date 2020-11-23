#!/usr/bin/env python

import os
import pathlib

from monty.serialization import loadfn

from mpcat.aggregate.database import CatDB
from mpcat.automate.generate_calcs import launch_jobs_from_queue

#TODO: Change to use pathlib


def main():
    mpcat_config = os.getenv("MPCAT_CONFIG")
    if mpcat_config:
        config_file = os.path.join(mpcat_config, "config.json")
    else:
        config_file = "config.json"

    config = loadfn(config_file)

    if mpcat_config:
        db = CatDB.from_db_file(os.path.join(mpcat_config, config["db_file"]))
    else:
        db = CatDB.from_db_file(config["db_file"])

    if str(config["num_launches"]).lower() == "infinite":
        while True:
            try:
                launch_jobs_from_queue(db, config["base_dir"],
                                       num_launches=1,
                                       by_priority=config.get("use_priority", False),
                                       query=config.get("query"),
                                       num_cores=int(config.get("num_cores")),
                                       host=config.get("host"),
                                       command_line_args={"nsubjobs": 1, "WAIT": None})
            except ValueError:
                break
    else:
        launch_jobs_from_queue(db, config["base_dir"],
                               num_launches=int(config["num_launches"]),
                               by_priority=config.get("use_priority", False),
                               query=config.get("query"),
                               num_cores=int(config.get("num_cores")),
                               host=config.get("host"),
                               command_line_args={"nsubjobs": 1, "WAIT": None})


if __name__ == "__main__":
    main()
