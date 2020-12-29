#!/usr/bin/env python

from mpcat.automate.generate_calcs import launch_jobs_from_queue
from mpcat.utils.config import load_from_config


def main():
    config, db = load_from_config()

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
