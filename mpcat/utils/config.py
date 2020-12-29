import os
from pathlib import Path

from monty.serialization import loadfn
from mpcat.aggregate.database import CatDB


def load_from_config():
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

    return config, db
