
# coding: utf-8

import os
from typing import Optional, Dict
import datetime

from pymongo import UpdateOne

from pymatgen.analysis.graphs import MoleculeGraph

from mpcat.aggregate.database import CatDB
from mpcat.automate.generate_calcs import AutoTSJob, launch_mass_jobs


def launch_jobs_from_queue(database: CatDB,
                           base_dir: str,
                           num_launches: int = 1,
                           by_priority: bool = False,
                           query: Optional[Dict] = None,
                           schrodinger_dir: Optional[str] = "SCHRODINGER",
                           num_cores: Optional[int] = 40,
                           host: Optional[str] = "localhost",
                           save_scratch: Optional[bool] = False,
                           command_line_args: Optional[Dict] = None):
    """
    Use a CatDB.queue_collection to automatically submit jobs to a job manager.

    Args:
        database (CatDB): Database to query and identify calculations to run
        base_dir (str): Base directory in which to launch the calculations
        num_launches (int): Number of jobs to launch
        by_priority (bool): If True (default False), then prioritize jobs with
            higher priority, and do not run jobs with no priority assigned
        query (dict): MongoDB query dictionary.
        schrodinger_dir (str): A path to the Schrodinger Suite of software.
            This is used to call AutoTS and other utilities. By default,
            this is "$SCHRODINGER", which should be an environment variable
            set at the time of installation.
        num_cores (int): How many cores should the program be parallelized
            over (default 40). When multiple subjobs need to be run
            simultaneously, AutoTS will distribute these cores automatically
            between subjobs
        host (str): Which host should the calculation be run on? By default,
            this is "localhost", which should generally mean that the
            calculation is run on the current node without using a queueing
            system
        save_scratch (bool): If True (default False), save a *.zip file
            containing the contents of the calculation scratch directory
        command_line_args (dict): A dictionary of flag: value pairs to be
            provided to the autots command-line interface. Ex:
            {"WAIT": None,
             "subdir": None,
             "nsubjobs": 5}
             will be interpreted as "-WAIT -subdir -nsubjobs 5"
    """

    # If query is None, just grab all possible calculations
    if query is None:
        initial_query = [e for e in database.queue_cllection.find({"state": "READY"}, {"_id": 0,
                                                                                       "rxnid": 1,
                                                                                       "priority": 1})]
    else:
        if "state" not in query:
            query["state"] = "READY"
        initial_query = [e for e in database.queue_cllection.find(query,
                                                                  {"_id": 0,
                                                                   "rxnid": 1,
                                                                   "priority": 1})]

    if by_priority:
        initial_query = sorted([i for i in initial_query if i["priority"] is not None],
                               key=lambda x: x["priority"],
                               reverse=True)

    if num_launches > len(initial_query):
        raise ValueError("num_launches too high! Only {} jobs available".format(len(initial_query)))

    rxn_ids = [r["rxnid"] for r in initial_query[0:num_launches]]
    to_calculate = list(database.queue_cllection.find({"rxnid": {"$in": rxn_ids}}))

    requests = list()

    for calc in to_calculate:
        reactants = [MoleculeGraph.from_dict(r) for r in calc["reactants"]]
        products = [MoleculeGraph.from_dict(p) for p in calc["products"]]
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        name = "_".join("launcher", calc["rxnid"], timestamp)
        job = AutoTSJob(reactants, products,
                        os.path.join(base_dir, name),
                        schrodinger_dir=schrodinger_dir,
                        num_cores=num_cores,
                        host=host,
                        save_scratch=save_scratch,
                        input_params=calc["input"])

        requests.append(UpdateOne({"rxnid": calc["rxnid"]}, {"state": "SUBMITTED"}, upsert=True))

        job.setup_calculation()
        job.run(command_line_args=command_line_args)

    if len(requests) > 0:
        database.queue_cllection.bulk_write(requests, ordered=False)