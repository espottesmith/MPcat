#!/usr/bin/env python

import argparse
import textwrap

from mpcat.aggregate.database import CatDB
from mpcat.automate.generate_calcs import launch_jobs_from_queue


def main():
    parser = argparse.ArgumentParser(
        description="Infinite submission script for AutoTS jobs using MPCat",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
                Example of use:
                --------------------------------
                launch.py -d PATH_TO_DB -b BASE_LAUNCH_DIRECTORY -q '{QUERY}' -n infinite
                """)
    )

    parser.add_argument("-d", "--database", help='Path to database JSON file',  required=True)
    parser.add_argument("-b", "--base_dir", help="Base launch directory", required=True)
    parser.add_argument("-n", "--num_calcs", help='Number of jobs to launch', required=False)
    parser.add_argument("-c", "--num_cores", type=int, default=40, help="Number of cores for each job", required=False)
    parser.add_argument("-q", "--query", help="Database query using MongoDB format", required=False)
    parser.add_argument("-H", "--host", help="Hostname for calculations", required=False)
    parser.add_argument("-P", "--use_priority", help="Flag to choose jobs by priority (default False)", required=False, action="store_true")

    args = parser.parse_args()

    db = CatDB.from_db_file(args.database)

    if args.num.lower() == "infinite":
        while True:
            try:
                launch_jobs_from_queue(db, args.base_dir,
                                       num_launches=1,
                                       by_priority=args.use_priority,
                                       query=eval(args.query),
                                       num_cores=args.num_cores,
                                       host=args.host,
                                       command_line_args={"nsubjobs": 1, "WAIT": None})
            except ValueError:
                break
    else:
        launch_jobs_from_queue(db, args.base_dir,
                               num_launches=args.num_calcs,
                               by_priority=args.use_priority,
                               query=eval(args.query),
                               num_cores=args.num_cores,
                               host=args.host,
                               command_line_args={"nsubjobs": 1, "WAIT": None})


if __name__ == "__main__":
    main()
