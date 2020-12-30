from typing import Dict, Optional, Union, List
import datetime

from pymatgen.core.structure import Molecule

from fireworks import LaunchPad, Workflow

from atomate.qchem.database import QChemCalcDb
from atomate.qchem.fireworks.core import FrequencyFlatteningTransitionStateFW
from atomate.qchem.workflows.base.FF_and_critic import get_wf_FFTSopt_and_critic
from atomate.qchem.workflows.base.ts_search import (
    get_workflow_reaction_path,
    get_workflow_reaction_path_with_ts,
    get_workflow_reaction_path_with_ts_and_critic
)
from atomate.vasp.powerups import add_tags

from mpcat.aggregate.database import CatDB


def optimize_successful_ts(
        db: CatDB,
        lp: LaunchPad,
        query: Optional[Dict] = None,
        num_run: Optional[int] = None,
        with_critic: Optional[bool] = False,
        qchem_cmd: Optional[str] = ">>qchem_cmd<<",
        max_cores: Optional[Union[str, int]] = ">>max_cores<<",
        multimode: Optional[str] = ">>multimode<<",
        qchem_input_params: Optional[Dict] = None,
        db_file: Optional[str] = ">>db_file<<",
        tags: Optional[Dict] = None
):
    if query is None:
        success_query = {"completed": True, "run_atomate": {"$ne": True}}
    else:
        success_query = query
        success_query["completed"] = True
        success_query["run_atomate"] = {"$ne": True}

    possible_entries = [e for e in db.database[db.data_collection].find(
        success_query, {"_id": 0, "rxnid": 1, "output": 1}
    )]

    if num_run is not None:
        if num_run < len(possible_entries):
            possible_entries = possible_entries[:num_run]

    for entry in possible_entries:
        rxnid = entry["rxnid"]
        entry_names = list(entry["output"]["optimized_structure_energies"].keys())
        ts_structures = [
            Molecule.from_dict(e) for i, e in enumerate(entry["output"]["path_molecules"])
            if "transition_state" in entry_names[i]
        ]

        if with_critic:
            for i, ts in enumerate(ts_structures):
                name = "rxn_{}:ts_{}".format(rxnid, i + 1)
                wf = get_wf_FFTSopt_and_critic(ts, name,
                                               qchem_input_params=qchem_input_params,
                                               db_file=db_file)
                if tags is not None:
                    wf = add_tags(wf, tags)
                lp.add_wf(wf)

        else:
            for i, ts in enumerate(ts_structures):
                name = "rxn_{}:ts_{}".format(rxnid, i + 1)
                fw = FrequencyFlatteningTransitionStateFW(
                    molecule=ts,
                    name=name,
                    qchem_cmd=qchem_cmd,
                    multimode=multimode,
                    max_cores=max_cores,
                    qchem_input_params=qchem_input_params,
                    linked=True,
                    freq_before_opt=True,
                    db_file=db_file
                )
                wf = Workflow([fw], name=name)
                if tags is not None:
                    wf = add_tags(wf, tags)
                lp.add_wf(wf)

        time_now = datetime.datetime.now(datetime.timezone.utc)
        db.database[db.data_collection].update_one({"rxnid": rxnid}, {"$set": {"run_atomate": True,
                                                                               "updated_on": time_now}})


def optimize_failed_ts(
        db: CatDB,
        lp: LaunchPad,
        query: Optional[Dict] = None,
        num_run: Optional[int] = None,
        allowed_calc_types: Optional[frozenset] = frozenset(["relax_ts", "qst"]),
        allow_failed_calcs: Optional[bool] = False,
        with_critic: Optional[bool] = False,
        qchem_cmd: Optional[str] = ">>qchem_cmd<<",
        max_cores: Optional[Union[str, int]] = ">>max_cores<<",
        multimode: Optional[str] = ">>multimode<<",
        qchem_input_params: Optional[Dict] = None,
        db_file: Optional[str] = ">>db_file<<",
        tags: Optional[Dict] = None
):
    if query is None:
        failed_query = {"completed": False, "run_atomate": {"$ne": True}}
    else:
        failed_query = query
        failed_query["completed"] = False
        failed_query["run_atomate"] = {"$ne": True}

    possible_entries = [e for e in db.database[db.data_collection].find(
        failed_query, {"_id": 0, "rxnid": 1, "calcs": 1}
    )]

    if num_run is not None:
        if num_run < len(possible_entries):
            possible_entries = possible_entries[:num_run]

    for entry in possible_entries:
        rxnid = entry["rxnid"]

        for calc in entry["calcs"]:
            if not calc["success"] and not allow_failed_calcs:
                continue
            name = calc["job_name"]
            if not any([e in name for e in allowed_calc_types]):
                continue

            if "qst" in name:
                if calc["output"].get("frequencies", [None])[0] is None:
                    continue
                elif calc["output"]["frequencies"][0] > 0:
                    continue

            if calc["output"].get("molecule") is None:
                continue
            else:
                molecule = Molecule.from_dict(calc["output"]["molecule"])

                if calc["success"]:
                    calc_status = "success"
                else:
                    calc_status = "failed"

                wf_name = "failed_rxn_{}:{}_{}".format(
                        rxnid, calc_status, name)

                if with_critic:
                    wf = get_wf_FFTSopt_and_critic(molecule, wf_name,
                                                   qchem_input_params=qchem_input_params,
                                                   db_file=db_file)
                    if tags is not None:
                        wf = add_tags(wf, tags)
                    lp.add_wf(wf)

                else:
                    fw = FrequencyFlatteningTransitionStateFW(
                        molecule=molecule,
                        name=wf_name,
                        qchem_cmd=qchem_cmd,
                        multimode=multimode,
                        max_cores=max_cores,
                        qchem_input_params=qchem_input_params,
                        linked=True,
                        freq_before_opt=True,
                        db_file=db_file
                    )
                    wf = Workflow([fw], name=wf_name)
                    if tags is not None:
                        wf = add_tags(wf, tags)
                    lp.add_wf(wf)

        time_now = datetime.datetime.now(datetime.timezone.utc)
        db.database[db.data_collection].update_one({"rxnid": rxnid}, {"$set": {"run_atomate": True,
                                                                               "updated_on": time_now}})


def reaction_path_successful_ts(
        db: CatDB,
        lp: LaunchPad,
        query: Optional[Dict] = None,
        num_run: Optional[int] = None,
        with_critic: Optional[bool] = False,
        perturb_scale: Optional[float] = 1.0,
        qchem_cmd: Optional[str] = ">>qchem_cmd<<",
        max_cores: Optional[Union[str, int]] = ">>max_cores<<",
        multimode: Optional[str] = ">>multimode<<",
        qchem_input_params: Optional[Dict] = None,
        db_file: Optional[str] = ">>db_file<<",
        tags: Optional[Dict] = None
):
    pass


def reaction_path_failed_ts(
        db: CatDB,
        lp: LaunchPad,
        query: Optional[Dict] = None,
        num_run: Optional[int] = None,
        allow_failed_calcs: Optional[bool] = False,
        with_critic: Optional[bool] = False,
        perturb_scale: Optional[float] = 1.0,
        qchem_cmd: Optional[str] = ">>qchem_cmd<<",
        max_cores: Optional[Union[str, int]] = ">>max_cores<<",
        multimode: Optional[str] = ">>multimode<<",
        qchem_input_params: Optional[Dict] = None,
        db_file: Optional[str] = ">>db_file<<",
        tags: Optional[Dict] = None
):
    pass
