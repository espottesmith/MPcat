from typing import Dict, Optional, Union

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


def reaction_path_successful_ts(
        db: CatDB,
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


def optimize_failed_ts(
        db: CatDB,
        query: Optional[Dict] = None,
        num_run: Optional[int] = None,
        with_critic: Optional[bool] = False,
        qchem_cmd: Optional[str] = ">>qchem_cmd<<",
        max_cores: Optional[Union[str, int]] = ">>max_cores<<",
        multimode: Optional[str] = ">>multimode<<",
        qchem_input_params: Optional[Dict] = None,
        db_file: Optional[str] = ">>db_file<<",
        allow_failed_calcs: Optional[bool] = False,
        tags: Optional[Dict] = None
):
    pass


def reaction_path_failed_ts(
        db: CatDB,
        query: Optional[Dict] = None,
        num_run: Optional[int] = None,
        with_critic: Optional[bool] = False,
        perturb_scale: Optional[float] = 1.0,
        qchem_cmd: Optional[str] = ">>qchem_cmd<<",
        max_cores: Optional[Union[str, int]] = ">>max_cores<<",
        multimode: Optional[str] = ">>multimode<<",
        qchem_input_params: Optional[Dict] = None,
        db_file: Optional[str] = ">>db_file<<",
        allow_failed_calcs: Optional[bool] = False,
        tags: Optional[Dict] = None
):
    pass
