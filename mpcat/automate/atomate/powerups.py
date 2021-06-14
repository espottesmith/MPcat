"""THIS CODE IS TAKEN FROM atomate"""

from fireworks import Workflow

from atomate.utils.utils import get_fws_and_tasks

__author__ = "Anubhav Jain, Kiran Mathew, Alex Ganose"
__email__ = "ajain@lbl.gov, kmathew@lbl.gov"


def add_tags(original_wf, tags_list):
    """
    Adds tags to all Fireworks in the Workflow, WF metadata, as well as
    additional_fields for the VaspDrone to track them later (e.g. all fireworks
    and vasp tasks related to a research project)

    Args:
        original_wf (Workflow)
        tags_list: list of tags parameters (list of strings)

    Returns:
       Workflow
    """

    # WF metadata
    if "tags" in original_wf.metadata:
        original_wf.metadata["tags"].extend(tags_list)
    else:
        original_wf.metadata["tags"] = tags_list

    # FW metadata
    for idx_fw in range(len(original_wf.fws)):
        if "tags" in original_wf.fws[idx_fw].spec:
            original_wf.fws[idx_fw].spec["tags"].extend(tags_list)
        else:
            original_wf.fws[idx_fw].spec["tags"] = tags_list

    # DB insertion tasks
    for constraint in ["VaspToDb", "BoltztrapToDb"]:
        idxs = get_fws_and_tasks(original_wf, task_name_constraint=constraint)
        for idx_fw, idx_t in idxs:
            if "tags" in original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"]:
                original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"][
                    "tags"
                ].extend(tags_list)
            else:
                original_wf.fws[idx_fw].tasks[idx_t]["additional_fields"][
                    "tags"
                ] = tags_list

    return original_wf