# coding: utf-8

import os
from pathlib import Path

from typing import Optional, Union


def store_reaction_template(path: Optional[Union[str, Path]] = None,
                            schrodinger_dir: Optional[str] = "$SCHRODINGER",
                            input_file: Optional[str] = None,
                            output_file: Optional[str] = None,
                            assume_irc_connections: Optional[bool] = False):
    """
    A helper function to be used with the Schrodinger store_reaction_template
        utility.

    Args:
        path (str): A root directory to search for path files. Default is None,
            meaning that the current working directory will be used.
        schrodinger_dir (str): A path to the Schrodinger Suite of software.
            This is used to call AutoTS and other utilities. By default,
            this is "$SCHRODINGER", which should be an environment variable
            set at the time of installation.
        input_file (str): Input file to search for in the path. If this is
            not provided (default is None), then the store_reaction_template
            utility will search for any files with the suffix "_full_path.mae"
        output_file (str): Path to a *.mae database of reaction templates. If
            this is not provided, then the templates will be added to
            $HOME/.schrodinger/autots_templates/autots_templates.mae
        assume_irc_connections (bool): If True (default False), then
            it will be assumed not just that the TS in the full path file
            exists, but that it links the reactants to the products also listed
            in that file. This should only be used if the user is absolutely
            certain that this is the case.

    Returns:
        None
    """

    cwd = Path.cwd()

    if path is not None:
        if isinstance(path, Path):
            os.chdir(path.as_posix())
        else:
            os.chdir(path)

    command = [schrodinger_dir + "/utilities/store_reaction_template"]
    if input_file is not None:
        command.append("-input_file")
        command.append(input_file)
    if output_file is not None:
        command.append("-output_file")
        command.append(output_file)
    if assume_irc_connections:
        command.append("-assume_irc_connections")

    os.system(" ".join(command))

    if path is not None:
        os.chdir(cwd.as_posix())

