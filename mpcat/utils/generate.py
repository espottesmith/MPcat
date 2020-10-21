# coding: utf-8

from typing import Optional, List, Dict, Type, Union

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.analysis.fragmenter import metal_edge_extender


def mol_to_mol_graph(molecule: Union[Molecule, MoleculeGraph]):
    """
    Convert a Molecule to a MoleculeGraph using a default connectivity
    algorithm.

    Args:
        molecule (Molecule): Molecule to be converted

    Returns:
        mol_graph: MoleculeGraph
    """

    if isinstance(molecule, MoleculeGraph):
        return molecule
    else:
        mol_graph = MoleculeGraph.with_local_env_strategy(molecule, OpenBabelNN())
        return metal_edge_extender(mol_graph)
