import copy
from itertools import combinations, product
from typing import Optional, List
from difflib import SequenceMatcher

import numpy as np

import networkx as nx
import networkx.algorithms.isomorphism as iso

from pymatgen.analysis.graphs import MoleculeGraph


#TODO: Should there actually be a ReactionGraph class that subclasses MoleculeGraph
# This class could include some additional data, like a charge_change property,
# or some information about the charges of the fragments (if there are multiple
# disconnected subgraphs

def get_reaction_graphs(mol1: MoleculeGraph,
                        mol2: MoleculeGraph,
                        allowed_form: Optional[int] = 1,
                        allowed_break: Optional[int] = 1,
                        stop_at_one: Optional[bool] = True):
    """
    Generate "reaction graphs", MoleculeGraphs that indicate which bonds are
        formed and which bonds are broken in the course of a chemical reaction.

    Args:
        mol1 (MoleculeGraph): First molecule to be compared
        mol2 (MoleculeGraph): Second molecule to be compared
        allowed_form (int): How many bonds, at most, can be added to create an
            isomorphism between the two MoleculeGraphs? Default is 1
        allowed_break (int): How many bonds, at most, can be removed to create
            an isomorphism between the two MoleculeGraphs? Default is 1
        stop_at_one (bool): If True (default), then this function will stop
            after finding one appropriate reaction graph. Otherwise, all valid
            reaction graphs will be returned

        Note: because this function involves creating a combinatorial number
        of MoleculeGraphs, it can be quite slow. For best performance, use with
        small molecules and only allow a single possible reaction graph to be
        constructed (use stop_at_one=True)

    :return:
        reaction_graphs (List of MoleculeGraph objects)

    """

    reaction_graphs = list()

    if mol1.isomorphic_to(mol2):
        return [copy.deepcopy(mol1)]
    elif allowed_form < 1 and allowed_break < 1:
        return list()

    species_1 = nx.get_node_attributes(mol1.graph, "specie")
    species_2 = nx.get_node_attributes(mol2.graph, "specie")

    if set(species_1.values()) != set(species_2.values()):
        return list()

    diff_edges = mol2.graph.size() - mol1.graph.size()
    if diff_edges > allowed_form or -1 * diff_edges > allowed_break:
        return list()

    current_edges = list(mol1.graph.edges())
    new_edges = list()
    for i in mol1.graph:
        for j in mol1.graph:
            if j > i:
                if (i, j) not in current_edges and (j, i) not in current_edges:
                    new_edges.append((i, j))

    forming = 0
    stop_now = False
    while forming <= allowed_form and forming <= len(new_edges) and not stop_now:
        combo_form = [x for x in combinations(new_edges, forming)]
        breaking = 0
        while breaking <= allowed_break and breaking <= len(current_edges) and not stop_now:
            combo_break = [x for x in combinations(current_edges, breaking)]

            # All combinations of bond changes accessible
            # using "forming" bond additions and "breaking" bond removals
            form_break_combos = [x for x in product(combo_form, combo_break)]
            for bonds_formed, bonds_broken in form_break_combos:
                change_copy = copy.deepcopy(mol1)

                # Add/remove all bonds
                if forming > 0:
                    for bond in bonds_formed:
                        change_copy.add_edge(bond[0], bond[1],
                                             edge_properties={"formed": True,
                                                              "broken": False})
                if breaking > 0:
                    for bond in bonds_broken:
                        change_copy.break_edge(bond[0], bond[1],
                                               allow_reverse=True)

                if change_copy.isomorphic_to(mol2):
                    # Need to re-add broken bonds, but tagged appropriately
                    for bond in bonds_broken:
                        change_copy.add_edge(bond[0], bond[1],
                                             edge_properties={"formed": False,
                                                              "broken": True})

                    reaction_graphs.append(change_copy)
                    if stop_at_one:
                        stop_now = True
                        break

            breaking += 1
        forming += 1

    return reaction_graphs


def union_molgraph(mol_graphs: List[MoleculeGraph],
                   validate_proximity: Optional[bool] = False):
    """
    Combine arbitrarily many MoleculeGraphs into one.

    This function might be useful to, for instance, combine multiple reactants
    or products to form complexes. Note that this function does NOT form any
    bonds between the various MoleculeGraphs to be united; this has to be done
    by the user.

    Args:
        mol_graphs (List of MoleculeGraph objects)
        validate_proximity (bool): If True (default False), then when new nodes
            are added, an error will be raised if they are too close. Problems
            with different molecules being too close to one another can be
            resolved by, for instance, offsetting the centers of mass of the
            molecules.
    Return
        combined (MoleculeGraph)
    """

    if len(mol_graphs) == 0:
        raise ValueError("No MoleculeGraph provided!")
    elif len(mol_graphs) == 1:
        return copy.deepcopy(mol_graphs[0])

    combined = copy.deepcopy(mol_graphs[0])
    for mg in mol_graphs[1:]:
        mapping = dict()
        for s, site in enumerate(mg.molecule):
            mapping[s] = len(combined.molecule)
            combined.insert_node(len(combined.molecule),
                                 site.specie,
                                 site.coords,
                                 validate_proximity=validate_proximity,
                                 site_properties=site.properties)
        for edge in mg.graph.edges(data=True):
            combined.add_edge(mapping[edge[0]],
                              mapping[edge[1]],
                              edge_properties=edge[2])

    total_charge = sum([e.molecule.charge for e in mol_graphs])
    combined.molecule.set_charge_and_spin(total_charge)

    combined.set_node_attributes()
    return combined


def get_atom_mappings(mol1: MoleculeGraph,
                      mol2: MoleculeGraph,
                      allowed_form: Optional[int] = 1,
                      allowed_break: Optional[int] = 1,
                      stop_at_one: Optional[bool] = True,
                      give_best: Optional[bool] = False):
    """
    Create atom mappings between two MoleculeGraphs

    Args:
        mol1 (MoleculeGraph): First molecule to be compared
        mol2 (MoleculeGraph): Second molecule to be compared
        allowed_form (int): How many bonds, at most, can be added to create an
            isomorphism between the two MoleculeGraphs? Default is 1
        allowed_break (int): How many bonds, at most, can be removed to create
            an isomorphism between the two MoleculeGraphs? Default is 1
        stop_at_one (bool): If True (default), then this function will stop
            after finding one appropriate mapping. Otherwise, all valid
            mappings will be returned
        give_best (bool): If True (default False), only return the "best"
            mapping, as defined by the distances between atoms.

    Returns:
        mappings (List of Dicts): The possible atom mappings. Note that the keys
            of such a dict are the atom indices in mol1, and the values are the
            atom indices in mol2
    """

    reaction_graphs = get_reaction_graphs(mol1, mol2,
                                          allowed_form=allowed_form,
                                          allowed_break=allowed_break,
                                          stop_at_one=stop_at_one)

    if len(reaction_graphs) == 0:
        return list()

    nm = iso.categorical_node_match("specie", "ERROR")
    gr2 = mol2.graph.to_undirected()

    mappings = list()
    for rg in reaction_graphs:
        # Break all necessary bonds to create the isomorphism
        for atom_1, atom_2, data in rg.graph.edges(data=True):
            if data.get("broken", False):
                rg.break_edge(atom_1, atom_2, allow_reverse=True)

        gr1 = rg.graph.to_undirected()

        matcher = iso.GraphMatcher(gr1, gr2, node_match=nm)

        # Sanity check; after this, rg should always be isomorphic to mol2
        if not matcher.is_isomorphic():
            raise ValueError("Problem with reaction graphs! No isomorphism "
                             "found!")

        for isomorphism in matcher.isomorphisms_iter():
            mappings.append(isomorphism)
            if stop_at_one:
                return mappings

    if give_best:
        # Ideal matching would have very similar atom distances
        d1 = get_ranked_atom_dists(mol1.molecule)
        d2 = get_ranked_atom_dists(mol2.molecule)

        ratios = list()
        dist_seq_2 = np.array([v for k, v in sorted(list(d2.items()), key=lambda x: x[0])]).flatten()
        for mapping in mappings:
            dist_seq_1 = np.array([v for k, v in sorted(list(d1.items()), key=lambda x: mapping[x[0]])]).flatten()
            matcher = SequenceMatcher(None, dist_seq_1, dist_seq_2)
            ratios.append(matcher.ratio())

        min_ratio = min(ratios)
        min_index = ratios.index(min_ratio)
        return [mappings[min_index]]
    else:
        return mappings


def get_ranked_atom_dists(mol):
    """

    :param mol:
    :return:
    """
    dist_matrix = mol.distance_matrix

    result = dict()
    for num, row in enumerate(dist_matrix):
        ranking = np.argsort(row)
        # The first member will always be the atom itself, which should be excluded
        result[num] = ranking[1:]
    return result


def get_reaction_template(reaction_graph: MoleculeGraph,
                          order: Optional[int] = 0):
    """
    Generate a reaction template graph based on the bonds that break and form
        in a reaction.

    Args:
        reaction_graph (MoleculeGraph): A MoleculeGraph representing a reaction,
            such as those produced by get_reaction_graphs. Edges that are broken
            in the reaction must have the property "broken" = True, and edges
            that are formed must have the property "formed" = True
        order (int): This integer indicates what order of neighbors are added
            to the template. Order 0 (default) indicates that only the atoms
            involved in the reaction will be included. Order 1 includes the
            nearest neighbors of those atoms, and so on.

    Returns:
        template (networkx.Graph)
    """

    bonds = list()
    for ind_from, ind_to, data in reaction_graph.graph.edges(data=True):
        if data.get("formed", False) or data.get("broken", False):
            bonds.append((ind_from, ind_to))

    indices = reaction_graph.extract_bond_environment(bonds, order=order)
    template = reaction_graph.graph.subgraph(list(indices)).copy().to_undirected()
    return template

