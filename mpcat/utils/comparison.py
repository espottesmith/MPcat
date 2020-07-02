from typing import Dict, Optional, List

import numpy as np

from pymatgen.core.structure import Molecule


def all_equal(a: Dict, b: Dict, exclude_keys: Optional[List] = None):
    """
    Ensure that two dictionaries are equal, in a recursive fashion.

    Args:
        a (dict):
        b (dict):
        exclude_keys (list): List of dictionary keys to not consider when
            determining equivalence

    Returns:
        bool

    """

    if exclude_keys is None:
        exclude = list()
    else:
        exclude = exclude_keys

    if set(a.keys()) != set(b.keys()):
        return False

    for key in a:
        if key in exclude:
            continue
        if isinstance(a[key], dict) and isinstance(b[key], dict):
            if not all_equal(a[key], b[key]):
                return False
        elif isinstance(a[key], list) and isinstance(b[key], list):
            if len(a[key]) != len(b[key]):
                return False
            for pair in zip(a[key], b[key]):
                if not all_equal({"a": pair[0]}, {"a": pair[1]}):
                    return False
        elif isinstance(a[key], np.ndarray) and isinstance(b[key], np.ndarray):
            if not np.array_equal(a[key], b[key]):
                return False
        elif isinstance(a[key], np.ndarray) and isinstance(b[key], list):
            if not np.array_equal(a[key], np.array(b[key])):
                return False
        elif isinstance(b[key], np.ndarray) and isinstance(a[key], list):
            if not np.array_equal(np.array(a[key]), b[key]):
                return False
        elif not isinstance(a[key], type(b[key])):
            return False
        else:
            try:
                if a[key] != b[key]:
                    return False
            except ValueError:
                return False

    return True


def compositions_equal(reactants: List[Molecule],
                       products: List[Molecule]):
    """
    Compare the compositions of the reactants and products of a reaction.

    Args:
        reactants (list of Molecules):
        products (list of Molecules):

    Return:
        bool
    """

    rct_comp_dict = dict()
    pro_comp_dict = dict()

    for reactant in reactants:
        for element, value in reactant.composition.as_dict().items():
            if element in rct_comp_dict:
                rct_comp_dict[element] += value
            else:
                rct_comp_dict[element] = value

    for product in products:
        for element, value in product.composition.as_dict().items():
            if element in pro_comp_dict:
                pro_comp_dict[element] += value
            else:
                pro_comp_dict[element] = value

    return all_equal(rct_comp_dict, pro_comp_dict)