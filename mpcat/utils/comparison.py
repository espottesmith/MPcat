from typing import Dict, Optional, List

import numpy as np


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

    for key in a:
        if key in exclude:
            continue
        if isinstance(a[key], dict) and isinstance(b[key], dict):
            if not all_equal(a[key], b[key]):
                print(key, a[key], b[key])
                return False
        elif isinstance(a[key], list) and isinstance(b[key], list):
            if len(a[key]) != len(b[key]):
                print(key, a[key], b[key])
                return False
            for pair in zip(a[key], b[key]):
                if not all_equal({"a": pair[0]}, {"a": pair[1]}):
                    print(key, pair[0], pair[1])
                    return False
        elif isinstance(a[key], np.ndarray) and isinstance(b[key], np.ndarray):
            if not np.allclose(a[key], b[key]):
                return False
        elif isinstance(a[key], np.ndarray):
            if not np.allclose(a[key], np.array(b[key])):
                return False
        elif isinstance(b[key], np.ndarray):
            if not np.allclose(np.array(a[key]), b[key]):
                return False
        else:
            try:
                if a[key] != b[key]:
                    print(key, a[key], b[key])
                    return False
            except ValueError:
                print(key, a[key], b[key])
                return False

    return True

