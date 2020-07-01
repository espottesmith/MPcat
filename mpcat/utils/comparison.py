from typing import Dict

import numpy as np


def all_equal(a: Dict, b: Dict):
    """
    Ensure that two dictionaries are equal, in a recursive fashion.

    Args:
        a (dict):
        b (dict):

    Returns:
        bool

    """
    for key in a:
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
            if not np.array_equal(a[key], b[key]):
                return False
        else:
            if a[key] != b[key]:
                print(key, a[key], b[key])
                return False

    return True

