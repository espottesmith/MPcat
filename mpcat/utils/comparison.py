from typing import Dict


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
            for pair in zip(a[key], b[key]):
                if pair[0] != pair[1]:
                    print(key, a[key], b[key])
                    return False
        else:
            if a[key] != b[key]:
                print(key, a[key], b[key])
                return False

    return True