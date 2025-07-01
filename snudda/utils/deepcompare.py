import numpy as np


def deep_compare(a, b, float_tol=1e-9):
    """Recursively compares two structures with support for NumPy arrays and floating point tolerance."""

    if isinstance(a, dict) and isinstance(b, dict):
        if a.keys() != b.keys():
            return False
        return all(deep_compare(a[k], b[k], float_tol) for k in a)

    elif isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        if len(a) != len(b) or type(a) != type(b):
            return False
        return all(deep_compare(x, y, float_tol) for x, y in zip(a, b))

    elif isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
        return a.shape == b.shape and np.allclose(a, b, atol=float_tol)

    elif isinstance(a, float) and isinstance(b, float):
        return abs(a - b) <= float_tol

    else:
        return a == b
