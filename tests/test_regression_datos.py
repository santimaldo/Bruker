"""
Regression tests for the legacy ``Datos`` class family.

Every case in :mod:`regression_cases` is re-run and compared, array by
array, with the reference stored in ``tests/golden/<case>.npz`` by
``tests/make_golden.py``. A failure means the refactored code no longer
reproduces the numbers the analysis scripts produced before.

Cases whose dataset is not on this machine are skipped, not failed.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from make_golden import golden_file  # noqa: E402
from regression_cases import CASES, missing_data  # noqa: E402

#: Relative tolerance for every comparison. The refactor must reproduce the
#: old numbers to floating-point noise, not just "approximately".
RTOL = 1e-9


def assert_matches_golden(name, result):
    """
    Compare a freshly computed case result with its stored reference.

    Parameters
    ----------
    name : str
        Case name; selects ``tests/golden/<name>.npz``.
    result : dict of str to array_like
        Output of the case function.

    Raises
    ------
    AssertionError
        If keys differ, shapes differ, or any value differs by more than
        :data:`RTOL` (relative to the largest magnitude in that array).
    """
    with np.load(golden_file(name)) as golden:
        expected = {k: golden[k] for k in golden.files if not k.startswith("_meta.")}
    assert set(result) == set(expected), (
        f"keys changed: added {sorted(set(result) - set(expected))}, "
        f"removed {sorted(set(expected) - set(result))}")
    for key, exp in expected.items():
        act = np.asarray(result[key])
        assert act.shape == exp.shape, f"{key}: shape {act.shape} != {exp.shape}"
        scale = float(np.max(np.abs(exp))) if exp.size else 0.0
        np.testing.assert_allclose(act, exp, rtol=RTOL, atol=RTOL * scale,
                                   err_msg=f"{name}: '{key}' changed")


@pytest.mark.data
@pytest.mark.parametrize("name", list(CASES))
def test_case_matches_golden(name):
    """Each regression case reproduces its stored reference result."""
    missing = missing_data(name)
    if missing:
        pytest.skip(f"dataset not on this machine: {missing}")
    if not os.path.exists(golden_file(name)):
        pytest.fail(f"no reference for '{name}': run `python tests/make_golden.py {name}`")
    assert_matches_golden(name, CASES[name][0]())
