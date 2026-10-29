"""
Create the reference ("golden") results used by the regression tests.

Runs every case of :mod:`regression_cases` with the code currently in the
working tree and stores the output in ``tests/golden/<case>.npz``, together
with the git commit it was produced from.

Usage
-----
From the repository root::

    python tests/make_golden.py               # all cases with data on disk
    python tests/make_golden.py peo_t1        # only some cases
    python tests/make_golden.py --force       # overwrite existing files

Existing reference files are never overwritten without ``--force``:
regenerating them silently would hide exactly the changes the tests exist
to catch. Only regenerate after deciding that a change in results is
intended, and say so in the commit message.
"""
import argparse
import datetime
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import matplotlib  # noqa: E402

matplotlib.use("Agg")

from regression_cases import CASES, missing_data  # noqa: E402

#: Folder that holds the reference ``.npz`` files.
GOLDEN_DIR = os.path.join(HERE, "golden")


def git_commit():
    """
    Return the current commit hash, marked ``-dirty`` if there are changes.

    Returns
    -------
    str
        Short hash, or ``"unknown"`` if git is not available.
    """
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE,
                             capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--", "*.py"], cwd=HERE,
                               capture_output=True, text=True).stdout.strip()
        return sha + ("-dirty" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def golden_file(name):
    """
    Path of the reference file of a case.

    Parameters
    ----------
    name : str
        Case name.

    Returns
    -------
    str
        ``tests/golden/<name>.npz``.
    """
    return os.path.join(GOLDEN_DIR, f"{name}.npz")


def main(argv=None):
    """
    Command-line entry point. See the module docstring for usage.

    Parameters
    ----------
    argv : list of str, optional
        Arguments; defaults to ``sys.argv[1:]``.

    Returns
    -------
    int
        Process exit code (0 on success, 1 if any case failed).
    """
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("cases", nargs="*", help="case names (default: all)")
    parser.add_argument("--force", action="store_true", help="overwrite existing files")
    args = parser.parse_args(argv)

    names = args.cases or list(CASES)
    unknown = [n for n in names if n not in CASES]
    if unknown:
        parser.error(f"unknown case(s): {unknown}. Available: {list(CASES)}")

    os.makedirs(GOLDEN_DIR, exist_ok=True)
    commit = git_commit()
    failed = False
    for name in names:
        target = golden_file(name)
        if os.path.exists(target) and not args.force:
            print(f"[skip]    {name}: reference exists (use --force to overwrite)")
            continue
        missing = missing_data(name)
        if missing:
            print(f"[no data] {name}: missing {missing}")
            continue
        func = CASES[name][0]
        t0 = time.time()
        try:
            result = func()
        except Exception as exc:  # report and continue with the other cases
            print(f"[ERROR]   {name}: {type(exc).__name__}: {exc}")
            failed = True
            continue
        meta = {
            "_meta.commit": np.array(commit),
            "_meta.created": np.array(datetime.datetime.now().isoformat(timespec="seconds")),
        }
        np.savez_compressed(target, **result, **meta)
        size_kb = os.path.getsize(target) / 1024
        print(f"[ok]      {name}: {len(result)} arrays, {size_kb:.0f} kB, {time.time() - t0:.1f} s")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
