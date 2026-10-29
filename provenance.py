"""
Run records: save, next to the results, what produced them.

An analysis script calls :func:`write_run_record` when it saves its output.
The record (``run_record.json``) holds the script's parameters, the time of
the run, the script name and the git commit of this repository, so any
result folder can be traced back to the exact code and settings, and the
analysis can be rerun.

If the script has uncommitted changes, the commit alone does not describe
the code that ran, so a copy of the script is saved as well.

Example
-------
>>> from provenance import write_run_record
>>> write_run_record(savepath, {"ppm_range": (300, 200), "model": "mono"}, __file__)
"""
import datetime
import json
import os
import shutil
import subprocess

import numpy as np

#: Name of the record file written in each result folder.
RECORD_NAME = "run_record.json"


def _git(args, cwd):
    """
    Run a git command and return its output, or ``None`` if it fails.

    Parameters
    ----------
    args : list of str
        Arguments after ``git``.
    cwd : str
        Folder to run in.

    Returns
    -------
    str or None
    """
    try:
        out = subprocess.run(["git", *args], cwd=cwd, capture_output=True,
                             text=True, check=True)
        return out.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _jsonable(value):
    """
    Convert values JSON cannot store (numpy, tuples, paths, dates) to plain types.

    Parameters
    ----------
    value : object
        Any parameter value.

    Returns
    -------
    object
        A JSON-serialisable equivalent; unknown objects become their ``str``.
    """
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, (datetime.date, datetime.datetime)):
        return value.isoformat()
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    return str(value)


def write_run_record(folder, parameters, script):
    """
    Write ``run_record.json`` (and, if needed, a copy of the script) to a folder.

    Parameters
    ----------
    folder : str
        Result folder; created if it does not exist.
    parameters : dict
        Everything that defines the run (paths, experiments, ranges, models).
        Values are converted to JSON-friendly types.
    script : str
        Path of the analysis script, normally ``__file__``.

    Returns
    -------
    dict
        The record that was written.

    Notes
    -----
    The ``git`` block holds ``commit`` (short hash) and ``dirty`` (whether
    the script file had uncommitted changes). When ``dirty`` is true, or
    git is unavailable, the script is copied to ``<script name>.snapshot``
    in ``folder`` and its name is stored under ``script_snapshot``.
    """
    os.makedirs(folder, exist_ok=True)
    script = os.path.abspath(script)
    repo = os.path.dirname(script)
    commit = _git(["rev-parse", "--short", "HEAD"], repo)
    status = _git(["status", "--porcelain", "--", os.path.basename(script)], repo)
    dirty = commit is None or bool(status)

    record = {
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
        "script": os.path.basename(script),
        "git": {"commit": commit, "dirty": dirty},
        "parameters": _jsonable(parameters),
    }
    if dirty:
        snapshot = os.path.basename(script) + ".snapshot"
        shutil.copyfile(script, os.path.join(folder, snapshot))
        record["script_snapshot"] = snapshot
    with open(os.path.join(folder, RECORD_NAME), "w", encoding="utf-8") as f:
        json.dump(record, f, indent=2, ensure_ascii=False)
    return record
