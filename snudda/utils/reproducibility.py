"""Utilities for recording reproducibility metadata alongside network builds."""

import datetime
import json
import os
import subprocess


def _git_info(path):
    """Return a dict with branch, commit, clean flag, and change summary for a git repo.

    Returns None if `path` is not inside a git repo.
    """
    def run(cmd):
        try:
            r = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=path)
            return r.stdout.strip() if r.returncode == 0 else None
        except Exception:
            return None

    toplevel = run("git rev-parse --show-toplevel")
    if not toplevel:
        return None

    branch = run("git rev-parse --abbrev-ref HEAD")
    commit = run("git rev-parse HEAD")
    short_status = run("git status --porcelain")

    info = {
        "path": toplevel,
        "branch": branch,
        "commit": commit,
        "clean": short_status == "",
    }

    if short_status:
        info["changed_files"] = short_status.splitlines()
        diff_stat = run("git diff --stat HEAD")
        if diff_stat:
            info["diff_stat"] = diff_stat

    return info


def write_network_meta(network_path, snudda_data=None):
    """Write a JSON metadata file into `network_path` recording reproducibility info.

    Records:
    - build timestamp
    - snudda_data path
    - Snudda git state (branch, commit, dirty files if any)
    - bgdata git state if snudda_data lives inside a git repo

    Returns the metadata dict.
    """
    import snudda as _snudda_pkg

    meta = {
        "build_time": datetime.datetime.now().isoformat(),
        "snudda_data": str(snudda_data) if snudda_data else None,
        "python": _python_version(),
        "snudda_version": getattr(_snudda_pkg, "__version__", None),
    }

    # Snudda source tree
    snudda_src = os.path.dirname(os.path.dirname(_snudda_pkg.__file__))
    snudda_git = _git_info(snudda_src)
    if snudda_git:
        meta["snudda_git"] = snudda_git

    # bgdata (or any repo that contains snudda_data)
    if snudda_data:
        bgdata_git = _git_info(str(snudda_data))
        if bgdata_git:
            meta["bgdata_git"] = bgdata_git
        else:
            meta["bgdata_git"] = None  # not a git repo; path already recorded above

    out_file = os.path.join(network_path, "network_build_meta.json")
    os.makedirs(network_path, exist_ok=True)
    with open(out_file, "w") as f:
        json.dump(meta, f, indent=2)

    return meta


def _python_version():
    import sys
    return sys.version


def print_network_meta(network_path):
    """Pretty-print the metadata file for a network, if it exists."""
    meta_file = os.path.join(network_path, "network_build_meta.json")
    if not os.path.exists(meta_file):
        print(f"No network_build_meta.json found in {network_path}")
        return
    with open(meta_file) as f:
        meta = json.load(f)
    print(json.dumps(meta, indent=2))
    return meta
