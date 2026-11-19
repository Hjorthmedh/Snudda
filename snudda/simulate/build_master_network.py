"""Build the master striatal network used by run_pair_recording_from_master.py.

Cell-type fractions and overall density are read at runtime from the snudda_data
neuron JSON files (neurons/striatum/<type>/<type>.json), so the master network
always reflects the canonical values in bgdata/snudda_data rather than duplicating
them here.

The master is a 500 µm × 500 µm × 400 µm slab (~8 500 neurons).  It is built once
and then each experiment sub-slices it in-process (no place/detect/prune per run).

Usage (local test):
    python build_master_network.py <network_path> [--snudda_data PATH]

Usage (Dardel, inside SLURM):
    python build_master_network.py <network_path> --n_workers $SLURM_NTASKS

Place/detect/prune use ipyparallel (not MPI), so run as plain python, not srun.
The network_path will contain network-synapses.hdf5 when done.
"""

import argparse
import json
import os

from snudda.core import Snudda
from snudda.utils.reproducibility import write_network_meta


MASTER_VOLUME_SIDE  = 500e-6   # metres
MASTER_VOLUME_DEPTH = 400e-6   # metres
NEURON_DENSITY      = 80500    # per mm³ (matches init_tiny default)

# Cell types to include (directory names under snudda_data/neurons/striatum/).
CELL_TYPES = ["dspn", "ispn", "fs", "chin", "lts", "th", "ngf"]

# Display names matching the JSON keys in each cell type's metadata file.
_DISPLAY_NAME = {
    "dspn": "dSPN", "ispn": "iSPN", "fs": "FS",
    "chin": "ChIN", "lts": "LTS", "th": "TH", "ngf": "NGF",
}

def _read_fractions(snudda_data):
    """Read cell-type fractions from neuron JSON files.

    For each cell type, looks for neurons/striatum/<type>/<type>.json in:
      1. snudda_data  (deployed copy, may not have all types)
      2. ~/BasalGangliaData/data  (source repo, always complete)

    The JSON has a top-level display-name key whose value holds 'fraction'.
    Returns {display_name: fraction}.
    """
    bgdata_root = os.path.expanduser("~/BasalGangliaData/data")
    search_roots = [snudda_data, bgdata_root]

    fractions = {}
    for ct in CELL_TYPES:
        display = _DISPLAY_NAME[ct]
        frac = None
        for root in search_roots:
            jf = os.path.join(root, "neurons", "striatum", ct, f"{ct}.json")
            if os.path.exists(jf):
                data = json.load(open(jf))
                top = next(iter(data.values())) if isinstance(next(iter(data.values())), dict) else data
                frac = top.get("fraction")
                if frac is not None:
                    fractions[display] = float(frac)
                    break
        if frac is None:
            raise ValueError(
                f"No fraction found for {display}: checked {search_roots}"
            )

    total = sum(fractions.values())
    if abs(total - 1.0) > 0.05:
        print(f"WARNING: cell-type fractions sum to {total:.4f} (expected ~1.0)")
    return fractions


def build(network_path, snudda_data=None, random_seed=None, parallel=True, n_workers=None):
    if snudda_data is None:
        snudda_data = os.environ.get("SNUDDA_DATA")
    if snudda_data is None:
        raise ValueError("snudda_data not supplied and SNUDDA_DATA env var is not set")

    fractions = _read_fractions(snudda_data)

    volume_mm3 = (MASTER_VOLUME_SIDE * 1e3) ** 2 * (MASTER_VOLUME_DEPTH * 1e3)
    neuron_counts = {name: max(1, round(frac * NEURON_DENSITY * volume_mm3))
                     for name, frac in fractions.items()}

    neuron_names = list(neuron_counts.keys())
    number_of_neurons = [neuron_counts[n] for n in neuron_names]
    _reverse = {v: k for k, v in _DISPLAY_NAME.items()}
    neuron_paths = [os.path.join(snudda_data, "neurons", "striatum", _reverse[n])
                    for n in neuron_names]
    connection_config = os.path.join(snudda_data, "connectivity", "striatum",
                                     "striatum-connectivity.json")

    print(f"Building master network in: {network_path}")
    print(f"snudda_data: {snudda_data}")
    print(f"Volume: {MASTER_VOLUME_SIDE*1e6:.0f} µm × "
          f"{MASTER_VOLUME_SIDE*1e6:.0f} µm × "
          f"{MASTER_VOLUME_DEPTH*1e6:.0f} µm  ({volume_mm3:.3f} mm³)")
    print(f"Density: {NEURON_DENSITY} /mm³")
    for n, c in zip(neuron_names, number_of_neurons):
        print(f"  {n}: {c}  (fraction {fractions[n]:.4f})")
    print(f"  Total: {sum(number_of_neurons)}")

    sn = Snudda(network_path=network_path, parallel=parallel, n_workers=n_workers)

    sn.init_tiny(neuron_paths=neuron_paths,
                 neuron_names=neuron_names,
                 number_of_neurons=number_of_neurons,
                 snudda_data=snudda_data,
                 connection_config=connection_config,
                 random_seed=random_seed,
                 density=NEURON_DENSITY,
                 volume_mesh="slice",
                 side_len=MASTER_VOLUME_SIDE,
                 slice_depth=MASTER_VOLUME_DEPTH)

    sn.create_network()

    meta = write_network_meta(network_path, snudda_data=snudda_data)
    if not meta.get("snudda_git", {}).get("clean", True):
        print("WARNING: Snudda repo has uncommitted changes — see network_build_meta.json")

    print(f"Master network complete: {network_path}/network-synapses.hdf5")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build master striatal network")
    parser.add_argument("network_path", type=str,
                        help="Output directory for master network")
    parser.add_argument("--snudda_data", type=str, default=None,
                        help="Path to snudda_data (default: $SNUDDA_DATA)")
    parser.add_argument("--random_seed", type=int, default=None)
    parser.add_argument("--no_parallel", action="store_true",
                        help="Disable ipyparallel (single-process, for testing)")
    parser.add_argument("--n_workers", type=int, default=None,
                        help="Number of ipyparallel workers (default: $N_WORKERS or auto)")
    args = parser.parse_args()

    n_workers = args.n_workers or (int(os.environ["N_WORKERS"])
                                   if "N_WORKERS" in os.environ else None)

    build(network_path=args.network_path,
          snudda_data=args.snudda_data,
          random_seed=args.random_seed,
          parallel=not args.no_parallel,
          n_workers=n_workers)
