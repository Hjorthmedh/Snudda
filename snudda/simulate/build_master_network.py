"""Build the master striatal network used by run_pair_recording_from_master.py.

The master is a 1 mm × 1 mm × 400 µm slab with all 7 cell types at realistic
densities (~85 000 neurons/mm³, fractions from snudda_data).  It is built once
and then each experiment sub-slices it in-process (no place/detect/prune per run).

Usage (local test):
    python build_master_network.py <network_path> [--snudda_data PATH]

Usage (Dardel, inside SLURM):
    srun -n <N> python build_master_network.py <network_path>

The network_path will contain network-synapses.hdf5 when done.
"""

import argparse
import os

from snudda.core import Snudda


# Realistic cell-type fractions and density for mouse striatum (snudda_data).
# Volume: 1 mm × 1 mm × 400 µm = 0.4 mm³, density ~85 000 /mm³.
MASTER_VOLUME_SIDE = 1e-3     # metres
MASTER_VOLUME_DEPTH = 400e-6  # metres
NEURON_DENSITY = 85000        # per mm³

CELL_FRACTIONS = {
    "dSPN": 0.475,
    "iSPN": 0.475,
    "FS":   0.013,
    "ChIN": 0.011,
    "LTS":  0.007,
    "TH":   0.005,
    "NGF":  0.0019,
}

VOLUME_MM3 = (MASTER_VOLUME_SIDE * 1e3) ** 2 * (MASTER_VOLUME_DEPTH * 1e3)

NEURON_COUNTS = {k: max(1, round(v * NEURON_DENSITY * VOLUME_MM3))
                 for k, v in CELL_FRACTIONS.items()}


def build(network_path, snudda_data=None, random_seed=None, parallel=True):
    if snudda_data is None:
        snudda_data = os.environ.get("SNUDDA_DATA")
    if snudda_data is None:
        raise ValueError("snudda_data not supplied and SNUDDA_DATA env var is not set")

    neuron_names = list(NEURON_COUNTS.keys())
    number_of_neurons = [NEURON_COUNTS[n] for n in neuron_names]
    neuron_paths = [os.path.join(snudda_data, "neurons", "striatum", n.lower())
                    for n in neuron_names]
    connection_config = os.path.join(snudda_data, "connectivity", "striatum",
                                     "striatum-connectivity.json")

    print(f"Building master network in: {network_path}")
    print(f"Volume: {MASTER_VOLUME_SIDE*1e3:.0f} mm × "
          f"{MASTER_VOLUME_SIDE*1e3:.0f} mm × "
          f"{MASTER_VOLUME_DEPTH*1e6:.0f} µm")
    for n, c in zip(neuron_names, number_of_neurons):
        print(f"  {n}: {c}")
    print(f"  Total: {sum(number_of_neurons)}")

    sn = Snudda(network_path=network_path, snudda_data=snudda_data,
                random_seed=random_seed)

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

    sn.create_network(parallel=parallel)

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
    args = parser.parse_args()

    build(network_path=args.network_path,
          snudda_data=args.snudda_data,
          random_seed=args.random_seed,
          parallel=not args.no_parallel)
