import lzma
import json
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os


def plot_opt_state(state_file, output_file=None):

    print(f"Loading {state_file}")
    with lzma.open(state_file, "rt") as f:
        state = json.load(f)

    yi = np.array(state["yi"])
    n = len(yi)
    iterations = np.arange(1, n + 1)
    min_so_far = np.minimum.accumulate(yi)

    print(f"Loaded {n} evaluations")
    print(f"Best error: {np.min(yi):.6f} at iteration {np.argmin(yi) + 1}")

    fig, ax = plt.subplots(figsize=(10, 5))

    ax.scatter(iterations, yi, s=6, alpha=0.4, color="steelblue", label="Error per evaluation")
    ax.plot(iterations, min_so_far, color="crimson", linewidth=2, label="Min error so far")

    ax.set_xlabel("Iteration")
    ax.set_ylabel("Error")
    ax.set_title(os.path.basename(state_file))
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_file is None:
        output_file = state_file.replace(".json.xz", "-plot.png")
        if output_file == state_file:
            output_file = state_file + "-error.png"

    plt.savefig(output_file, dpi=150)
    print(f"Saved plot to {output_file}")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot optimisation error from saved state file.")
    parser.add_argument("state_file", type=str, help="Path to the .json.xz opt state file")
    parser.add_argument("--output", type=str, default=None, help="Output image path (optional)")
    args = parser.parse_args()

    plot_opt_state(args.state_file, args.output)
