import lzma
import json
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

PARAM_NAMES = ["U", "tauR", "tauF", "tauRatio", "cond"]

def plot_opt_params(state_file, output_file=None):
    print(f"Loading {state_file}")
    with lzma.open(state_file, "rt") as f:
        state = json.load(f)

    xi = np.array(state["xi"])   # shape: (n_evals, n_params)
    yi = np.array(state["yi"])   # shape: (n_evals,)
    n, n_params = xi.shape

    param_names = PARAM_NAMES[:n_params]
    iterations = np.arange(1, n + 1)
    min_so_far = np.minimum.accumulate(yi)

    # Colour-code points by whether they beat the running minimum
    is_new_best = np.zeros(n, dtype=bool)
    current_best = np.inf
    for i, y in enumerate(yi):
        if y < current_best:
            is_new_best[i] = True
            current_best = y

    print(f"Loaded {n} evaluations, {n_params} parameters")
    print(f"Best error: {np.min(yi):.6f} at iteration {np.argmin(yi) + 1}")
    print(f"Best params: {dict(zip(param_names, xi[np.argmin(yi)]))}")

    fig, axes = plt.subplots(n_params, 1, figsize=(12, 3 * n_params), sharex=True)
    if n_params == 1:
        axes = [axes]

    for ax, name, col in zip(axes, param_names, range(n_params)):
        values = xi[:, col]

        # All evaluations (faint)
        ax.scatter(iterations[~is_new_best], values[~is_new_best],
                   s=6, alpha=0.3, color="steelblue")

        # New-best evaluations (highlighted)
        ax.scatter(iterations[is_new_best], values[is_new_best],
                   s=30, alpha=0.9, color="crimson", zorder=3, label="New best")

        # Running best value as a step line
        best_value_so_far = np.where(is_new_best, values, np.nan)
        # Forward-fill
        last = values[0]
        running_best_param = np.empty(n)
        for i in range(n):
            if is_new_best[i]:
                last = values[i]
            running_best_param[i] = last

        ax.step(iterations, running_best_param, where="post",
                color="crimson", linewidth=1.5, alpha=0.7)

        ax.set_ylabel(name)
        ax.grid(True, alpha=0.3)
        if col == 0:
            ax.legend(loc="upper right", fontsize=8)

    axes[-1].set_xlabel("Iteration")
    fig.suptitle(f"Parameter evolution — {os.path.basename(state_file)}", fontsize=12)
    plt.tight_layout()

    if output_file is None:
        output_file = state_file.replace(".json.xz", "-params.png")
        if output_file == state_file:
            output_file = state_file + "-params.png"

    plt.savefig(output_file, dpi=150)
    print(f"Saved plot to {output_file}")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot parameter evolution from saved optimisation state.")
    parser.add_argument("state_file", type=str, help="Path to the .json.xz opt state file")
    parser.add_argument("--output", type=str, default=None, help="Output image path (optional)")
    args = parser.parse_args()
    plot_opt_params(args.state_file, args.output)
