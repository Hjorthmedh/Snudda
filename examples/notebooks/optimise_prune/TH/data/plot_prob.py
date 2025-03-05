import numpy as np
import matplotlib.pyplot as plt

def plot_prob(filename, title):

    data = np.genfromtxt(filename)

    distances = np.linalg.norm(data[:,:3], axis=1)
    connected = data[:,3]

    bins = np.arange(0, 550, 50)
    bin_indices = np.digitize(distances, bins) - 1
    bin_counts = np.bincount(bin_indices, minlength=len(bins))
    bin_connected_counts = np.bincount(bin_indices, weights=connected, minlength=len(bins))  # Connected counts per bin

    fraction_connected = np.divide(bin_connected_counts, bin_counts, 
                                   where=bin_counts > 0, 
                                   out=np.zeros_like(bin_counts, dtype=float))
    
    
    bin_centers = (bins[:-1] + bins[1:]) / 2  # Midpoints for plotting

    plt.figure()
    plt.bar(bin_centers, fraction_connected[:-1], width=18, align='center', alpha=0.7, edgecolor='k')
    
    idx1 = connected == 1
    idx0 = connected == 0
    plt.scatter(distances[idx1], np.random.uniform(size=(np.sum(idx1),)), marker='*', color="red")
    plt.scatter(distances[idx0], np.random.uniform(size=(np.sum(idx0),)), marker='o', color="grey")
    
    
    plt.xlabel("Distance Range")
    plt.ylabel("Fraction Connected")
    plt.title("Fraction of Connected vs Distance")
    plt.xticks(bins)
    plt.ylim(0, 1)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.title(f"{title} {np.sum(bin_connected_counts)}/{np.sum(bin_counts)}")
    plt.ion()
    plt.show()

    plot_fig = f"{filename}.png"
    plt.savefig(plot_fig, dpi=150)
    print(f"Writing figure {plot_fig}")
    


plot_prob("ChIN_to_TH.csv", "ChIN to TH")
plot_prob("TH_to_ChIN.csv", "TH to ChIN")
