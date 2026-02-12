import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# -------------------------
# USER SETTINGS
# -------------------------

# This code plots the domain percolation length and net spread (of the spider web)
# Youl will need to add the directory path of results 
#
results_dir = Path("./")   # directory with binary files
iterations = 100                # must match C++ Iterations
dtype = np.float64

# -------------------------
# HELPERS
# -------------------------
def load_binary(fname):
    data = np.fromfile(fname, dtype=dtype)
    if data.size != iterations:
        raise ValueError(f"{fname.name} has wrong size: {data.size}")
    return data

def plot_with_stats(data, ylabel, title, outfile):
    mean = np.mean(data)
    std = np.std(data)

    print(f"{title}")
    print(f"  Mean = {mean:.3f}")
    print(f"  Std  = {std:.3f}\n")

    x = np.arange(len(data))

    plt.figure(figsize=(7,4))
    plt.plot(x, data, 'o-', alpha=0.7, label="Microstructures")
    plt.axhline(mean, color='red', linestyle='--', linewidth=2, label=f"Mean = {mean:.2f}")
    plt.fill_between(
        x,
        mean - std,
        mean + std,
        color='red',
        alpha=0.2,
        label=f"±1σ = {std:.2f}"
    )

    plt.xlabel("Microstructure index")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)
    plt.close()

# -------------------------
# LOAD FILES
# -------------------------
pl_files = sorted(f for f in results_dir.glob("PL_*") if not f.name.startswith("PL_tot_"))

pl_tot_files = sorted(results_dir.glob("PL_tot_*"))

if len(pl_files) == 0 or len(pl_tot_files) == 0:
    raise RuntimeError("No percolation files found")

# -------------------------
# PROCESS EACH DATASET
# -------------------------
for f in pl_files:
    data = load_binary(f)
    label = f.name.replace("PL_", "")

    plot_with_stats(
        data,
        ylabel="Domain percolation length (Mean)",
        title=f"Domain Percolation length (Longest Thread)",
        outfile=f"percolation_max_{label}.png"
    )

for f in pl_tot_files:
    data = load_binary(f)
    label = f.name.replace("PL_tot_", "")

    plot_with_stats(
        data,
        ylabel="Domain Percolation Spread (Mean)",
        title=f"Domain percolation spread (Total Spread) of percolating domains",
        outfile=f"percolation_spread_{label}.png"
    )
