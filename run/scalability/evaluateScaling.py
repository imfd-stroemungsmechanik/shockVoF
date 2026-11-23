import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator

plt.rcParams.update({
    "font.size": 6,
    "lines.linewidth": 0.8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "grid.linewidth": 0.4,
})



def get_times(file_path, time_type):
    """
    Time of an OpenFOAM log file from the end of the first to the end of the
    last time step, which excludes the initialization.
    """
    pattern = re.compile(
        rf"{time_type}\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    )

    times = []

    with open(file_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            match = pattern.search(line)
            if match:
                times.append(float(match.group(1)))

    return times[-1] - times[0]


# Create plot
fig, axes = plt.subplots(
    1,
    2,
    figsize=(6.3, 2.5),
    layout="constrained",
    sharey=True
)

# Array of all node configurations
nodes = np.array([1, 2, 4, 8, 16, 32, 64])

# Time steps covered by the timing: the interval between the first and the
# last time step of the 500-step runs (see get_times)
n_steps = 499



#
# Strong scaling
#

# Lists for diagrams
execution_times_strong_scaling = []

# Loop over all strong scaling cases
for node in nodes:

    # Get slurm files in case folders
    folder = Path("strongScaling/run_nodes_" + f"{node:02d}")
    files = sorted(folder.glob("slurm-*.out"))

    # List of execution times
    execution_times_per_case = []

    # Loop over all slurm files
    for file in files:

        # Execution time
        execution_times_per_case.append(get_times(file, "ClockTime"))

    # Median execution time
    execution_times_strong_scaling.append(np.median(np.array(execution_times_per_case)))

execution_times_strong_scaling = np.array(execution_times_strong_scaling)

# Achieved speedup
speedup = execution_times_strong_scaling[0] / execution_times_strong_scaling

# Ideal speedup
ideal_speedup = nodes

# Parallel efficiency
efficiency_strong = speedup / ideal_speedup * 100

# Auxiliary calculation
total_cells = 256**3
cores_per_node = 104
cores = nodes * cores_per_node
cells_per_core = total_cells / cores

# Performance index: core time per cell and time step [s]
PID_strong = execution_times_strong_scaling * cores / (total_cells * n_steps)

# Plot data
axes[0].plot(nodes, efficiency_strong, marker="o")
axes[0].axhline(y=100, color="black", linewidth=0.75, linestyle=(0, (5, 5)))
axes[0].grid(alpha=0.5)
axes[0].set_ylabel("Parallel efficiency [%]")
axes[0].set_ylim(34, 126)
axes[0].set_xscale("log")
axes[0].set_xticks(nodes, labels=[f"{n}\n({c})" for n, c in zip(nodes, cores)])
axes[0].xaxis.set_minor_locator(NullLocator())
axes[0].set_xlabel("Nodes\n(CPU cores)")

# Top axis on the strong-scaling panel
top = axes[0].secondary_xaxis("top")
top.set_ticks(nodes, labels=[f"{c/1e3:.1f}k" for c in cells_per_core])
top.xaxis.set_minor_locator(NullLocator())
top.set_xlabel("Cells per core [-]")


#
# Weak scaling
#

# Lists for diagrams
execution_times_weak_scaling = []

# Loop over all weak scaling cases
for node in nodes:

    # Get slurm files in case folders
    folder = Path("weakScaling/run_nodes_" + f"{node:02d}")
    files = sorted(folder.glob("slurm-*.out"))

    # List of execution times
    execution_times_per_case = []

    # Loop over all slurm files
    for file in files:

        # Execution time
        execution_times_per_case.append(get_times(file, "ClockTime"))

    # Median execution time
    execution_times_weak_scaling.append(np.median(np.array(execution_times_per_case)))

execution_times_weak_scaling = np.array(execution_times_weak_scaling)

# Ideal runtime
ideal_runtime = np.full(len(nodes), execution_times_weak_scaling[0])

# Parallel efficiency
efficiency_weak = ideal_runtime / execution_times_weak_scaling * 100

# Plot data
axes[1].plot(nodes, efficiency_weak, marker="o")
axes[1].axhline(y=100, color="black", linewidth=0.75, linestyle=(0, (5, 5)))
axes[1].grid(alpha=0.5)
axes[1].set_xlabel("Nodes\n(CPU cores)")
axes[1].set_xscale("log")

# Node counts with the core counts as tick labels
axes[1].set_xticks(nodes, labels=[f"{n}\n({c})" for n, c in zip(nodes, cores)])
axes[1].xaxis.set_minor_locator(NullLocator())

# Auxiliary calculation
cells_per_node = 128**3

# Performance index: core time per cell and time step [s]
PID_weak = execution_times_weak_scaling * cores / (nodes * cells_per_node * n_steps)

# Top axis on the weak-scaling panel
top = axes[1].secondary_xaxis("top")
top.set_ticks(nodes, labels=[f"{node*cells_per_node/1e6:.1f}M" for node in nodes])
top.xaxis.set_minor_locator(NullLocator())
top.set_xlabel("Total cell count [-]")

fig.savefig("results_scalability.png", dpi=450)


# Execution times, efficiencies and performance indices
print(f"{'Nodes':>5} | {'T strong [s]':>12} {'E strong [%]':>12} | {'T weak [s]':>10} {'E weak [%]':>10} | {'PID strong [us]':>15} {'PID weak [us]':>13}")

for row in zip(nodes, execution_times_strong_scaling, efficiency_strong, execution_times_weak_scaling, efficiency_weak, PID_strong, PID_weak):

    print(f"{row[0]:5d} | {row[1]:12.1f} {row[2]:12.1f} | {row[3]:10.1f} {row[4]:10.1f} | {row[5]*1e6:15.2f} {row[6]*1e6:13.2f}")
