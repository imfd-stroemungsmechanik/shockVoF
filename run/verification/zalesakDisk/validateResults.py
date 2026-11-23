import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

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

labels = ["No sharpening", r"$\rho$-THINC"]
modes = ["none", "THINC"]
resolutions = np.array([32, 64, 128, 256])



# alpha = 0.5 isocontour of the disk, all resolutions, with / without sharpening

# Slotted-disk geometry
diskX, diskY = 0.5, 0.75
radius = 0.15
slotHalfWidth = 0.025
slotTop = 0.85

# Exact solution
GX, GY = np.meshgrid(np.linspace(0.3, 0.7, 400), np.linspace(0.55, 0.95, 400))
inDisk = (GX - diskX)**2 + (GY - diskY)**2 <= radius**2
inSlot = (np.abs(GX - diskX) <= slotHalfWidth) & (GY <= slotTop)
alpha_ex = np.where(inDisk & ~inSlot, 1.0, 0.0)

fig, axes = plt.subplots(
    2,
    len(resolutions),
    figsize=(6.3, 2.75),
    layout="constrained",
    sharex=True,
    sharey=True
)

# Loop over all modes (rows) and mesh resolutions (columns)
for row, mode in enumerate(modes):

    for col, resolution in enumerate(resolutions):

        ax = axes[row, col]

        # Assemble case name
        case_name = "run_" + mode + "_" + str(resolution)

        # Exact interface
        ax.contour(GX, GY, alpha_ex, levels=[0.5], colors="black", linewidths=0.75,
                   linestyles=[(0, (5, 5))])

        # Simulated interface from the sampled cutting plane
        x, y, alpha = np.loadtxt(case_name + "/postProcessing/sampleSurface/0.125/disk.xy", usecols=(0, 1, 3), unpack=True)
        ax.tricontour(x, y, alpha, levels=[0.5], colors="#4292C6", linewidths=1.25)

        ax.set_aspect("equal")
        ax.set_xlim(0.25, 0.75)
        ax.set_ylim(0.55, 0.95)
        ax.grid(alpha=0.5)
        ax.xaxis.set_major_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(0.1))

        # Mesh resolution as column header
        if row == 0:
            ax.set_title(str(resolution) + r"$^2$", size=7)

    # Sharpening mode as row label
    axes[row, 0].set_ylabel(labels[row])

fig.supxlabel(r"$x$-coordinate [m]", size=6)
fig.supylabel(r"$y$-coordinate [m]", size=6)

handles = [Line2D([0], [0], color="black", label="Exact", linestyle=(0, (5, 5)), linewidth=0.75),
           Line2D([0], [0], color="#4292C6", label="Simulation", linewidth=1.25)]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.02),
           ncol=2, frameon=False)

plt.savefig("results_zalesakDisk_isoAlpha.png", dpi=450, bbox_inches="tight")





# alpha distribution along y = 0.7501 (through the slot), with / without sharpening

# Exact top-hat: disk minus slot at the sampling height y = 0.7501
xe = np.linspace(0.3, 0.7, 2001)
half = np.sqrt(radius**2 - (0.7501 - diskY)**2)
alpha_ex = np.where((np.abs(xe - diskX) <= half)
                    & (np.abs(xe - diskX) > slotHalfWidth), 1.0, 0.0)

# One color per mesh resolution (coarse -> fine)
res_colors = ["#A6CEE3", "#4292C6", "#0072B2", "#08306B"]

fig, axes = plt.subplots(
    1,
    2,
    figsize=(6.3, 2.4),
    sharey=True,
    layout="constrained")

# Loop over all modes (columns)
for col, mode in enumerate(modes):

    ax = axes[col]

    # Exact interface
    ax.plot(xe, alpha_ex, color="black", linewidth=0.75, linestyle=(0, (5, 5)))

    # Simulated profiles for all mesh resolutions
    for index, resolution in enumerate(resolutions):

        # Assemble case name
        case_name = "run_" + mode + "_" + str(resolution)

        # Read the sampled line (columns x alpha) and plot
        x, alpha = np.loadtxt(case_name + "/postProcessing/slotProfile/0.125/slot.xy", unpack=True)
        ax.plot(x, alpha, color=res_colors[index])

    ax.set_title(labels[col], size=7)
    ax.set_xlabel(r"$x$-coordinate [m]")
    ax.set_xlim(0.3, 0.7)
    ax.grid(alpha=0.5)

axes[0].set_ylabel("Phase fraction [-]")

# Single legend: exact reference plus the four mesh resolutions
handles = [Line2D([0], [0], color="black", label="Exact", linestyle=(0, (5, 5)), linewidth=0.75)]
handles += [Line2D([0], [0], color=res_colors[i], label=str(r) + r"$^2$")
            for i, r in enumerate(resolutions)]
fig.legend(handles=handles, loc="outside upper center", ncol=5)

plt.savefig("results_zalesakDisk_slotProfile.png", dpi=450, bbox_inches="tight")
