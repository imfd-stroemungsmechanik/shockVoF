#!/usr/bin/env python3
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

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

# Colors for the plots
colors = ["#A6CEE3", "#4292C6", "#0072B2", "#08306B", "#000000"]

# Case parameters
d_0 = 0.00076

# The EPL peak is located by a parabola fitted to the points above this
# fraction of the maximum
peak_fraction = 0.8

# Sampling time: the latest time written by the trajectory function object
time = max(os.listdir("postProcessing/trajectory"), key=float)

# Heights y above the wall [mm] of the EPL profiles. The experimental files
# are named after the vertical coordinate z of Lin et al. (2020), which is y
# here
PROFILE_HEIGHTS = [0.5, 1.5, 2.5, 3.0, 4.0]

# Interpolation grid of the horizontal (x-z) planes
x_bounds_xz = (-1.0e-3, 6.0e-3)
z_bounds = (-3.0e-3, 3.0e-3)
nx_xz, nz = 600, 600

x_grid_xz = np.linspace(*x_bounds_xz, nx_xz)
z_grid = np.linspace(*z_bounds, nz)
mesh_grid_xz_x, mesh_grid_xz_z = np.meshgrid(x_grid_xz, z_grid, indexing="ij")


def epl_profile(file_name):
    """Equivalent path length (EPL) along x_grid_xz on a horizontal plane:
    the liquid volume fraction integrated along the z-axis"""

    x, z, alpha = np.loadtxt(file_name, usecols=(0, 2, 3), unpack=True)

    # Only the points around the interpolation grid, which speeds up the
    # triangulation of planes through the refined wall region
    margin = 0.5e-3
    near = (
        (x > x_bounds_xz[0] - margin) & (x < x_bounds_xz[1] + margin)
        & (z > z_bounds[0] - margin) & (z < z_bounds[1] + margin)
    )

    alpha_g = griddata(
        (x[near], z[near]),
        alpha[near],
        (mesh_grid_xz_x, mesh_grid_xz_z),
        method="linear",
        fill_value=0.0
    )

    return np.trapezoid(alpha_g, z_grid, axis=1)


def peak_location(x, f):
    """Location of the maximum of f(x): vertex of a parabola fitted to the
    contiguous points around it above peak_fraction of the maximum"""

    i = np.argmax(f)

    below = np.flatnonzero(f < peak_fraction * f[i])
    lo = below[below < i].max(initial=-1) + 1
    hi = below[below > i].min(initial=len(f)) - 1

    # No peak within the sampled range, e.g., no liquid on the plane
    if lo == 0 or hi == len(f) - 1:
        return np.nan

    a, b, _ = np.polyfit(x[lo:hi + 1] - x[i], f[lo:hi + 1], 2)

    return np.clip(x[i] - b / (2 * a), x[lo], x[hi]) if a < 0 else x[i]


# EPL profiles on all horizontal planes, which are named y_<height in mm>
epl = {
    float(os.path.basename(file_name)[2:-3]): epl_profile(file_name)
    for file_name in glob.glob(f"postProcessing/trajectory/{time}/y_*.xy")
}

heights = np.array(sorted(epl))

# Create plot
fig, axes = plt.subplots(
    1,
    2,
    figsize=(6.3, 2.25),
    layout="constrained"
)



# Jet trajectory

# Experimental reference data (non-dimensionalized)
EXPERIMENT = np.array([
    [0.027041797388928247, 0.1351282757428977],
    [0.14055257565023785, 0.6600031319693067],
    [0.12271362547650373, 1.3238942258365867],
    [0.30787707047194957, 1.9864893324045423],
    [0.6151873966872449, 2.6489548402425656],
    [0.8490638599332052, 2.9720876735407993],
    [1.1049623857234423, 3.309476367131602],
    [1.756475390480046, 3.9516380739468753],
    [2.388638127045171, 4.627236253084719],
    [3.5738083189960532, 5.278599469725196],
    [4.511446209897784, 5.934757839373174]
])

# Trajectory: location of the maximum EPL at each height
x_t = np.array([peak_location(x_grid_xz, epl[height]) for height in heights])


axes[0].scatter(
    EXPERIMENT[:, 0], EXPERIMENT[:, 1],
    color="black",
    s=5, marker="s",
    label="Lin et al. (2020)",
)

axes[0].plot(
    x_t / d_0,
    heights / (d_0 * 1e3),
    label="Simulation",
    color="#4292C6"
)

axes[0].grid(alpha=0.5)

axes[0].set_xlabel(r"Dimensionless $x$-coordinate [-]")
axes[0].set_ylabel(r"Dimensionless $y$-coordinate [-]")
axes[0].legend(loc="upper left")



# EPL profiles

for i, y_mm in enumerate(PROFILE_HEIGHTS):

    x_exp, epl_exp = np.loadtxt(
        "experiments/z_" + str(y_mm) + ".csv",
        delimiter=",",
        unpack=True
    )

    axes[1].scatter(
        x_exp, epl_exp,
        color=colors[i],
        s=5, marker="s",
    )

    axes[1].plot(
        x_grid_xz * 1e3, epl[y_mm] * 1e3,
        color=colors[i],
        linewidth=0.75,
        label=f"y = {y_mm} mm",
    )

axes[1].set_xlabel(r"$x$-coordinate [mm]")
axes[1].set_ylabel("EPL [mm]")
axes[1].grid(alpha=0.5)
axes[1].legend(loc="upper right")

fig.savefig("results_Lin2020.png", dpi=450)
