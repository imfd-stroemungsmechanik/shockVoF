#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

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

# Case parameters
d0     = 4.8e-3
rho_l  = 1000.0
rho_g  = 2.182
U_g    = 225.9

# Ranger-Nicholls non-dimensional time scale
eps      = rho_g/rho_l
sqrt_eps = np.sqrt(eps)
T_norm   = U_g * sqrt_eps / d0

# Time of shock impact
t0 = 4.76e-6


# Load numerical data
t, _, width, leading_edge = np.loadtxt(
    "postProcessing/dropletCharacteristics/0/dropletCharacteristics",
    unpack=True,
)

# Non-dimensional quantities

# Time adjusted for shock impact
t_shifted = t - t0
T_star    = t_shifted * T_norm

mask = (T_star > 0) & (T_star < 0.8)

# Leading edge position (non-dimensionalized)
leading_edge_position = (leading_edge - leading_edge[0]) / d0

# Centerline width (non-dimensionalized)
centerline_width = width / d0


# Experimental data

# Leading edge drift
leading_edge_drift_exp = np.loadtxt(
    "experiments/leading_edge_position.csv",
    delimiter=","
)

# Centerline width
centerline_width_exp = np.loadtxt(
    "experiments/centerline_width.csv",
    delimiter=","
)


# Plots
fig, axes = plt.subplots(
    2,
    1,
    figsize=(2.95, 4),
    layout="constrained",
    sharex=True
)

# Leading edge drift
ax = axes[0]
ax.scatter(
    leading_edge_drift_exp[:, 0],
    leading_edge_drift_exp[:, 1],
    color="black",
    s=5, marker="s",
    label="Igra & Takayama (2001)"
)
ax.plot(
    T_star[mask],
    leading_edge_position[mask],
    color="#4292C6",
    label="Simulation"
)
ax.set_ylabel("Dimensionless leading edge drift [-]")
ax.grid(alpha=0.5)
ax.legend(loc="upper left")

# Centerline width
ax = axes[1]
ax.scatter(
    centerline_width_exp[:, 0],
    centerline_width_exp[:, 1],
    color="black",
    s=5, marker="s",
)
ax.plot(
    T_star[mask],
    centerline_width[mask],
    color="#4292C6",
)
ax.set_xlabel("Dimensionless time [-]")
ax.set_ylabel("Dimensionless centerline width [-]")
ax.grid(alpha=0.5)

plt.savefig("results_Igra2001.png", dpi=450)
