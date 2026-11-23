import sys
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

# Case parameters (must match the OpenFOAM setup)
p0, u0 = 101325.0, 100.0
t_end = 0.01

# Read OpenFOAM data
x, p_of, U_of, alpha_of = np.loadtxt("postProcessing/sample/" + str(t_end) + "/data.xy", unpack=True)

# Exact solution: uniform p, T, u; the alpha step returns to its start
x_ex = np.linspace(0, 1, 2001)
alpha_ex = np.where((x_ex >= 0.25) & (x_ex <= 0.75), 1.0, 0.0)

# Oscillation metrics
dp = p_of - p0
dU = U_of - u0

# Plot
fig, axes = plt.subplots(1, 3, figsize=(6.3, 1.75), layout="constrained")

axes[0].axhline(0.0, color="black", linewidth=0.75, linestyle=(0, (5, 5)), label="Exact")
axes[0].plot(x, dp*1e5, color="#4292C6", label="Simulation")
axes[0].set_ylabel(r"Pressure error / $10^{-5}$ [Pa]")
axes[0].legend(loc="lower left")

axes[1].axhline(0.0, color="black", linewidth=0.75, linestyle=(0, (5, 5)))
axes[1].plot(x, dU*1e9, color="#4292C6")
axes[1].set_ylabel(r"Velocity error / $10^{-9}$ [m/s]")

axes[2].plot(x_ex, alpha_ex, color="black", linewidth=0.75, linestyle=(0, (5, 5)))
axes[2].plot(x, alpha_of, color="#4292C6")
axes[2].set_ylabel("Phase fraction [-]")

for ax in axes:

    ax.set_xlabel(r"$x$-coordinate [m]")
    ax.grid(alpha=0.5)

plt.savefig("results_interfaceAdvection.png", dpi=450, bbox_inches="tight")

# Test: p and U must stay uniform to about nine significant digits
tolerance = 1e-9

errors = {
    "p": np.max(np.abs(dp))/p0,
    "U": np.max(np.abs(dU))/u0,
}

for name, error in errors.items():
    print(f"Maximum relative {name} error: {error:.3e} (tolerance {tolerance:.0e})")

sys.exit(0 if all(error <= tolerance for error in errors.values()) else 1)
