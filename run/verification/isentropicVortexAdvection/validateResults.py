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

# List of mesh resolutions
resolutions = [32, 64, 128, 256, 512]

# Dictionary for the results
results = {
    "L1_rho": [],
    "L1_U": [],
    "L1_alpha": []
}

# Read the L1 norms
for resolution in resolutions:

    # Case folder name
    case_name = "run_n" + str(resolution)

    # Read raw data
    time, L1_rho, L1_U, L1_alpha = np.loadtxt(case_name + "/postProcessing/convergence/0/errors.dat", unpack=True)

    # Sort data into dictionary
    results["L1_rho"].append(float(L1_rho))
    results["L1_U"].append(float(L1_U))
    results["L1_alpha"].append(float(L1_alpha))

# Field names and their y-labels
fields = [("L1_rho", r"$L_1$ error in $\rho$ [kg/m$^3$]"),
          ("L1_U", r"$L_1$ error in $|U|$ [m/s]"),
          ("L1_alpha", r"$L_1$ error in $\alpha$ [-]")]

# Cell size
dx = 0.1 / np.array(resolutions)

# Observed order of convergence
orders = {}

for key, _ in fields:

    L1 = np.array(results[key])

    # Average order: least-squares slope of ln(error) vs ln(dx)
    p = np.polyfit(np.log(dx), np.log(L1), 1)[0]

    orders[key] = p

    print(f"{key:9s} p = {p:.2f}")

# Create plot
fig, axes = plt.subplots(1, 3, figsize=(6.3, 1.75), layout="constrained")

# Plot data per field
for ax, (key, ylabel) in zip(axes, fields):

    # Simulation results
    L1 = np.array(results[key])
    ax.plot(np.log2(resolutions), L1, "o-", color="#4292C6", markersize=4)

    # Second-order reference line (anchored to the finest solution)
    ref2 = L1[-1] * (dx / dx[-1])**2
    ax.plot(np.log2(resolutions), ref2, color="black", linewidth=0.75, linestyle=(0, (5, 5)))

    # Order of accuracy
    ax.text(0.05, 0.08, rf"$p = {orders[key]:.2f}$", transform=ax.transAxes)

    # Formatting
    ax.set_yscale("log")
    ax.set_xticks(np.log2(resolutions))
    ax.set_xticklabels(resolutions)
    ax.set_xlabel("Cells per direction [-]")
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.5)

plt.savefig("results_isentropicVortexAdvection.png", dpi=450, bbox_inches="tight")
