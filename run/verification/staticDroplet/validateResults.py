import glob
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
sigma = 8.0              # surface tension
mu    = 2.0e-3           # liquid dynamic viscosity
rho   = 1000.0           # liquid density
R     = 0.5e-3           # droplet radius
p_out = 101325.0         # outer pressure
p_in  = p_out + sigma/R  # analytical Laplace jump
t_end = 0.0025           # evaluation time
t_sigma = np.sqrt(rho*R**3/sigma)

# Mesh resolutions
resolutions = [30, 60, 120]
labels = [r"$N = 30$", r"$N = 60$", r"$N = 120$"]
colors = ["#A6CEE3", "#4292C6", "#08306B"]


# Exact pressure step
x = np.linspace(-1.5e-3, 1.5e-3, 1000)
p_exact = np.where(np.abs(x) < R, p_in, p_out) - p_out

# Create plot
fig, axes = plt.subplots(2, 1, figsize=(2.95, 4), layout="constrained")

# Reference solution for pressure
axes[0].plot(x*1e3, p_exact/1000, label="Exact", color="black", linewidth=0.75, linestyle=(0, (5, 5)))

# Read and plot data
for index, resolution in enumerate(resolutions):

    case_name = "run_lowMach_on_smoothed_curvature_n" + str(resolution)

    x, p = np.loadtxt(case_name + "/postProcessing/sample/" + str(t_end) + "/centerline.xy", unpack=True)
    t, Umax = np.loadtxt(case_name + "/postProcessing/maxFieldValues/0/volFieldValue.dat", unpack=True)

    Ca = mu * Umax / sigma

    axes[0].plot(x*1e3, (p - p_out)/1000, color=colors[index])
    axes[1].plot(t[::100]/t_sigma, Ca[::100], label=labels[index], color=colors[index])


# Formatting
axes[0].set_xlabel(r"$x$-coordinate [mm]")
axes[0].set_ylabel("Pressure jump [kPa]")

axes[1].set_xlabel("Dimensionless time [-]")
axes[1].set_ylabel("Maximum capillary number [-]")
axes[1].set_yscale("log")

for ax in axes:
    ax.grid(alpha=0.5)

# Show legend
axes[1].legend()

plt.savefig("results_staticDroplet.png", dpi=450, bbox_inches="tight")



# Time-averaged capillary number over the second half of each run
for case_name in glob.glob("run*"):

    t, Umax = np.loadtxt(case_name + "/postProcessing/maxFieldValues/0/volFieldValue.dat", unpack=True)

    Ca = mu * Umax / sigma

    half = t >= 0.5*t[-1]

    Ca_mean = np.trapezoid(Ca[half], t[half]) / (t[half][-1] - t[half][0])

    print(f"{case_name}: mean Ca = {Ca_mean:.3e}")
