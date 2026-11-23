import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

plt.rcParams.update({
    "font.size": 6,
    "lines.linewidth": 0.8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "grid.linewidth": 0.4,
})

# Mach numbers
mach_numbers = [0.04, 0.02, 0.01, 0.005]

# Colors
colors = ["#A6CEE3", "#4292C6", "#0072B2", "#08306B", "#000000"]


# Case parameters
R     = 8.333333e-5      # unperturbed droplet radius R = 1/12 mm [m]
sigma = 0.072            # surface tension [N/m]
dz    = 2.0e-6           # pseudo-2D slab thickness [m]

# Semi-axes of the initial ellipse
A = 1.25*R
B = 0.8*R

# Perimeters of the ellipse (Ramanujan's approximation) and the circle
perimeter_ellipse = np.pi*(3.0*(A + B) - np.sqrt((3.0*A + B)*(A + 3.0*B)))
perimeter_circle = 2.0*np.pi*R
# Stored surface energy: the excess interface length times the surface
# tension coefficient and the slab thickness
energy_reference = sigma*(perimeter_ellipse - perimeter_circle)*dz

# Phase densities
rho_l = 1000.0
rho_g = 10.0

# Analytical period
period_theory = 2.0*np.pi*np.sqrt((rho_l + rho_g)*R**3/(6*sigma))


def read_case(case_name):
    """Kinetic-energy history, the loss per half-period [%] and the period [s]."""

    time, energy = np.loadtxt(
        case_name + "/postProcessing/kineticEnergy/0/kineticEnergy.dat",
        unpack=True
    )

    # One peak per half-period
    peaks, _ = find_peaks(
        energy,
        prominence=0.1e-12
    )

    # Exponential fit through the peaks
    n = np.arange(len(peaks))
    slope, intercept = np.polyfit(n, np.log(energy[peaks]), 1)

    # Linear fit through the peak times: t_n = T/4 + n*T/2
    half_period, t0 = np.polyfit(n, time[peaks], 1)

    return time, energy, 100*(1 - np.exp(slope)), 2.0*half_period


fig, (ax1, ax2) = plt.subplots(
    1,
    2,
    figsize=(6.3, 1.75),
    layout="constrained"
)

# Lists for losses per half-period, with and without the low-Mach correction
losses_lowMach = []
losses_hllc_plain = []

# Lists for the measured periods
periods_lowMach = []
periods_hllc_plain = []

for index, mach_number in enumerate(mach_numbers):

    case_name = "run_Ma" + str(mach_number) + "_on"

    # With the low-Mach correction: kinetic energy and decay
    time, energy, loss_on, period_on = read_case(case_name)
    losses_lowMach.append(loss_on)
    periods_lowMach.append(period_on)

    # Without it, for the decay comparison only
    case_name = "run_Ma" + str(mach_number) + "_off"
    time_off, energy_off, loss_off, period_off = read_case(case_name)
    losses_hllc_plain.append(loss_off)
    periods_hllc_plain.append(period_off)

    # Kinetic energy at the lowest Mach number (the last loop iteration)
    if index == len(mach_numbers) - 1:

        ax1.plot(
            time*1e6,
            energy/energy_reference*100,
            label="HLLC low-Mach",
            color=colors[index]
        )

# Plain HLLC reference at the lowest Mach number (the last loop iteration)
ax1.plot(
    time_off*1e6,
    energy_off/energy_reference*100,
    label="HLLC",
    color=colors[0],
    linestyle="-"
)

# Reference energy
ax1.axhline(
    100,
    linestyle="--",
    color="black"
)

ax1.set_ylim(-5, 145)
ax1.set_xlabel(r"Time [µs]")
ax1.set_ylabel("Relative kinetic energy [%]")
ax1.legend()
ax1.grid(alpha=0.5)


# Loss over Mach number, with and without the low-Mach correction
ax2.plot(
    mach_numbers,
    losses_hllc_plain,
    "s--",
    markersize=4,
    color=colors[0],
    label="HLLC"
)
ax2.plot(
    mach_numbers,
    losses_lowMach,
    "o-",
    markersize=4,
    color="#08306B",
    label="HLLC low-Mach"
)
ax2.set_xscale("log")
ax2.set_yscale("log")

# Ticks only at the Mach numbers that have data
ax2.set_xticks(
    [0.005, 0.01, 0.02, 0.04],
    labels=[r"$5\times10^{-3}$", r"$10^{-2}$", r"$2\times10^{-2}$", r"$4\times10^{-2}$"]
)
ax2.xaxis.set_minor_locator(plt.NullLocator())

ax2.set_xlabel(r"Mach number Ma")
ax2.set_ylabel("Decay per half-period [%]")
ax2.legend()
ax2.grid(alpha=0.5, which="both")

plt.savefig("results_oscillatingDroplet.png", dpi=450)


# Evaluate the period length
print(f"\nAnalytical mode-2 period: {period_theory*1e3:.3f} ms\n")
print(f"{'Ma':>6} | {'T [ms]':>10} {'err [%]':>8}")

for Ma, T in zip(mach_numbers, periods_lowMach):

    err_on = 100*(T/period_theory - 1)

    print(f"{Ma:6.3f} | {T*1e3:10.4f} {err_on:8.2f}")
