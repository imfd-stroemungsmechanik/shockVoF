# Static Droplet

A circular liquid droplet is placed in stationary gas initialized with the exact Laplace pressure jump across the interface. The well-balanced surface-tension scheme must maintain the pressure jump while keeping the spurious currents small. The case compares the centerline pressure against the exact distribution and tracks the spurious capillary number:
$$
\text{Ca} = \frac{\mu_\text{liq} \, |U|_\text{max}}{\sigma}
$$
over a mesh refinement.

## Setup

- Two-dimensional square domain of size 6R x 6R with non-reflecting boundary conditions. The simulation covers 20 capillary times.
- Circular droplet in the center of the domain with a diffuse interface consistent with the $\rho$-THINC interface sharpening. Pressure initialization follows the same profile with the exact curvature.
- HLLC flux, two-stage SSP Runge-Kutta, vanLeer (MUSCL) reconstruction, viscous, $\rho$-THINC interface sharpening, smoothed curvature model.
- `./Allrun` sweeps over the mesh resolution with N = 30, 60, 120 cells in each direction and runs two additional cases at N = 60 without the low-Mach correction, one with the smoothed curvature model and one with the constant exact curvature $\kappa = 1/R$; `validateResults.py` plots the centerline pressure profiles and the maximum capillary number over time for the mesh refinement sweep and prints the mean capillary number over the second half of each run.

| Parameter        | Symbol                                                       | Value     |
|------------------|--------------------------------------------------------------|-----------|
| Droplet radius   | $R$                                                          | 0.5 mm    |
| Surface tension  | $\sigma$                                                     | 8 N/m     |
| Pressure jump    | $\Delta p$                                                   | 16000 Pa  |
| Laplace number   | $\text{La} = \rho_\text{liq} \, \sigma R / \mu_\text{liq}^2$ | 1e6       |
| Ohnesorge number | $\text{Oh} = 1/\sqrt{\text{La}}$                             | 1e-3      |
| Density ratio    | $\rho_\text{liq} / \rho_\text{gas}$                          | 1000      |
| Viscosity ratio  | $\mu_\text{liq} / \mu_\text{gas}$                            | 100       |
| Capillary time   | $t_\sigma = \sqrt{ \rho_\text{liq} R^3/\sigma}$              | 1.25e-4 s |
