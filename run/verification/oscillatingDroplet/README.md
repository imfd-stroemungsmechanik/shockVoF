# Oscillating Droplet

An elliptical droplet oscillates under surface tension according to Shukla (2014) and Garrick et al. (2017) with the analytical period length:
$$
T = 2 \pi \sqrt{\frac{(\rho_\text{liq} + \rho_\text{gas}) R^3}{6 \, \sigma}}~.
$$
Since the flow is inviscid, the decay of the kinetic-energy peaks is due to numerical dissipation. The case is run at four acoustic Mach numbers on a fixed mesh, with and without the low-Mach correction of the HLLC flux.

## Setup

- Two-dimensional square domain of size 4R x 4R with non-reflecting boundary conditions. The structured mesh has 120 x 120 cells in the x- and y-directions.
- Elliptical droplet initialized by
    $$
    1 = \frac{x^2}{a^2} + \frac{y^2}{b^2}
    $$
    with $a = 5/4 R$ and $b = 4/5 R$ and a smoothed interface. The pressure is the ambient pressure outside and the Laplace pressure inside.
- HLLC flux (with and without low-Mach correction), two-stage SSP Runge-Kutta, vanLeer (MUSCL) reconstruction, inviscid, $\rho$-THINC interface sharpening, smoothed curvature model.
- The speed of sound is varied by scaling the ambient pressure, the temperature and the pressure $p + p_{\infty,\text{liq}}$ inside the droplet by
    $$
    k = \left( \frac{u_\sigma}{\text{Ma} \, c_\text{ref}} \right)^2 ~,
    $$
    with the Mach number $\text{Ma} = u_\sigma / c$ and capillary velocity $u_\sigma = \sqrt{\sigma / (\rho_\text{liq} R)}$. The phase densities (1000 kg/m3 inside the droplet, 10 kg/m3 outside), the speeds of sound $\sqrt{k} \, c_\text{ref}$ of both phases, the surface tension coefficient and the equivalent droplet radius are then the same for all Mach numbers, and with them the density ratio, the analytical period and the stored surface energy.
- `./Allrun` sweeps over the Mach numbers and the flux scheme settings into the eight run directories; `validateResults.py` plots the kinetic-energy history and the decay per half-period against Ma, both with and without the correction, and prints the measured oscillation periods.

| Parameter                | Symbol                              | Value        |
|--------------------------|-------------------------------------|--------------|
| Droplet radius           | $R$                                 | 1/12 mm      |
| Surface tension          | $\sigma$                            | 0.072 N/m    |
| Capillary velocity       | $u_\sigma$                          | 0.9295 m/s   |
| Reference speed of sound | $c_\text{ref}$                      | 155.54 m/s   |
| Density ratio            | $\rho_\text{liq} / \rho_\text{gas}$ | 100          |
| Viscosity                | $\mu$                               | inviscid     |
| Mach numbers             | Ma                                  | 0.005 - 0.04 |
| Period length            | $T$                                 | 0.231 ms     |

## References

1. Garrick, D. P., Owkes, M. & Regele, J. D. (2017). A finite-volume HLLC-based scheme for compressible interfacial flows with surface tension. *Journal of Computational Physics*, 339, 46-67.
2. Shukla, R. K. (2014). Nonlinear preconditioning for efficient and accurate interface capturing in simulation of multicomponent compressible flows. *Journal of Computational Physics*, 276, 508-540.
