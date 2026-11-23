# Isentropic Vortex Advection

An isentropic vortex superimposed on a uniform background flow on a doubly-periodic square (Shu 1999, Spiegel et al. 2015). After one period the exact solution returns to the initial condition, so the initial data are the analytical reference. Since the flow is smooth, the case can be used to measure the formal spatial order of the reconstruction scheme. Two identical ideal-gas phases are used and the phase fraction is a smooth Gaussian that is passively advected.

## Setup

- Two-dimensional doubly-periodic square with $L$ = 0.1 m resolved by a structured mesh with N x N cells.
- The initial vortex has a strength of $\beta$ = 100 m/s and core radius $R$ = 0.01 m on a uniform background flow. End time is exactly one period.
- HLLC flux, two-stage SSP Runge-Kutta, unlimited `linearUpwind` reconstruction, inviscid, no interface sharpening.
- `./Allrun` runs the grid-convergence sweep N = 32, 64, 128, 256, 512; `validateResults.py` plots the $L_1$-error convergence of density, velocity and phase fraction.

| Parameter                | Symbol         | Value            |
|--------------------------|----------------|------------------|
| Background flow velocity | $u$            | 100 m/s          |
| Vortex strength          | $\beta$        | 100 m/s          |
| Core radius              | $R$            | 0.01 m           |
| Phase fraction           | $\alpha$       | Gaussian profile |
| End time                 | $t_\text{end}$ | 0.001 s          |

## References

1. Shu, C.-W. (1999). High order ENO and WENO schemes for computational fluid dynamics. In *High-Order Methods for Computational Physics*, Lecture Notes in Computational Science and Engineering, 9, Springer, 439-582.
2. Spiegel, S. C., Huynh, H. T. & DeBonis, J. R. (2015). A survey of the isentropic Euler vortex problem using high-order methods. *AIAA Paper* 2015-2444.
