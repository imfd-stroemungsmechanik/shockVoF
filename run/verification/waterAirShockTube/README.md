# Water-Air Shock Tube

Two-phase shock tube with stiffened-gas water and ideal-gas air according to Saurel et al. (2009). It tests the stiffened-gas EoS under extreme conditions and oscillation-free transport of a material interface across a large density and impedance jump.

## Setup

- One-dimensional domain in [0, 1] m with 800 cells and zero-gradient boundaries in the x-direction.
- Flow initialized at rest with the high-pressure water on the left and the low-pressure air on the right. The density ratio between both phases is 20.
- HLLC flux, two-stage SSP Runge-Kutta, vanLeer (MUSCL) reconstruction, inviscid, $\rho$-THINC interface sharpening.
- `./Allrun` runs a single simulation; `validateResults.py` plots the pressure, velocity, density and phase-fraction profiles against the exact Riemann solution.

| State | Phase | Phase fraction | $\rho$ [kg/m3] | $p$ [Pa] | $T$ [K] | $u$ [m/s] |
|-------|-------|----------------|----------------|----------|---------|-----------|
| Left  | Water | 1.0            | 1000           | 1e9      | 781.6   | 0         |
| Right | Air   | 0.0            | 50             | 1e5      | 6.97    | 0         |

## References

1. Saurel, R., Petitpas, F. & Berry, R. A. (2009). Simple and efficient relaxation methods for interfaces separating compressible fluids, cavitating flows and shocks in multiphase mixtures. *Journal of Computational Physics*, 228(5), 1678-1712.
