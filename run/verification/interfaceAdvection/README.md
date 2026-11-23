# Interface Advection

A water-air material interface advected at uniform velocity, pressure and temperature (Abgrall 1996). The exact solution keeps p, T and u uniform for all time while only the volume fraction and density translate. Any pressure or velocity deviation is a spurious oscillation from incorrect coupling between the non-conservative volume-fraction equation and the conservative phasic system.

## Setup

- One-dimensional domain in [0, 1] m with 200 cells and cyclic boundaries in the x-direction.
- Flow initialized as uniform with a volume-fraction discontinuity between 0.25 < $x$ < 0.75 m. With $u$ = 100 m/s the interface completes exactly one period at 0.01 s and returns to its start.
- HLLC flux, two-stage SSP Runge-Kutta, vanLeer (MUSCL) reconstruction, inviscid, $\rho$-THINC interface sharpening.
- `./Allrun` runs a single simulation; `validateResults.py` plots the pressure and velocity errors and the phase-fraction profile.

| Parameter      | Symbol   | Value                           |
|----------------|----------|---------------------------------|
| Pressure       | $p$      | 101325 Pa                       |
| Temperature    | $T$      | 293.15 K                        |
| Velocity       | $u$      | 100 m/s                         |
| Phase fraction | $\alpha$ | 1 for 0.25 < $x$ < 0.75, else 0 |

## References

1. Abgrall, R. (1996). How to prevent pressure oscillations in multicomponent flow calculations: A quasi conservative approach. *Journal of Computational Physics*, 125(1), 150-160.
