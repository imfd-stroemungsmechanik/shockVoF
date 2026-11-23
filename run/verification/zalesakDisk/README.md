# Zalesak's Slotted Disk

A slotted disk is rotated by a solid-body-rotation velocity field (Zalesak 1979). After one period the exact solution returns to the initial condition, so the initial data are the analytical reference and any smearing or corner rounding is directly measurable. Here the case isolates the effect of $\rho$-THINC interface sharpening by advecting the disk with and without it.

## Setup

- Two-dimensional square with $L$ = 1 m resolved by a structured mesh with N x N cells.
- Solid-body rotation about the center (0.5, 0.5) at $\omega$ = 50.27 rad/s. One revolution takes 0.125 s.
- Slotted disk centered at (0.5, 0.75), radius 0.15 m, with a slot of width 0.05 m cut from the bottom of the disk, toward the rotation center, up to y = 0.85 m.
- The flow is frozen: the prescribed velocity field is kept constant and only the volume fraction is advected.
- HLLC flux, two-stage SSP Runge-Kutta, vanLeer (MUSCL) reconstruction, inviscid, with and without $\rho$-THINC interface sharpening.
- `./Allrun` runs the mesh sweep N = 32, 64, 128, 256 with and without $\rho$-THINC; `validateResults.py` plots the $\alpha = 0.5$ isocontours of the disk and the phase-fraction profile through the slot against the exact solution.

## References

1. Zalesak, S. T. (1979). Fully multidimensional flux-corrected transport algorithms for fluids. *Journal of Computational Physics*, 31(3), 335-362.
