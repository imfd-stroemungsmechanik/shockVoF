# Shock-Water-Column Interaction at Ma = 1.47

Shock-induced deformation of a 2D water column behind a planar shock in air at Ma = 1.47 according to the experiments of Igra & Takayama (2001).

## Setup

- Two-dimensional square domain of size 40 $D_0$ x 40 $D_0$. The graded 1200 x 1200 structured mesh resolves the droplet at 120 cells per diameter and coarsens toward the boundaries. Left and right are non-reflecting boundaries, top and bottom are symmetry planes.
- The water column is initialized at rest in the center of the domain with a diameter of $D_0$ = 4.8 mm and a diffuse interface consistent with the $\rho$-THINC interface sharpening, which spans the same number of cells on every mesh. The shock is located 4.8 mm upstream of the water column center with post-shock conditions following the Rankine-Hugoniot solution for a Ma = 1.47 shock in air.
- HLLC flux, three-stage SSP Runge-Kutta, vanLeer reconstruction, viscous, $\rho$-THINC interface sharpening, smoothed curvature model.
- `./Allrun` creates the mesh, decomposes and runs the case in parallel; `validateResults.py` plots the leading-edge drift and the centerline width against the experimental data.

| Parameter               | Symbol            | Value      |
|-------------------------|-------------------|------------|
| Shock Mach number       | Ma                | 1.47       |
| Initial column diameter | $D_0$             | 4.8 mm     |
| Post-shock gas velocity | $u_g$             | 225.9 m/s  |
| Water density           | $\rho_\text{liq}$ | 1000 kg/m3 |
| Surface tension         | $\sigma$          | 0.0774 N/m |
| Weber number            | We                | 6900       |
| Reynolds number         | Re                | 1.07e5     |

The surface tension is set so the crossflow Weber number matches the experimental Weber number of 6900. Time is non-dimensionalized by the Ranger & Nicholls (1969) scaling

$$ t^* = (t - t_0)\,\frac{u_g}{D_0}\sqrt{\frac{\rho_g}{\rho_l}} $$

with the post-shock gas velocity $u_g$, initial diameter $D_0$ and post-shock gas and liquid densities $\rho_g$ and $\rho_l$.

## References

1. Igra, D. & Takayama, K. (2001). Investigation of aerodynamic breakup of a cylindrical water droplet. *Atomization and Sprays*, 11(2), 167-185.
2. Ranger, A. A. & Nicholls, J. A. (1969). Aerodynamic shattering of liquid drops. *AIAA Journal*, 7(2), 285-290.
