# Liquid Jet Atomization in a Supersonic Crossflow at Ma = 2.0

Atomization of a round water jet injected into a Mach 2.0 air crossflow according to the experiments of Lin et al. (2020).

## Setup

- Cubic domain with a length of 24 $d_0$ and the injector centered at the bottom wall 12 $d_0$ downstream of the gas inlet. The background mesh has N = 60 cells in each direction and is refined by four levels using `snappyHexMesh` around the injector, the bottom wall and the primary atomization zone. The jet is resolved with 40 cells per diameter ($\Delta x = 19\,\text{µm}$).
- The subgrid-scale turbulence of both phases is modeled by a single mixture WALE LES model (Nicoud & Ducros 1999) with Spalding's wall function at the bottom wall.
- HLLC flux, three-stage SSP Runge-Kutta, vanLeer reconstruction, viscous, $\rho$-THINC interface sharpening, smoothed curvature model.
- The gas inlet carries a uniform velocity with a 1/7th power-law turbulent boundary layer of 1.5 mm thickness at the bottom wall. The water inlet is a 1/7th power-law pipe profile with superimposed divergence-free synthetic turbulence (Smirnov et al. 2001).
- `./Allrun` creates the mesh, decomposes and runs the case in parallel, and samples the horizontal planes of `system/trajectory` at the latest time; from these planes, `validateResults.py` plots the trajectory and the equivalent path length (EPL) profiles against the experimental data. As for the experimental data, the trajectory is the location of the maximum EPL at each height.
- `system/configDict` enables the debug switches of `fluxPredictor` and `reconstructionScheme`, which write the diagnostic fields `fluxLimiter` and `positivityLimiter`.

## Flow Conditions

| Parameter             | Symbol | Value     |
|-----------------------|--------|-----------|
| Jet diameter          | $d_0$  | 0.76 mm   |
| Crossflow Mach number | Ma     | 2.0       |
| Momentum flux ratio   | q      | 4.09      |
| Gas inlet pressure    | $p_g$  | 39.2 kPa  |
| Gas inlet temperature | $T_g$  | 166.7 K   |
| Gas inlet velocity    | $u_g$  | 517.8 m/s |
| Liquid inlet velocity | $u_l$  | 30 m/s    |

## References

1. Lin, K.-C., Kastengren, A. L., Hammack, S. & Carter, C. (2020). Exploration of water jets in supersonic crossflow using X-ray diagnostics. *Atomization and Sprays*, 30(5), 331-350.
2. Nicoud, F. & Ducros, F. (1999). Subgrid-scale stress modelling based on the square of the velocity gradient tensor. *Flow, Turbulence and Combustion*, 62(3), 183-200.
3. Smirnov, A., Shi, S. & Celik, I. (2001). Random flow generation technique for large eddy simulations and particle-dynamics modeling. *Journal of Fluids Engineering*, 123(2), 359-371.
