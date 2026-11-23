# Parallel Scaling: 3D Shock-Droplet Interaction

Strong and weak scaling of `shockVoF` for a three-dimensional shock-droplet interaction on 1 to 64 nodes with 104 cores each.

## Setup

- Water droplet of diameter $D_0$ = 4.8 mm at rest in air, hit by a planar Ma = 1.47 shock with the post-shock state of the [`Igra2001`](../validation/Igra2001) case. The shock starts at the leading edge of the droplet, 2.4 mm upstream of its center. Each droplet occupies a cube with an edge length of 4 $D_0$. The domain has non-reflecting boundaries in the x-direction and symmetry planes on the four lateral sides.
- **Strong scaling:** a single droplet on a fixed uniform mesh of 256 x 256 x 256 = 16.8 M cells (64 cells per $D_0$).
- **Weak scaling:** 128 x 128 x 128 = 2.1 M cells (32 cells per $D_0$) and one droplet per node. Each doubling of the node count doubles the domain and the number of droplets.
- HLLC flux, three-stage SSP Runge-Kutta, vanLeer reconstruction, viscous, $\rho$-THINC interface sharpening, smoothed curvature model.
- Every run advances 500 time steps of fixed size $\Delta t$ = 5 ns without function objects or field output, so only the solution algorithm and the parallel communication are timed. The mesh is decomposed with `scotch`, with one MPI rank per physical core.
- `./Allpre` creates and decomposes the cases in `strongScaling` and `weakScaling`; the cluster-specific job scripts are not included. `evaluateScaling.py` reads the timings from the `slurm-*.out` logs in the run directories and plots the parallel efficiencies.

| Parameter      | Strong scaling      | Weak scaling        |
|----------------|---------------------|---------------------|
| Nodes (cores)  | 1 - 64 (104 - 6656) | 1 - 64 (104 - 6656) |
| Total cells    | 16.8 M              | 2.1 M - 134.2 M     |
| Cells per core | 161.3 k - 2.5 k     | 20.2 k              |
| Droplets       | 1                   | 1 - 64              |
| Time steps     | 500                 | 500                 |
| Time on 1 node | 1317 s              | 137 s               |

## Performance Metrics

The execution time $T_N$ on $N$ nodes is the wall-clock time (`ClockTime`) after the last time step minus that after the first, so it excludes the initialization and covers $n_\text{steps}$ = 499 time steps; for repeated runs the median is taken. With the single-node run as reference, the strong and weak scaling efficiencies are

$$ E_\text{strong} = \frac{T_1}{N \, T_N}, \qquad E_\text{weak} = \frac{T_1}{T_N} ~. $$

The performance index (PID) is the core time per cell and time step:

$$ \text{PID} = \frac{T_N \, n_\text{cores}}{n_\text{cells} \, n_\text{steps}} ~. $$

`evaluateScaling.py` prints $T_N$, the efficiencies and the PID for each node count.
