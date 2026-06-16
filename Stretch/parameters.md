# Stretch Parameters

Parameters used for the Stretch stage in Barrientos et al. (2026) *Cell*.

This stage continues from a restart file produced by the Equilibration stage. A biaxial equibiaxial strain is applied to the network in incremental steps, then the network is held at fixed strain and allowed to relax. The main readout is the network stress time series.

## Simulation timing

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Time_Eq` | 0 | No additional NVT before stretching (system is already equilibrated) |
| `Time_Str` | 7.5×10⁴ | NVT steps for the stretch phase (divided across 20 incremental steps) |
| `Time_Rel` | 1×10⁶ | NVT steps for stress relaxation at fixed strain |
| `NumFrames` | 100 | Trajectory frames saved |
| `Nevery` | 5000 | Timesteps between bond reaction attempts |
| `tstep` | 0.01 | Integration timestep (LJ units) |

## Strain protocol

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Xstretch` | 50 | Total box length increase applied in each direction (LJ units) |
| `Nstretch` | 20 | Number of incremental stretch steps |

The stretch is applied biaxially: the 20 steps alternate between expanding the box in x and expanding in y. Each step uses `fix deform` to shift the box boundaries by ±`Xstretch / (2 × Nstretch)` = ±1.25 LJ units on each side, for a total extension of 50 per direction. After all 20 steps the box is held fixed and the relaxation run begins.

## Bond kinetics and GCE parameters

Identical to the Equilibration stage. GCE remains active throughout so the monomer pool is maintained at the equilibrium level.

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MD` / `MD7s` | 0.65 | NC1 / 7S bond formation distance cutoff |
| `BD` / `BD7s` | 1.35 / 0.95 | NC1 / 7S bond breaking distance cutoff |
| `MP` | 0.12, 0.18, 0.27 | NC1 bond formation probability |
| `FactorMult` | 0.66 | Break/make probability ratio |
| `KangNC1` | 4.0 | NC1 angular stiffness |
| `kappa_den` | 0.01 | GCE density penalty strength |
| `N_preferred` | 8, 8, 8 | Target free monomer count for MP = 0.12, 0.18, 0.27 |

## Replicates

Three independent seeds (1, 2, 3) were run for each `MP` value, giving 9 simulations in total.
