# Assembly Parameters

Parameters used to generate the assembled networks in Barrientos et al. (2026) *Cell*.

## Simulation timing

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Time_Eq` | 2×10⁷ | NVT timesteps for GCE assembly |
| `NumFrames` | 100 | Trajectory frames saved |
| `Nevery` | 500 | Timesteps between bond reaction attempts |
| `tstep` | 0.01 | Integration timestep (LJ units) |

## Bond kinetics

Both NC1 (end-to-end) and 7S (lateral) bonds use a distance-gated stochastic make/break scheme via `fix bond/react`. The break probability is `FactorMult × MakeProb`.

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MD` | 0.65 | NC1 bond formation distance cutoff |
| `BD` | 1.35 | NC1 bond breaking distance cutoff |
| `MD7s` | 0.65 | 7S bond formation distance cutoff |
| `BD7s` | 0.95 | 7S bond breaking distance cutoff |
| `MP` | 0.12, 0.18, 0.27 | NC1 bond formation probability (three runs) |
| `MP7s` | = MP | 7S bond formation probability (same as NC1) |
| `FactorMult` | 0.66 | Break/make probability ratio |

## Mechanical parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `KangNC1` | 4.0 | NC1 angular stiffness (harmonic, degrees from 180°) |
| `k_nc1` | 6.0 | NC1 bond spring constant |
| `r0_nc1` | 0.5 | NC1 bond rest length |
| `k_7s` | 6.0 | 7S bond spring constant |
| `r0_7s` | 0.5 | 7S bond rest length |
| `max_bonds_7s` | 3 | Maximum 7S bonds per bead |

## Grand canonical ensemble (GCE)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `N_preferred` | 0 | Target monomer count (0 = fixed number) |
| `OvaR` | 3.5 | Overlap exclusion radius for monomer insertion |

## Replicates

Three independent seeds (1, 2, 3) were run for each `MP` value, giving 9 simulations in total.

## Initial configuration

`input/dataNucType9_dense` — a pre-equilibrated random configuration of collagen IV protomers at the dense (physiological) density, with no bonds formed.
