# Bond Sampling Parameters

Parameters used for the Bond Sampling stage in Barrientos et al. (2026) *Cell*.

This stage continues from a restart file produced by the Equilibration stage and runs a long NVT simulation with the same bond/react fix. Rather than equilibrating the network, the purpose is to accumulate statistics on how many bond formation and breaking events occur over time — i.e. the molecular turnover rate of the network.

## Simulation timing

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Time_Eq` | 2×10⁷ | NVT timesteps for the bond sampling run |
| `NumFrames` | 50 | Trajectory frames saved |
| `Nevery` | 500 | Timesteps between bond reaction attempts |
| `tstep` | 0.01 | Integration timestep (LJ units) |

## Bond kinetics and GCE parameters

Identical to the Equilibration stage. GCE (Nucleation and Death reactions) remains active during bond sampling so that the monomer pool is maintained at the equilibrium level.

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MD` / `MD7s` | 0.65 | NC1 / 7S bond formation distance cutoff |
| `BD` / `BD7s` | 1.35 / 0.95 | NC1 / 7S bond breaking distance cutoff |
| `MP` | 0.12, 0.18, 0.27 | NC1 bond formation probability |
| `FactorMult` | 0.66 | Break/make probability ratio |
| `KangNC1` | 4.0 | NC1 angular stiffness |
| `kappa_den` | 0.01 | GCE density penalty strength |
| `N_preferred` | 8, 8, 8 | Target free monomer count for MP = 0.12, 0.18, 0.27 |

## Bond counting

The key output of this stage is cumulative bond event counts, tracked via `f_freact_creation` indices and written to `thermo.dat`:

- **`NumRelease`**: total bond-breaking events (NC1 + 7S), counting each broken bond pair once.
  Built from breaking-reaction counters for all NC1 variants (2N–8N, loop) and 7S variants.
- **`NumIncorp`**: total bond-forming events (NC1 + 7S), counting each formed bond pair once.
  Built from formation-reaction counters for all NC1 and 7S variants.

These cumulative counts grow monotonically. The turnover rate is extracted from their slope over time.

## Replicates

Three independent seeds (1, 2, 3) were run for each `MP` value, giving 9 simulations in total.
