# Growth Parameters

Parameters used for the Growth stage in Barrientos et al. (2026) *Cell*.

This stage continues from a restart file produced by the Equilibration stage. A constant target stress is applied to the x and y box dimensions via `fix press/berendsen`, allowing the network area to evolve freely (grow or shrink) while stress is maintained. Three experiments are run from a single launcher (`Growth.sh`) by switching the active parameter block.

## Simulation timing

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Time_Eq` | 1×10⁷ | NVT steps for the growth run |
| `NumFrames` | 100 | Trajectory frames saved |
| `Nevery` | 500 | Timesteps between bond reaction attempts |
| `tstep` | 0.01 | Integration timestep (LJ units) |

## Stress control

The box is coupled to a target pressure in x and y independently via `fix press/berendsen`. The z dimension is fixed. Growth is measured as the change in box area (Lx × Ly) over time.

| Parameter | Description |
|-----------|-------------|
| `Press` | Target stress (negative = tensile; units LJ pressure) |
| `Pdamp` = 100.0 | Barendsen damping timescale |

## Experiment A — Vary MP at fixed stress

Probes how bond formation probability affects the growth rate.

| Parameter | Value |
|-----------|-------|
| `MPs` | 0.0001, 0.0002, 0.0004, 0.0008 |
| `Press` | −0.1 |
| `N_preferred` | 5 (all MPs) |
| `Build script` | `build_Growth.py` (GCE + bond breaking active) |

## Experiment B — Vary stress at fixed MP

Probes how the applied tensile stress affects the growth rate.

| Parameter | Value |
|-----------|-------|
| `MP` | 0.0004 |
| `Press` | −0.1, −0.11, −0.12, −0.13, −0.15 |
| `N_preferred` | 5 |
| `Build script` | `build_Growth.py` (GCE + bond breaking active) |

## Experiment C — Three regimes

Compares three mechanistic conditions at the same stress and MP to isolate the roles of bond breaking and monomer addition in growth.

| Regime | Label | `MP` | GCE | Bond breaking | Build script |
|--------|-------|------|-----|---------------|--------------|
| i | No bond breakage | 0.0 | off | off (`BreakProb = MP × FactorMult = 0`) | `build_Growth_NoGCE.py` |
| ii | Bond breakage, no new material | 0.01 | off | on | `build_Growth_NoGCE.py` |
| iii | Bond breakage + new material | 0.01 | on | on | `build_Growth.py` |

`Press = −0.1`, `N_preferred = 5` for all three regimes.

In regime i, setting `MP = 0.0` makes `BreakProb = FactorMult × MP = 0`, so all bond-breaking reactions fire with zero probability — no bonds ever break. In regime ii, breaking is active but GCE is disabled, so no new monomers can be inserted. Regime iii is the full model.

## Common bond kinetics and GCE parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MD` / `MD7s` | 0.65 | NC1 / 7S bond formation distance cutoff |
| `BD` / `BD7s` | 1.35 / 0.95 | NC1 / 7S bond breaking distance cutoff |
| `FactorMult` | 0.66 | Break/make probability ratio |
| `KangNC1` | 4.0 | NC1 angular stiffness |
| `kappa_den` | 0.01 | GCE density penalty strength (Experiments A, B, and regime iii only) |
| `N_preferred` | 5 | Target free monomer count |
| `InitialVol` | 3600 | Reference box volume |

## Replicates

Three independent seeds (1, 2, 3) were run for each condition.
