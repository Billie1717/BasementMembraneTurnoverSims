# Equilibration Parameters

Parameters used for the Equilibration with GCE stage in Barrientos et al. (2026) *Cell*.

This stage continues from a restart file produced by the Assembly stage and introduces active monomer exchange via the Grand Canonical Ensemble (GCE).

## Simulation timing

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Time_Eq` | 3×10⁶ | NVT timesteps for GCE equilibration |
| `NumFrames` | 100 | Trajectory frames saved |
| `Nevery` | 500 | Timesteps between bond reaction attempts |
| `tstep` | 0.01 | Integration timestep (LJ units) |

## Bond kinetics

Identical to the Assembly stage.

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MD` | 0.65 | NC1 bond formation distance cutoff |
| `BD` | 1.35 | NC1 bond breaking distance cutoff |
| `MD7s` | 0.65 | 7S bond formation distance cutoff |
| `BD7s` | 0.95 | 7S bond breaking distance cutoff |
| `MP` | 0.12, 0.18, 0.27 | NC1 bond formation probability (three runs) |
| `MP7s` | = MP | 7S bond formation probability |
| `FactorMult` | 0.66 | Break/make probability ratio |

## Mechanical parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `KangNC1` | 4.0 | NC1 angular stiffness (harmonic, rest angle 180°) |

## Grand Canonical Exchange (GCE)

Monomers (type-9 free NC1 beads) are inserted and deleted stochastically each timestep via `bond/react` Nucleation and Death reactions. The acceptance probability follows a soft harmonic density penalty:

$$\phi = \kappa_\rho \left(C_\mathrm{tot} V - C_\mathrm{tot}^\mathrm{target} V\right)^2$$

where $C_\mathrm{tot}$ is the current total protomer concentration (bound + free), $C_\mathrm{tot}^\mathrm{target}$ is fixed by `InitialVol` (the box volume at the start of Assembly), and $\kappa_\rho$ controls how tightly the monomer number is regulated.

Acceptance probabilities:

$$\mathrm{Acc}(N \to N+1) = \min\!\left(1,\, \frac{C_H V}{N_H+1}\exp\!\left[\frac{\phi(N_H)-\phi(N_H+1)}{k_B T}\right]\right)$$

$$\mathrm{Acc}(N \to N-1) = \min\!\left(1,\, \frac{N_H}{C_H V}\exp\!\left[\frac{\phi(N_H)-\phi(N_H-1)}{k_B T}\right]\right)$$

where $C_H$ is the preferred free-monomer concentration (`N_preferred / InitialVol`) and $N_H$ is the current free-monomer count.

| Parameter | Value | Description |
|-----------|-------|-------------|
| `kappa_den` | 0.01 | Density penalty strength ($\kappa_\rho$) |
| `InitialVol` | 12500 | Reference box volume (dense density) |
| `N_preferred` | 8, 8, 8 | Target free monomer count for MP = 0.12, 0.18, 0.27 |

`N_preferred` is chosen as the mean number of free monomers (type 9) at the end of the Assembly stage, read from the `NumMons` column of `thermo.dat` once the assembly has plateaued. This anchors the GCE to the equilibrium monomer pool already established by assembly, so the GCE exchange fluctuates around a physically meaningful baseline rather than driving the system to an arbitrary density.

## Replicates

Three independent seeds (1, 2, 3) were run for each `MP` value.
