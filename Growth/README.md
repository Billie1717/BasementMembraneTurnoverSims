# Growth

Scripts for the constant-stress growth stage of the collagen IV basement membrane model.

Starting from a restart file produced by the Equilibration stage, this stage applies a constant tensile stress to the x and y box dimensions via `fix press/berendsen`. The box area grows freely to maintain the target stress, producing indefinite network growth. Three experiments are run from a single launcher (`Growth.sh`) by switching the active parameter block.

## Experiments

**Experiment A — Vary MP at fixed stress.** The NC1 bond formation probability (`MP`) is varied across four values at a fixed applied stress (`Press = −0.1`). This probes how bond kinetics set the growth rate.

**Experiment B — Vary stress at fixed MP.** The applied stress is varied across five values at a fixed `MP = 0.0004`. This probes how the mechanical driving force sets the growth rate.

**Experiment C — Three regimes.** At fixed stress and fixed MP, three mechanistic conditions are compared to isolate the roles of bond breaking and monomer addition:

| Regime | Condition | `MP` | GCE | Bond breaking |
|--------|-----------|------|-----|---------------|
| i | No bond breakage | 0.0 | off | off |
| ii | Bond breakage, no new material | 0.01 | off | on |
| iii | Bond breakage + new material | 0.01 | on | on |

In regime i, `MP = 0.0` makes `BreakProb = FactorMult × MP = 0`, so no bonds ever break. In regime ii, breaking is active but the GCE is disabled so no new monomers can enter. Regime iii is the full model.

## Build scripts

Two build scripts are provided:

- **`build_Growth.py`** — GCE (Nucleation + Death) active; bond counting (`NumRelease`, `NumIncorp`) included. Used for Experiments A, B, and Experiment C regime iii.
- **`build_Growth_NoGCE.py`** — GCE disabled (Nucleation and Death reactions commented out). Used for Experiment C regimes i and ii.

For a detailed walkthrough of both scripts, see `build_Growth_annotated.py`. For the exact parameters used in the paper, see `parameters.md`.

## Files

| File | Description |
|------|-------------|
| `Growth.sh` | Top-level launcher covering all three experiments via comments |
| `build_Growth.py` | Writes `collagen.in` and `runscript.sh` with GCE active |
| `build_Growth_NoGCE.py` | Writes `collagen.in` and `runscript.sh` with GCE disabled |
| `build_Growth_annotated.py` | Non-runnable annotated walkthrough of both build scripts |
| `parameters.md` | Parameters used for the paper |
| `templates/` | `bond/react` template files (same set as Equilibration) |
| `analysis/Histogram_Types.py` | Counts beads of each type at each trajectory frame |

## Running the simulations

**Requirements:** Python ≥ 3.9, LAMMPS `lammps-15Jun2023` built with the custom `REACTION` package (see top-level README). The Equilibration stage must be complete.

1. Edit `Growth.sh` to select the experiment (A, B, or C) by commenting/uncommenting the appropriate parameter block.
2. Edit `Growth.sh` to set the path to your Equilibration `runs/` directory (default: `../Equilibration/runs`).
3. Edit the relevant build script (Section 11) to set the path to your LAMMPS executable inside `runscript.sh`.
4. From this directory, run:
   ```bash
   bash Growth.sh
   ```
   The script automatically finds the latest restart file from each Equilibration run and copies it as the starting configuration. Output is written to `runs/`.
5. Each job is submitted via `sbatch`. On completion, each run folder contains:
   - `dumplin/` — particle trajectory (`.lammpstrj` files)
   - `dumplin_bonds/` — bond topology at each frame
   - `thermo.dat` — thermodynamic time series including box dimensions and stress
   - `restart/` — restart files
   - `data.postmakebonds` — final configuration with bond topology
6. One can also directly run using `mpirun -np 1 /PATH/lammps-15Jun2023/src/lmp_mpi -in collagen.in` from within the run folder.

## Output columns of `thermo.dat`

### `build_Growth.py` (Experiments A, B, regime iii)

| Column | Description |
|--------|-------------|
| step | Timestep |
| etot | Total energy |
| ke | Kinetic energy |
| peBond | Bond potential energy |
| peAngle | Angle potential energy |
| temp | Temperature |
| press | Pressure |
| pxx / pyy | Per-direction pressure |
| stressX/Y/Z | Total virial stress (all atoms) |
| stressX_N/Y_N/Z_N | Network-only virial stress |
| sigma_yy / sigma_xx | Network stress normalised by volume (primary observable) |
| FtotX / FtotY | Total correction force applied to network |
| AvBondForce | Mean bond force magnitude |
| NumMons | Number of free monomers (type 9) |
| BoxLength | Box length in x |
| BoxVol | Box volume (primary growth observable — grows under tensile stress) |
| CurrentConc | Current free-monomer concentration |
| ConcPref | Target free-monomer concentration |
| P_add | GCE insertion acceptance probability |
| P_remov | GCE deletion acceptance probability |
| NreactNuc | Cumulative successful insertions |
| NreactDeath | Cumulative successful deletions |
| T_den_add | Boltzmann factor for insertion |
| T_den_remov | Boltzmann factor for deletion |
| T_id_add | Ideal-gas factor for insertion |
| T_id_remov | Ideal-gas factor for deletion |
| NumTot | Total protomer count (bound + free) |
| **NumRelease** | Cumulative bond-breaking events (all types, weighted) |
| **NumIncorp** | Cumulative bond-forming events (all types, weighted) |

### `build_Growth_NoGCE.py` (Experiment C regimes i and ii)

Same columns as above except `NumRelease` and `NumIncorp` are absent. The GCE columns (`P_add`, `P_remov`, etc.) are still written but are unused since GCE reactions are disabled.

The primary observable is `BoxVol` — its increase over time gives the area growth rate.
