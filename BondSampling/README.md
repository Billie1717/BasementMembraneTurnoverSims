# Bond Sampling

Scripts for measuring molecular turnover in the collagen IV basement membrane model.

Starting from a restart file produced by the Equilibration stage, this stage runs a long NVT simulation with the same bond/react fix and GCE as Equilibration. The key addition is a set of cumulative bond event counters that track how many bonds form and break over the course of the run — quantifying the molecular turnover rate of the network at steady state.

## Model overview

The physics, force field, and GCE are identical to the Equilibration stage. The only additions are LAMMPS variables that accumulate counts of bond formation (`NumIncorp`) and bond breaking (`NumRelease`) events across all reaction types. These are built from the `f_freact_creation` per-reaction counters exposed by `fix bond/react`.

The run is longer (2×10⁷ steps) and saves fewer trajectory frames (50) since the readout is a time series of cumulative event counts rather than structural snapshots.

For a detailed walkthrough of the counting variables, see `build_BondSampling_annotated.py`. For parameters, see `parameters.md`.

## Files

| File | Description |
|------|-------------|
| `BondSampling.sh` | Top-level launcher; locates Equilibration restart files and submits jobs |
| `build_BondSampling.py` | Writes `collagen.in` (LAMMPS input) and `runscript.sh` (SLURM script) |
| `build_BondSampling_annotated.py` | Non-runnable annotated walkthrough focused on the bond counting variables |
| `parameters.md` | Parameters used for the paper |
| `templates/` | `bond/react` template files (same set as Equilibration) |
| `analysis/Histogram_Types.py` | Counts beads of each type at each trajectory frame |
| `analysis/HistogrammingTypesBondSampling.sh` | Launcher for the histogramming analysis |

## Running the simulations

**Requirements:** Python ≥ 3.9, LAMMPS `lammps-15Jun2023` built with the custom `REACTION` package (see top-level README). The Equilibration stage must be complete.

1. Edit `BondSampling.sh` to set the path to your Equilibration `runs/` directory (default: `../Equilibration/runs`).
2. Edit `build_BondSampling.py` Section 10 to set the path to your LAMMPS executable inside `runscript.sh`.
3. From this directory, run:
   ```bash
   bash BondSampling.sh
   ```
   The script automatically finds the latest restart file from each Equilibration run and copies it as the starting configuration. Output is written to `runs/`.
4. Each job is submitted via `sbatch`. On completion, each run folder contains:
   - `dumplin/` — particle trajectory (`.lammpstrj` files)
   - `dumplin_bonds/` — bond topology at each frame
   - `thermo.dat` — thermodynamic time series including bond event counts
   - `restart/` — restart files
   - `data.postmakebonds` — final configuration with bond topology
5. One can also directly run using `mpirun -np 1 /PATH/lammps-15Jun2023/src/lmp_mpi -in collagen.in` from within the run folder.

## Analysis: histogramming particle types

After simulations complete, run the histogramming script to count bead types at each frame.

```bash
bash analysis/HistogrammingTypesBondSampling.sh
```

Output `.txt` files are written to `analysis/Data/`. These are read by the plotting notebooks (see `analysis/Plotting/`).

## Output columns of `thermo.dat`

All columns from the Equilibration stage are present, plus two additional columns at the end:

| Column | Description |
|--------|-------------|
| step | Timestep |
| etot | Total energy |
| ke | Kinetic energy |
| peBond | Bond potential energy |
| peAngle | Angle potential energy |
| temp | Temperature |
| press | Pressure |
| stressX/Y/Z | Total virial stress in each direction |
| AvBondForce | Mean bond force magnitude |
| NumMons | Number of free monomers (type 9) |
| BoxLength | Box length in x |
| BoxVol | Box volume |
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
| **NumRelease** | Cumulative bond-breaking events (all NC1 and 7S types, weighted) |
| **NumIncorp** | Cumulative bond-forming events (all NC1 and 7S types, weighted) |

`NumRelease` and `NumIncorp` grow monotonically. The molecular turnover rate is extracted from the slope of either quantity versus time (at steady state they are equal).
