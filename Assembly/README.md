# Assembly

Scripts for the self-assembly stage of the collagen IV basement membrane model.

Starting from a random configuration of unbonded protomers, the network assembles through stochastic bond formation and breaking at both NC1 and 7S domains until it reaches a percolating, equilibrated state. The assembled network is the input for all subsequent simulation stages (Equilibration with GCE, Stretch, Growth).

## Model overview

Each collagen IV protomer is represented as a two-bead rod:

- **NC1 bead** (type 1): end-to-end domain. Bonds pairwise with other NC1 beads. The NC1–NC1 bond has a harmonic angle potential that enforces near-linearity (rest angle 180°).
- **7S bead** (types 2–8): lateral domain. Bonds with up to three other 7S beads, forming dimers, trimers, and tetramers. The bead type encodes how many 7S bonds it currently has (type 2 = 1 bond, type 3 = 2 bonds, etc.), which is why multiple template files are needed.
- **Type 9**: free (unbound) NC1 bead, i.e. a protomer not yet incorporated into the network.

Bonds form and break stochastically via LAMMPS' `fix bond/react` package.

For a full walkthrough of the LAMMPS input deck, see `build_Assembly_annotated.py`. For the exact parameters used in the paper, see `parameters.md`.

## Files

| File | Description |
|------|-------------|
| `Assemble.sh` | Top-level launcher; loops over parameters and submits jobs |
| `build_Assembly.py` | Writes `collagen.in` (LAMMPS input) and `runscript.sh` (SLURM script) |
| `build_Assembly_annotated.py` | Non-runnable annotated walkthrough of `build_Assembly.py` |
| `parameters.md` | Parameters used for the paper |
| `input/dataNucType9_dense` | Initial configuration (random, no bonds) at physiological density |
| `templates/` | `bond/react` template files for all allowed bond topology changes |
| `analysis/Histogram_Types.py` | Counts beads of each type at each trajectory frame |
| `analysis/HistogrammingTypesAssemble.sh` | Launcher for the histogramming analysis |

## Running the simulations

**Requirements:** Python ≥ 3.9, LAMMPS `lammps-15Jun2023` built with the custom `REACTION` package (see top-level README).

1. Edit `Assemble.sh` to set the path to your LAMMPS executable inside `runscript.sh` (written by `build_Assembly.py`, see Section 10 of the annotated script).
2. From this directory, run:
   ```bash
   bash Assemble.sh
   ```
   This creates a `runs/` subdirectory containing one folder per parameter combination, each with `collagen.in`, `runscript.sh`, and a copy of the input data file.
3. Each job is submitted via `sbatch`. On completion, each run folder contains:
   - `dumplin/` — particle trajectory (`.lammpstrj` files)
   - `dumplin_bonds/` — bond topology at each frame
   - `thermo.dat` — thermodynamic time series
   - `restart/` — restart files
   - `data.postmakebonds` — final configuration with bond topology
4. One can also directly run using mpirun -np 1 /PATH/lammps-15Jun2023/src/lmp_mpi -in collagen.in from within the folder

## Analysis: histogramming particle types

After simulations complete, run the histogramming script to count bead types at each frame. This gives a time series of network connectivity (number of free monomers, dimers, trimers, etc.) used to assess equilibration.

```bash
bash analysis/HistogrammingTypesAssemble.sh
```

Output `.txt` files are written to `analysis/Data/`. These are read by the plotting notebooks (see `analysis/Plotting/`).

### Output columns of `thermo.dat`

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
| CurrentConc | Current monomer concentration |
| ConcPref | Reference monomer concentration |
