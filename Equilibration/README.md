# Equilibration with GCE

Scripts for the GCE equilibration stage of the collagen IV basement membrane model.

Starting from a restart file produced by the Assembly stage, this stage runs MD under NVT dynamics coupled with MC with an active Grand Canonical Exchange (GCE): type-9 free monomers are inserted and deleted stochastically to bring the system to a prescribed target monomer concentration. Equilibrated configurations from this stage are the input for the Stretch and Growth stages.

## Model overview

The force field is identical to the Assembly stage (see `Assembly/README.md`). The key addition is the GCE, which operates via two `fix bond/react` reactions:

- **Nucleation**: inserts a new type-9 free monomer at a randomly chosen overlap-free site (determined by the neighbour list cutoff via `modify_create nuc yes`).
- **Death**: removes an existing type-9 free monomer.

Each insertion or deletion is accepted or rejected stochastically using a soft harmonic density penalty (κ_ρ) that drives the total protomer count toward a target value `N_preferred`. See `parameters.md` for the full acceptance equations.

For a detailed walkthrough of the LAMMPS input deck, see `build_Equilibration_annotated.py`. For the exact parameters used in the paper, see `parameters.md`.

## Files

| File | Description |
|------|-------------|
| `Equilibrate.sh` | Top-level launcher; locates Assembly restart files and submits jobs |
| `build_Equilibration.py` | Writes `collagen.in` (LAMMPS input) and `runscript.sh` (SLURM script) |
| `build_Equilibration_annotated.py` | Non-runnable annotated walkthrough of `build_Equilibration.py` |
| `parameters.md` | Parameters used for the paper |
| `templates/` | `bond/react` template files (same set as Assembly) |
| `analysis/Histogram_Types.py` | Counts beads of each type at each trajectory frame |
| `analysis/HistogrammingTypesEquib.sh` | Launcher for the histogramming analysis |

## Running the simulations

**Requirements:** Python ≥ 3.9, LAMMPS `lammps-15Jun2023` built with the custom `REACTION` package (see top-level README). The Assembly stage must be complete.

1. Edit `Equilibrate.sh` to set the path to your Assembly `runs/` directory (default: `../Assembly/runs`).
2. Edit `build_Equilibration.py` Section 10 to set the path to your LAMMPS executable inside `runscript.sh`.
3. From this directory, run:
   ```bash
   bash Equilibrate.sh
   ```
   The script automatically finds the latest restart file from each Assembly run and copies it as the starting configuration. Output is written to `runs/`.
4. Each job is submitted via `sbatch`. On completion, each run folder contains:
   - `dumplin/` — particle trajectory (`.lammpstrj` files)
   - `dumplin_bonds/` — bond topology at each frame
   - `thermo.dat` — thermodynamic time series
   - `restart/` — restart files
   - `data.postmakebonds` — final configuration with bond topology
5. One can also directly run using `mpirun -np 1 /PATH/lammps-15Jun2023/src/lmp_mpi -in collagen.in` from within the run folder.

## Analysis: histogramming particle types

After simulations complete, run the histogramming script to count bead types at each frame.

```bash
bash analysis/HistogrammingTypesEquib.sh
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
