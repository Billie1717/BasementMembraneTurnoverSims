# Stretch

Scripts for the stretch-and-relax stage of the collagen IV basement membrane model.

Starting from a restart file produced by the Equilibration stage, this stage applies a biaxial strain to the network in 20 incremental steps, then holds the strain fixed and runs a long NVT + MC relaxation. The main readout is the network stress time series in `thermo.dat`, which captures how the network stress evolves as bonds rearrange at the imposed strain.

## Protocol

The simulation proceeds in three phases:

**Stretch phase** (`Time_Str = 7.5×10⁴` steps): The simulation box is expanded biaxially in 20 incremental steps, alternating between the x and y directions. Each step uses `fix deform` to extend the box by 2.5 LJ units (1.25 on each side), for a total extension of 50 LJ units in both x and y. Atoms move affinely with the box (`remap x`).

**Short relaxation** (`tInitRel = Time_Str/2` steps): The bond/react fix is reset to clear internal counters accumulated during the stretch, then a brief NVT run re-equilibrates the bond network at the new geometry.

**Long relaxation** (`Time_Rel = 1×10⁶` steps): NVT run at fixed strain. This is the main measurement window — the network stress decays as bonds break and reform.

The GCE (Nucleation and Death reactions) remains active throughout, maintaining the monomer pool at the equilibrium level.

For a detailed walkthrough of the LAMMPS input deck, see `build_Stretch_annotated.py`. For the exact parameters used in the paper, see `parameters.md`.

## Files

| File | Description |
|------|-------------|
| `Stretch.sh` | Top-level launcher; locates Equilibration restart files and submits jobs |
| `build_Stretch.py` | Writes `collagen.in` (LAMMPS input) and `runscript.sh` (SLURM script) |
| `build_Stretch_annotated.py` | Non-runnable annotated walkthrough of `build_Stretch.py` |
| `parameters.md` | Parameters used for the paper |
| `templates/` | `bond/react` template files (same set as Equilibration) |
| `analysis/Histogram_Types.py` | Counts beads of each type at each trajectory frame |

## Running the simulations

**Requirements:** Python ≥ 3.9, LAMMPS `lammps-15Jun2023` built with the custom `REACTION` package (see top-level README). The Equilibration stage must be complete.

1. Edit `Stretch.sh` to set the path to your Equilibration `runs/` directory (default: `../Equilibration/runs`).
2. Edit `build_Stretch.py` Section 11 to set the path to your LAMMPS executable inside `runscript.sh`.
3. From this directory, run:
   ```bash
   bash Stretch.sh
   ```
   The script automatically finds the latest restart file from each Equilibration run and copies it as the starting configuration. Output is written to `runs/`.
4. Each job is submitted via `sbatch`. On completion, each run folder contains:
   - `dumplin/` — particle trajectory (`.lammpstrj` files)
   - `dumplin_bonds/` — bond topology at each frame
   - `thermo.dat` — thermodynamic and stress time series
   - `restart/` — restart files
   - `data.poststretch` — configuration after the stretch phase
   - `data.postrelax` — configuration after the relaxation phase
5. One can also directly run using `mpirun -np 1 /PATH/lammps-15Jun2023/src/lmp_mpi -in collagen.in` from within the run folder.

## Output columns of `thermo.dat`

| Column | Description |
|--------|-------------|
| step | Timestep |
| etot | Total energy |
| ke | Kinetic energy |
| peBond | Bond potential energy |
| peAngle | Angle potential energy |
| temp | Temperature |
| press | Pressure |
| stressX/Y/Z | Total virial stress (all atoms) |
| stressX_N/Y_N/Z_N | Network-only virial stress (NC1 + 7S beads, excluding free monomers) |
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

The primary observables are `stressX_N` and `stressY_N` — the network contribution to the virial stress normalised by volume. Their time evolution during the relaxation phase gives the stress relaxation curve.
