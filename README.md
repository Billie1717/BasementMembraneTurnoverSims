# Basement membrane turnover controls tissue shape — simulation code

Code to reproduce the coarse-grained MD simulations in:

> Barrientos et al. (2026) "Basement membrane turnover controls tissue shape." *Cell.*

---

## Model overview

Each collagen IV protomer is represented as a two-bead rod: one NC1 bead and one 7S bead. NC1 beads bind pairwise (end-to-end); 7S beads bind laterally to up to three partners. Bond formation and breaking are implemented via the LAMMPS `bond/react` package, which requires one pre/post/map template triple per allowed change in NC1 or 7S valency — hence the large `templates/` directories. A Grand Canonical Ensemble (GCE) controls free monomer exchange via stochastic insertion (Nucleation) and deletion (Death) reactions.

## Requirements

**LAMMPS** `lammps-15Jun2023` built with the custom `REACTION` package. Replace `src/REACTION/fix_bond_react.{cpp,h}` in the LAMMPS source tree with the versions in `Lammps/` before building.

**Python ≥ 3.9** with `numpy`.

## Simulation stages

Each stage folder contains a launcher `*.sh`, a build script `build_*.py` that writes the LAMMPS input (`collagen.in`) and SLURM submission script (`runscript.sh`), a non-runnable annotated walkthrough (`build_*_annotated.py`), and a `parameters.md` with the exact values used in the paper.

| Folder | Description |
|--------|-------------|
| `Assembly/` | Self-assembly of the network from a lattice initial configuration |
| `Equilibration/` | GCE equilibration — stochastic monomer insertion/deletion to reach target concentration |
| `BondSampling/` | Measures molecular turnover rate (cumulative bond formation and breaking events) |
| `Stretch/` | Biaxial strain applied in 20 steps, then stress relaxation at fixed strain |
| `Growth/` | Constant tensile stress applied via `press/berendsen`; network area grows freely |

Stages run in order: Assembly → Equilibration → {BondSampling, Stretch, Growth}.

## Running

Each stage follows the same pattern:

```bash
cd <Stage>/
bash <Stage>.sh        # builds collagen.in + runscript.sh and submits via sbatch
```

To run without SLURM:

```bash
cd <run_folder>/
mpirun -np 1 /PATH/lammps-15Jun2023/src/lmp_mpi -in collagen.in
```

**For questions/comments please contact:** Billie Meadowcroft (billiemead@hotmail.co.uk)  


See each stage's `README.md` for prerequisites and output file descriptions.
