#!/bin/bash
# Assemble.sh
# Launcher for collagen IV network self-assembly simulations.
#
# For each combination of parameters, this script:
#   1. Creates an output directory
#   2. Copies the initial configuration file
#   3. Calls build_Assembly.py to write the LAMMPS input deck (collagen.in)
#      and a SLURM submission script (runscript.sh)
#   4. Submits the job via sbatch
#
# Parameters used for the paper are set as defaults below.
# Adjust the loop ranges to run parameter sweeps.
#
# Usage: bash Assemble.sh
# Requires: Python 3, LAMMPS built with custom REACTION package (see README)

# ── Simulation timing ────────────────────────────────────────────────────────
Time_Eq=2e7       # NVT steps for bond equilibration
Time_Str=0        # NVT steps for stretch (0 = no stretch in this stage)
Time_Rel=0        # NVT steps for relaxation
NumFrames=100     # Number of trajectory frames to save
Nevery=500        # How often bond reactions are attempted (timesteps)
tstep=0.01        # Integration timestep

# ── Bond parameters ───────────────────────────────────────────────────────────
# NC1 domain bonds
MD=0.65           # Make distance for NC1 bonds
# BD set in loop below

# 7S domain bonds
MD7s=0.65         # Make distance for 7S bonds
# BD7s set in loop below

# ── Other parameters ──────────────────────────────────────────────────────────
KangNC1=4.0       # Angular stiffness of NC1 bonds
FactorMult=0.66   # Break probability = FactorMult * Make probability
OvaR=3.5          # Overlap radius for monomer insertion (GCE)
N_preferred=0     # Preferred number of monomers (0 = no preference)
Xstretch=51       # Box length in x (used if stretch stage is enabled)

# ── Initial configuration ─────────────────────────────────────────────────────
input='input/dataNucType9_dense'

# ── Parameter loops ───────────────────────────────────────────────────────────
for seed in 1 #2 3
do
for BD in 1.35
do
for BD7s in 0.95
do
for MP in 0.12 #0.18 0.27
do

MP7s=${MP}

foldernameadd='runs/runDense_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_KangNC1'${KangNC1}'_MP'${MP}'_MP'${MP7s}'_FactMu'${FactorMult}'_MD'${MD}'_BD'${BD}'_MD'${MD7s}'_BD'${BD7s}'_seed'${seed}

mkdir -p ${foldernameadd}
cp ${input} ${foldernameadd}'/data'

echo "Submitting: ${foldernameadd}"

python3 build_Assembly.py ${foldernameadd} ${Time_Eq} ${Time_Str} ${Time_Rel} \
    ${NumFrames} ${Nevery} ${tstep} ${MD} ${BD} ${MD7s} ${BD7s} \
    ${MP} ${MP7s} ${FactorMult} ${KangNC1} ${N_preferred} ${Xstretch} ${OvaR} ${seed}

cd ${foldernameadd}
sbatch runscript.sh
cd ../..

done
done
done
done
