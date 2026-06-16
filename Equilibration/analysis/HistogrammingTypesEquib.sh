#!/bin/bash
# HistogrammingTypesEquib.sh
# Post-processing script for equilibration simulations.
#
# For each run, calls Histogram_Types.py to count the number of beads
# of each type (bond valency state) at each saved trajectory frame.
# This gives a time series used to verify that the GCE has driven the
# network to the target monomer count and that bond topology has equilibrated.
#
# Run this from the Equilibration/ directory after simulations are complete.
# Output .txt files are written to analysis/Data/.
#
# Usage: bash analysis/HistogrammingTypesEquib.sh

# ── Parameters (must match those used in Equilibrate.sh) ─────────────────────
Time_Eq=3e6
NumFrames=100
Nevery=500
tstep=0.01
MD=0.65
MD7s=0.65
KangNC1=4.0
FactorMult=0.66
kappa_den=0.01

mkdir -p analysis/Data

for seed in 1 2 3
do
for BD in 1.35
do
for BD7s in 0.95
do
for i in 0 1 2
do

MPs=(0.12 0.18 0.27)
NumMons=(8 8 8)
MP="${MPs[$i]}"
N_preferred="${NumMons[$i]}"
MP7s=${MP}

filepattern='runs/runEquib_kappaD'${kappa_den}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_KangNC1'${KangNC1}'_MP'${MP}'_MP'${MP7s}'_FactMu'${FactorMult}'_MonsPref'${N_preferred}'_MD'${MD}'_BD'${BD}'_MD'${MD7s}'_BD'${BD7s}'_seed'${seed}'/'

fileOut='analysis/Data/HistEquib_kappaD'${kappa_den}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_KangNC1'${KangNC1}'_MP'${MP}'_MP'${MP7s}'_FactMu'${FactorMult}'_MonsPref'${N_preferred}'_MD'${MD}'_BD'${BD}'_MD'${MD7s}'_BD'${BD7s}'_seed'${seed}'.txt'

echo "Processing: ${filepattern}"
python3 analysis/Histogram_Types.py ${filepattern} ${fileOut}

done
done
done
done
