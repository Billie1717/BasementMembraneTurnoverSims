#!/bin/bash
# HistogrammingTypesBondSampling.sh
# Post-processing script for bond sampling simulations.
#
# For each run, calls Histogram_Types.py to count the number of beads
# of each type (bond valency state) at each saved trajectory frame.
#
# Run this from the BondSampling/ directory after simulations are complete.
# Output .txt files are written to analysis/Data/.
#
# Usage: bash analysis/HistogrammingTypesBondSampling.sh

# ── Parameters (must match those used in BondSampling.sh) ────────────────────
Time_Eq=2e7
NumFrames=50
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

filepattern='runs/runBondSample_kappaD'${kappa_den}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_KangNC1'${KangNC1}'_MP'${MP}'_MP'${MP7s}'_FactMu'${FactorMult}'_MonsPref'${N_preferred}'_MD'${MD}'_BD'${BD}'_MD'${MD7s}'_BD'${BD7s}'_seed'${seed}'/'

fileOut='analysis/Data/HistBondSample_kappaD'${kappa_den}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_KangNC1'${KangNC1}'_MP'${MP}'_MP'${MP7s}'_FactMu'${FactorMult}'_MonsPref'${N_preferred}'_MD'${MD}'_BD'${BD}'_MD'${MD7s}'_BD'${BD7s}'_seed'${seed}'.txt'

echo "Processing: ${filepattern}"
python3 analysis/Histogram_Types.py ${filepattern} ${fileOut}

done
done
done
done
