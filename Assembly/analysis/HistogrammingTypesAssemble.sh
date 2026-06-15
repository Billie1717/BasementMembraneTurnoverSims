#!/bin/bash
# HistogrammingTypesAssemble.sh
# Post-processing script for assembly simulations.
#
# For each run, calls Histogram_Types.py to count the number of beads
# of each type (i.e. bond valency state) at each saved trajectory frame.
# This gives a time series showing the progress of network self-assembly.
#
# Run this from the Assembly/ directory after simulations are complete.
# Output .txt files are written to analysis/Data/ for use by the
# plotting notebooks.
#
# Usage: bash analysis/HistogrammingTypesAssemble.sh

# ── Parameters (must match those used in Assemble.sh) ────────────────────────
Time_Eq=2e7
NumFrames=100
Nevery=500
tstep=0.01
MD=0.65
MD7s=0.65
KangNC1=4.0
FactorMult=0.66

mkdir -p analysis/Data

for seed in 1 2 3
do
for BD in 1.35
do
for BD7s in 0.95
do
for MP in 0.12 0.18 0.27
do

MP7s=${MP}

filepattern='runs/runDense_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_KangNC1'${KangNC1}'_MP'${MP}'_MP'${MP7s}'_FactMu'${FactorMult}'_MD'${MD}'_BD'${BD}'_MD'${MD7s}'_BD'${BD7s}'_seed'${seed}'/'

fileOut='analysis/Data/HistAssemble_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_KangNC1'${KangNC1}'_MP'${MP}'_MP'${MP7s}'_FactMu'${FactorMult}'_MD'${MD}'_BD'${BD}'_MD'${MD7s}'_BD'${BD7s}'_seed'${seed}'.txt'

echo "Processing: ${filepattern}"
python3 analysis/Histogram_Types.py ${filepattern} ${fileOut}

done
done
done
done
