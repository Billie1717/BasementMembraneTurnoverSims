#!/bin/bash
# BondSampling.sh
# Launcher for the Bond Sampling stage.
#
# Reads a restart file from the Equilibration stage and runs a long NVT
# simulation that counts cumulative bond formation and breaking events
# (incorporation and release) for all NC1 and 7S bond types. This measures
# molecular turnover in the assembled network and can be applied to networks
# from any stage (equilibrated, stretched, or grown).
#
# Prerequisites: Equilibration stage must be complete. The Equilibration
# runs/ directory is expected at ../Equilibration/runs/ relative to this script.
#
# Usage: bash BondSampling.sh
# Requires: Python 3, LAMMPS lammps-15Jun2023 with custom REACTION package

# ── Helper: find the latest restart file in a directory ──────────────────────
get_latest_restart() {
  local d="$1"
  local by_name=""
  by_name=$(
    shopt -s nullglob
    for f in "$d"/run_T*.restart; do
      local base=${f##*/}
      if [[ $base =~ ^run_T([0-9]+)(\.([0-9]+))?\.restart$ ]]; then
        local int=${BASH_REMATCH[1]}
        local frac=${BASH_REMATCH[3]}
        local key="${int}${frac}"
        printf '%s %s\n' "$key" "$f"
      fi
    done | LC_ALL=C sort -k1,1n | awk '{print $2}' | tail -n 1
  )
  if [[ -n "${by_name}" ]]; then
    printf '%s\n' "${by_name}"
    return 0
  fi
  # Fallback: newest by modification time
  local by_mtime=""
  if command -v find >/dev/null 2>&1; then
    by_mtime=$(find "$d" -maxdepth 1 -type f -name 'run_T*.restart' -printf '%T@ %p\n' 2>/dev/null \
      | LC_ALL=C sort -nr \
      | awk '{ $1=""; sub(/^ /,""); print }' \
      | head -n 1 || true)
  fi
  if [[ -n "${by_mtime}" ]]; then
    printf '%s\n' "${by_mtime}"
    return 0
  fi
  return 1
}

# ── Simulation timing ─────────────────────────────────────────────────────────
Time_Eq=2e7       # NVT steps for bond sampling run
Time_Str=0
Time_Rel=0
NumFrames=50      # Number of trajectory frames to save
Nevery=500        # How often bond reactions are attempted (timesteps)
tstep=0.01        # Integration timestep

# ── Bond parameters ───────────────────────────────────────────────────────────
MD=0.65           # NC1 bond formation distance cutoff
BD=1.35           # NC1 bond breaking distance cutoff
MD7s=0.65         # 7S bond formation distance cutoff
BD7s=0.95         # 7S bond breaking distance cutoff
KangNC1=4.0       # NC1 angular stiffness
FactorMult=0.66   # BreakProb = FactorMult * MakeProb
OvaR=3.5          # Passed to build script (unused in collagen.in — see annotated script)

# ── GCE parameters (carried over from Equilibration; GCE remains active) ─────
kappa_den=0.01
InitialVol=12500

# ── Paper MP values and paired N_preferred ────────────────────────────────────
MPs=(0.12 0.18 0.27)
NumMons=(8 8 8)

# ── Equilibration runs directory (relative to this script) ───────────────────
equib_runs="../Equilibration/runs"

# ── Parameter loops ───────────────────────────────────────────────────────────
for seed in 1 2 3
do
for BD in 1.35
do
for BD7s in 0.95
do
for kappa_den in 0.01
do
for OvaR in 3.5
do
for i in "${!MPs[@]}"; do
MP="${MPs[$i]}"
N_preferred="${NumMons[$i]}"
MP7s=${MP}

# Locate the Equilibration restart file for this parameter combination
base="${equib_runs}/runEquib_kappaD${kappa_den}_teq3e6_Frames100_tstep${tstep}_Nev500_KangNC1${KangNC1}_MP${MP}_MP${MP7s}_FactMu${FactorMult}_MonsPref${N_preferred}_MD${MD}_BD${BD}_MD${MD7s}_BD${BD7s}_seed${seed}"
restart_dir="${base}/restart"
if ! restart_path="$(get_latest_restart "${restart_dir}")"; then
  echo "ERROR: No restart files found in ${restart_dir}" >&2
  continue
fi
echo "Using restart: ${restart_path}"

foldernameadd='runs/runBondSample_kappaD'${kappa_den}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_KangNC1'${KangNC1}'_MP'${MP}'_MP'${MP7s}'_FactMu'${FactorMult}'_MonsPref'${N_preferred}'_MD'${MD}'_BD'${BD}'_MD'${MD7s}'_BD'${BD7s}'_seed'${seed}

mkdir -p ${foldernameadd}
cp ${restart_path} ${foldernameadd}'/data'

python3 build_BondSampling.py ${foldernameadd} ${Time_Eq} ${Time_Str} ${Time_Rel} ${NumFrames} ${Nevery} ${tstep} ${MD} ${BD} ${MD7s} ${BD7s} ${MP} ${MP7s} ${FactorMult} ${KangNC1} ${N_preferred} 51 ${OvaR} ${kappa_den} ${InitialVol} ${seed}

cd ${foldernameadd}
sbatch runscript.sh
cd ../..

done
done
done
done
done
done
