#!/bin/bash
# Growth.sh
# Launcher for the Growth stage: constant-stress (press/berendsen) NVT runs
# starting from Equilibration restart files. The box expands or contracts
# freely in x and y while the stress is held at a target value, giving
# indefinite growth (or shrinkage) of the network area.
#
# THREE EXPERIMENTS are run from this launcher by switching the active
# parameter block and build script (see comments below):
#
#   Experiment A — Vary MP at fixed stress
#     Build script : build_Growth.py (GCE + bond breaking active)
#     MPs          : 0.0001 0.0002 0.0004 0.0008
#     Press        : -0.1
#     N_preferred  : 5 (all MPs)
#
#   Experiment B — Vary stress at fixed MP
#     Build script : build_Growth.py (GCE + bond breaking active)
#     MPs          : 0.0004
#     Press        : -0.1 -0.11 -0.12 -0.13 -0.15
#     N_preferred  : 5
#
#   Experiment C — Three regimes at fixed stress and fixed MP
#     Regime i   (no bond breakage):
#       Build script : build_Growth_NoGCE.py, MP=0.0  (BreakProb=0, GCE off)
#     Regime ii  (bond breakage, no new material):
#       Build script : build_Growth_NoGCE.py, MP=0.01 (breaking active, GCE off)
#     Regime iii (bond breakage + new material):
#       Build script : build_Growth.py,       MP=0.01 (breaking active, GCE on)
#     Press        : -0.1
#     N_preferred  : 5
#
# To switch experiments: comment/uncomment the MPs, NumMons, Press loops,
# and python3 call in the parameter section below.
#
# Prerequisites: Equilibration stage must be complete. The Equilibration
# runs/ directory is expected at ../Equilibration/runs/ relative to this script.
#
# Usage: bash Growth.sh
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
Time_Eq=1e7       # NVT steps for the growth run
Time_Str=0
Time_Rel=0
NumFrames=100     # Number of trajectory frames to save
Nevery=500        # Timesteps between bond reaction attempts
tstep=0.01        # Integration timestep

# ── Bond parameters ───────────────────────────────────────────────────────────
MD=0.65
BD=1.35
MD7s=0.65
BD7s=0.95
KangNC1=4.0
FactorMult=0.66
OvaR=3.5
Xstretch=51       # Passed to build script (unused in collagen.in)

# ── GCE parameters ────────────────────────────────────────────────────────────
kappa_den=0.01
InitialVol=3600   # Box volume for this simulation density

# ── Equilibration runs directory ─────────────────────────────────────────────
equib_runs="../Equilibration/runs"

# ═════════════════════════════════════════════════════════════════════════════
# EXPERIMENT A: Vary MP at fixed stress
# Uncomment this block and comment out Experiments B and C
# ═════════════════════════════════════════════════════════════════════════════
MPs=(0.0001 0.0002 0.0004 0.0008)
NumMons=(5 5 5 5)
BUILD_SCRIPT="build_Growth.py"
RUN_PREFIX="runGrowth_VaryMP"

# ═════════════════════════════════════════════════════════════════════════════
# EXPERIMENT B: Vary stress at fixed MP
# Comment out Experiment A above and uncomment this block
# ═════════════════════════════════════════════════════════════════════════════
# MPs=(0.0004)
# NumMons=(5)
# BUILD_SCRIPT="build_Growth.py"
# RUN_PREFIX="runGrowth_VaryPress"

# ═════════════════════════════════════════════════════════════════════════════
# EXPERIMENT C: Three regimes
# Comment out Experiments A and B above and uncomment ONE of these blocks.
#
# Regime i — no bond breakage (MP=0.0 → BreakProb=0, GCE off):
# ═════════════════════════════════════════════════════════════════════════════
# MPs=(0.0)
# NumMons=(5)
# BUILD_SCRIPT="build_Growth_NoGCE.py"
# RUN_PREFIX="runGrowth_NoGCE_NoBreak"
#
# Regime ii — bond breakage, no new material (GCE off):
# MPs=(0.01)
# NumMons=(5)
# BUILD_SCRIPT="build_Growth_NoGCE.py"
# RUN_PREFIX="runGrowth_NoGCE_Break"
#
# Regime iii — bond breakage + new material (GCE on):
# MPs=(0.01)
# NumMons=(5)
# BUILD_SCRIPT="build_Growth.py"
# RUN_PREFIX="runGrowth_BreakAndAdd"

# ── Parameter loops ───────────────────────────────────────────────────────────
for seed in 1 2 3
do
for BD in 1.35
do
for BD7s in 0.95
do
for Press in -0.1        # Experiment A: use -0.1
                         # Experiment B: use -0.1 -0.11 -0.12 -0.13 -0.15
                         # Experiment C: use -0.1
do
for Pdamp in 100.0
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

foldernameadd='runs/'${RUN_PREFIX}'_kappaD'${kappa_den}'_teq'${Time_Eq}'_Frames'${NumFrames}'_tstep'${tstep}'_Nev'${Nevery}'_KangNC1'${KangNC1}'_MP'${MP}'_MP'${MP7s}'_FactMu'${FactorMult}'_MonsPref'${N_preferred}'_MD'${MD}'_BD'${BD}'_MD'${MD7s}'_BD'${BD7s}'_Press'${Press}'_Pdamp'${Pdamp}'_seed'${seed}

mkdir -p ${foldernameadd}
cp ${restart_path} ${foldernameadd}'/data'

python3 ${BUILD_SCRIPT} ${foldernameadd} ${Time_Eq} ${Time_Str} ${Time_Rel} ${NumFrames} ${Nevery} ${tstep} ${MD} ${BD} ${MD7s} ${BD7s} ${MP} ${MP7s} ${FactorMult} ${KangNC1} ${N_preferred} ${Xstretch} ${OvaR} ${kappa_den} ${InitialVol} ${Press} ${Pdamp} ${seed}

cd ${foldernameadd}
sbatch runscript.sh
cd ../..

done
done
done
done
done
done
done
done
