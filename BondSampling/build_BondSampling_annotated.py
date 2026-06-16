"""
build_BondSampling_annotated.py
================================
THIS SCRIPT IS NON-RUNNABLE. It is a heavily commented walkthrough of
build_BondSampling.py, focused on what differs from build_Equilibration.py.
For the actual simulation script, see build_BondSampling.py.

Overview
--------
This stage continues from a restart file from the Equilibration stage.
The LAMMPS input deck (collagen.in) is nearly identical to the Equilibration
stage — same force field, same GCE reactions, same bond/react fix. The only
additions are a set of variables that accumulate cumulative counts of bond
formation and breaking events across all reaction types. These are written to
thermo.dat and give the molecular turnover rate of the network.

The run is longer (2×10^7 steps vs 3×10^6) and saves fewer frames (50 vs 100)
since the goal is statistics, not structural snapshots.
"""


# =============================================================================
# SECTION 1: PARAMETERS (differences from Equilibration)
# =============================================================================
# Timing differences:
#   tEq      = 2e7   (vs 3e6 in Equilibration — longer run for statistics)
#   Frames   = 50    (vs 100 — fewer trajectory frames needed)
#   Nevery   = 500   (same — bond reactions attempted every 500 steps)
#
# All other parameters (bond distances, probabilities, GCE, force field)
# are identical to Equilibration.


# =============================================================================
# SECTIONS 2–6: IDENTICAL TO EQUILIBRATION
# =============================================================================
"""
Initialisation (read_restart), force field, molecule templates, GCE
acceptance probability variables, and the bond/react fix are all identical
to build_Equilibration.py. See build_Equilibration_annotated.py for details.
"""


# =============================================================================
# SECTION 7: BOND COUNTING VARIABLES (new vs Equilibration)
# =============================================================================
"""
After the bond/react fix is defined, LAMMPS exposes per-reaction event
counters as f_freact_creation[N], where N is the reaction index (1-based,
in the order reactions appear in the fix command):

  [1]  Nucleation   (GCE insertion)
  [2]  Death        (GCE deletion)
  [3]  Rc1_2N_01   (NC1 bond formation, first variant)
  [4]  Rc1_2N_02
  ...
  [13] MonoToDimer_01  (7S bond formation, first variant)
  ...
  [30] BRc1_2N_01  (NC1 bond breaking, first variant)
  ...
  [40] BMonoToDimer_01 (7S bond breaking)
  ...

Each counter increments by 1 each time that reaction fires. To get the number
of bonds formed/broken (rather than reaction events), multi-bond reactions
(e.g. the first NC1 or 7S bond formed from a free monomer, which creates 2
bonds simultaneously) are weighted by 2; all others by 1.

Cumulative release (breaking) variables:
  NumRelNC1  = 2*[3] + [4] + [5] + ... + [12]   # NC1 breaking events
  NumRel7S_1 = 2*[13] + [14] + ... + [19]        # 7S Monomer→Dimer breaking
  NumRel7S_2 = [20] + [21] + ... + [29]          # 7S higher-order breaking
  NumRelease = NumRelNC1 + NumRel7S_1 + NumRel7S_2

Cumulative incorporation (formation) variables:
  NumIncNC1  = 2*[30] + [31] + ... + [39]        # NC1 formation events
  NumInc7S_1 = 2*[40] + [41] + ... + [46]        # 7S Monomer→Dimer formation
  NumInc7S_2 = [47] + ... + [56]                 # 7S higher-order formation
  NumIncorp  = NumIncNC1 + NumInc7S_1 + NumInc7S_2

These variables are cumulative — they grow monotonically over the run.
The bond turnover rate is extracted from the slope of NumRelease or NumIncorp
versus time (they should be equal at steady state).
"""


# =============================================================================
# SECTION 8: OUTPUT (additional thermo columns vs Equilibration)
# =============================================================================
"""
thermo.dat has all Equilibration columns plus two more at the end:

  ... T_id_add T_id_remov NumTot NumRelease NumIncorp

NumRelease  -- cumulative bond-breaking events (weighted, all types)
NumIncorp   -- cumulative bond-forming events (weighted, all types)
"""


# =============================================================================
# SECTION 9: RUN AND OUTPUT
# =============================================================================
"""
run         ${tbonds}           # 2e7 steps
write_data  data.postmakebonds
write_restart restart.postmakebonds
"""


# =============================================================================
# SECTION 10: SLURM SUBMISSION SCRIPT
# =============================================================================
"""
Identical to Equilibration. Replace module load commands and LAMMPS path
with your local equivalents.
"""
