"""
build_Growth_annotated.py
=========================
THIS SCRIPT IS NON-RUNNABLE. It is a heavily commented walkthrough of
build_Growth.py (and build_Growth_NoGCE.py), focused on what differs from
build_Equilibration.py. For the actual simulation scripts, see those files.

Overview
--------
This stage continues from a restart file from the Equilibration stage.
The force field and bond/react reactions are identical to Equilibration,
with two key changes:

  1. STRESS CONTROL: instead of a fixed box, `fix press/berendsen` is
     applied in x and y, coupling the box dimensions to a target stress.
     The box area grows or shrinks freely to maintain the target.

  2. BOND COUNTING: as in BondSampling, cumulative NumRelease and NumIncorp
     variables are tracked (in build_Growth.py only; absent in
     build_Growth_NoGCE.py).

Two build scripts exist for the three-regime experiment:

  build_Growth.py       — GCE (Nucleation + Death) active. Used for
                          Experiments A, B, and regime iii of Experiment C.
  build_Growth_NoGCE.py — GCE disabled (Nucleation and Death reactions
                          commented out). Used for regimes i and ii of
                          Experiment C.
"""


# =============================================================================
# SECTION 1: PARAMETERS (differences from Equilibration)
# =============================================================================
# New parameters vs Equilibration:
#   Press     — target stress applied in x and y (negative = tensile)
#   Pdamp     — Barendsen damping timescale (100.0)
#
# Timing:
#   Time_Eq   = 1e7    (longer than Equilibration's 3e6)
#   NumFrames = 100
#   Nevery    = 500
#
# Note: InitialVol = 3600 here (different simulation density from the
# Assembly/Equilibration/BondSampling/Stretch stages which use 12500).


# =============================================================================
# SECTIONS 2–6: IDENTICAL TO EQUILIBRATION (build_Growth.py)
#               OR IDENTICAL BUT WITHOUT NUCLEATION/DEATH (build_Growth_NoGCE.py)
# =============================================================================
"""
Initialisation (read_restart), force field, molecule templates, GCE
acceptance probability variables are defined in the same way as
build_Equilibration.py.

In build_Growth_NoGCE.py the Nucleation and Death reactions are commented
out of the bond/react fix, so no monomers are inserted or deleted.
The GCE probability variables (P_add, P_remov etc.) are still written and
printed to thermo.dat but are never actually used by any reaction.
"""


# =============================================================================
# SECTION 7: THREE-REGIME LOGIC (Experiment C only)
# =============================================================================
"""
The three regimes are produced purely by parameter choice — no code
branching is needed:

  Regime i  — No bond breakage:
      MP = 0.0  →  BreakProb = FactorMult * MP = 0.0
      All breaking reactions fire with prob = 0, so no bonds ever break.
      build_Growth_NoGCE.py used (GCE also off).

  Regime ii — Bond breakage, no new material:
      MP = 0.01 →  BreakProb = FactorMult * 0.01 = 0.0066
      Breaking reactions active; Nucleation/Death commented out.
      build_Growth_NoGCE.py used.

  Regime iii — Bond breakage + new material (full model):
      MP = 0.01 →  BreakProb = 0.0066
      All reactions active including Nucleation and Death.
      build_Growth.py used.
"""


# =============================================================================
# SECTION 8: STRESS CONTROL FIX
# =============================================================================
"""
After the bond/react fix and thermo print fix are defined, two Berendsen
barostat fixes are applied:

  fix fDeformY all press/berendsen y ${PressTarget} ${PressTarget} ${Pdamp}
  fix fDeformX all press/berendsen x ${PressTarget} ${PressTarget} ${Pdamp}

`press/berendsen` rescales the box dimension each step to drive the
instantaneous pressure toward `PressTarget`. With a negative (tensile)
target, the box expands in x and y, stretching the network and causing
growth. The z dimension is not coupled and remains fixed.

`PressTarget` is written as a LAMMPS variable from the `Press` argument.
`Pdamp` controls how quickly the box responds to stress deviations.

The stress measured in thermo.dat columns `sigma_xx` / `sigma_yy` gives
the actual network stress at each timestep; their convergence toward
`PressTarget` confirms the barostat is working.
"""


# =============================================================================
# SECTION 9: BOND COUNTING (build_Growth.py only)
# =============================================================================
"""
Identical to BondSampling: cumulative NumRelease and NumIncorp are built
from f_freact_creation indices and appended to thermo.dat.

These are absent in build_Growth_NoGCE.py (the reaction list is shorter
so the index mapping would differ, and turnover is not the focus of the
no-GCE regimes).
"""


# =============================================================================
# SECTION 10: RUN AND OUTPUT
# =============================================================================
"""
run         ${tbonds}           # 1e7 steps
write_data  data.postmakebonds
write_restart restart.postmakebonds
"""


# =============================================================================
# SECTION 11: SLURM SUBMISSION SCRIPT
# =============================================================================
"""
Identical to Equilibration. Replace module load commands and LAMMPS path
with your local equivalents.
"""
