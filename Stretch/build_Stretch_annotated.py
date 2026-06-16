"""
build_Stretch_annotated.py
==========================
THIS SCRIPT IS NON-RUNNABLE. It is a heavily commented walkthrough of
build_Stretch.py, focused on what differs from build_Equilibration.py.
For the actual simulation script, see build_Stretch.py.

Overview
--------
This stage continues from a restart file from the Equilibration stage.
The force field, GCE, and bond/react fix are identical to Equilibration.
The key difference is the run protocol: instead of a single long NVT run,
the simulation proceeds in three phases:

  1. Stretch phase    (tStr = 7.5e4 steps):
       Box is deformed biaxially in 20 incremental steps.
  2. Short relaxation (tInitRel = tStr/2 steps):
       Bond reactions re-initialised; short NVT at fixed strain.
  3. Long relaxation  (tRel = 1e6 steps):
       NVT at fixed strain — this is the main measurement window.

The readout is the network stress time series in thermo.dat.
"""


# =============================================================================
# SECTION 1: PARAMETERS (differences from Equilibration)
# =============================================================================
# Timing:
#   Time_Eq   = 0        (no pre-stretch NVT — system already equilibrated)
#   Time_Str  = 7.5e4    (stretch phase, divided into 20 steps)
#   Time_Rel  = 1e6      (relaxation phase at fixed strain)
#   NumFrames = 100
#   Nevery    = 5000     (vs 500 in Equilibration — reactions less frequent
#                         during the fast stretch, matches cluster runs)
#
# Strain:
#   Xstretch  = 50       (total box length increase per direction, LJ units)
#   Nstretch  = 20       (number of incremental steps; hardcoded in build script)
#
#   Each step expands one dimension by +/-(Xstretch/2)/Nstretch = ±1.25 on
#   each side (so +2.5 per step), alternating x then y. Total expansion:
#   50 LJ units in both x and y.


# =============================================================================
# SECTIONS 2–6: IDENTICAL TO EQUILIBRATION
# =============================================================================
"""
Initialisation (read_restart), force field, molecule templates, GCE
acceptance probability variables, and the bond/react fix are all identical
to build_Equilibration.py. See build_Equilibration_annotated.py for details.
"""


# =============================================================================
# SECTION 7: STRETCH PHASE
# =============================================================================
"""
After the initial bond/react fix is defined, the stretch is applied as a
Python loop that writes 20 pairs of LAMMPS commands:

  for i in range(20):
      fix fDeformX all deform 1 x delta -1.25 +1.25 remap x
      run  ${tstretchcycl}        # = tStr / 20 = 3750 steps
      unfix fDeformX

      fix fDeformY all deform 1 y delta -1.25 +1.25 remap x
      run  ${tstretchcycl}
      unfix fDeformY

  write_data   data.poststretch
  write_restart restart.poststretch

`fix deform x delta lo hi` shifts the lower boundary by `lo` and the upper
by `hi` over the duration of the run. The sign convention means each step
shrinks the low side by 1.25 and expands the high side by 1.25, for a net
box length increase of 2.5 per step. After 20 steps: +50 total in x and y.

`remap x` keeps atom fractional coordinates constant during the deformation
(atoms move affinely with the box).

After the 20 steps the box dimensions are frozen — no further deform fixes.
"""


# =============================================================================
# SECTION 8: SHORT RELAXATION AND BOND REACTION RE-INITIALISATION
# =============================================================================
"""
After the stretch phase:
  1. `run ${tInitRel}` (= tStr/2 = 37500 steps) — short NVT at the new strain.
     During this run the original bond/react fix (defined before the stretch)
     is still active.
  2. `unfix freact_creation` — the bond/react fix is removed.
  3. A new identical bond/react fix is defined with the same reactions.
     This reset is needed because internal LAMMPS counters from the stretch
     phase can cause issues; re-defining the fix clears them cleanly.

The two-step reset (run briefly, unfix, refix) ensures bond kinetics
re-equilibrate at the new box geometry before the main relaxation.
"""


# =============================================================================
# SECTION 9: LONG RELAXATION PHASE (main measurement window)
# =============================================================================
"""
  run  ${trelax}          # = 1e6 steps

  write_data   data.postrelax
  write_restart restart.postrelax

This is the primary measurement window. The box is held at the stretched
dimensions (no deform fix) and the network stress relaxes as bonds rearrange.
The stress time series is read from thermo.dat columns stressX_N / stressY_N
(network contribution to the virial stress, normalised by volume).
"""


# =============================================================================
# SECTION 10: THERMO OUTPUT (differences from Equilibration)
# =============================================================================
"""
thermo.dat has a different column set from Equilibration — no NumRelease /
NumIncorp (those are only in BondSampling), but adds stress columns:

  step etot ke peBond peAngle temp press
  stressX stressY stressZ         # total virial stress (all atoms)
  stressX_N stressY_N stressZ_N  # network-only virial stress
  AvBondForce
  NumMons BoxLength BoxVol CurrentConc ConcPref
  P_add P_remov NreactNuc NreactDeath
  T_den_add T_den_remov T_id_add T_id_remov
  NumTot

stressX_N / stressY_N are computed from the NC1+7S network beads only
(excluding free type-9 monomers), giving the elastic network stress.
"""


# =============================================================================
# SECTION 11: SLURM SUBMISSION SCRIPT
# =============================================================================
"""
Identical to Equilibration. Replace module load commands and LAMMPS path
with your local equivalents.
"""
