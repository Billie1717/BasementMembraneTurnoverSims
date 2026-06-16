"""
build_Equilibration_annotated.py
=================================
THIS SCRIPT IS NON-RUNNABLE. It is a heavily commented walkthrough of
build_Equilibration.py, intended to explain what the LAMMPS input deck
(collagen.in) does and why it differs from the Assembly stage.
For the actual simulation script, see build_Equilibration.py.

Overview
--------
This stage continues from a restart file produced by the Assembly stage.
The model and force field are identical to Assembly; the only change is that
the Grand Canonical Exchange (GCE) is now active: type-9 free collagen molecules are
inserted and deleted stochastically via `fix bond/react` Nucleation and Death
reactions. The acceptance probability follows a soft harmonic density penalty
(see parameters.md for the equations).

Sections below track collagen.in section by section, focusing on what is
new relative to build_Assembly_annotated.py.
"""


# =============================================================================
# SECTION 1: COMMAND-LINE PARAMETERS (new vs Assembly)
# =============================================================================
# Parameters identical to Assembly are not repeated here. New parameters:
#
#   PrefMonNum   -- N_preferred: target free monomer count for GCE
#   OvaR         -- passed to the script and written as a LAMMPS variable,
#                   but not referenced anywhere in collagen.in; effectively
#                   unused. The Nucleation reaction uses `modify_create nuc yes`
#                   which finds an overlap-free insertion site using the
#                   neighbour list cutoff, not OvaR.
#   kappa_density-- kappa_rho: strength of the harmonic density penalty
#   VolInit      -- InitialVol: reference box volume (fixes C_target)


# =============================================================================
# SECTION 2: LAMMPS INITIALISATION (difference from Assembly)
# =============================================================================
"""
The only difference here is that this stage starts from a RESTART FILE
rather than a fresh data file:

  read_restart ${fname}   # restart file from the last Assembly checkpoint

This reads the full system state (atom positions, velocities, bonds, angles)
so that the equilibration continues seamlessly from the assembled network.
The `reset_timestep 0` line present in Assembly is omitted here.
"""


# =============================================================================
# SECTION 3: FORCE FIELD
# =============================================================================
"""
Identical to Assembly (see build_Assembly_annotated.py, Section 3).
Bond coefficients, angle coefficients, and pair style are unchanged.
"""


# =============================================================================
# SECTION 4: MOLECULE TEMPLATES
# =============================================================================
"""
All Assembly templates are loaded (NC1_2N through NC1_8N, NC1_loop, 7S
topology changes). In addition, two NEW template pairs are active in this stage:

  molecule nucleation_pre  ../../templates/Nucleation/pre_Nucleation.txt
  molecule nucleation_post ../../templates/Nucleation/post_Nucleation.txt
  molecule death_pre       ../../templates/Nucleation/pre_Death.txt
  molecule death_post      ../../templates/Nucleation/post_Death.txt

  Nucleation: inserts a new type-9 free monomer at a randomly chosen empty site.
  Death: removes an existing type-9 free monomer.

These templates were present in the Assembly input but their reactions were not
registered in the bond/react fix. Here they are active.
"""


# =============================================================================
# SECTION 5: GCE ACCEPTANCE PROBABILITY VARIABLES
# =============================================================================
"""
This is the core addition relative to Assembly. After the computes and
integrator are set up, a set of LAMMPS variables implement the GCE acceptance
criterion.

Key variables:

  conc_preferred    = N_preferred / volInit     # target free-monomer concentration
  conc_pref_tot     = 3125 / volInit            # target TOTAL protomer concentration
                                                 # (fixed by the starting box volume)

  NumType9          = count of current type-9 atoms (free monomers)
  NumType1          = count of current type-1 atoms (bound NC1)
  NumTot            = NumType1 + NumType9        # total protomers

  CurrentConc       = NumType9 / box_Vol         # current free-monomer concentration
  CurrentConc_a     = (NumType9+1) / box_Vol     # concentration after a hypothetical insertion

Energy penalty (phi) before and after insertion:
  En_o_add = KappaD * (NumTot       - conc_pref_tot * box_Vol)^2
  En_n_add = KappaD * (NumTot + 1   - conc_pref_tot * box_Vol)^2

Energy penalty before and after deletion:
  En_o_remov = KappaD * (NumTot     - conc_pref_tot * box_Vol)^2
  En_n_remov = KappaD * (NumTot - 1 - conc_pref_tot * box_Vol)^2

Boltzmann factors:
  TermDen_add   = exp(-(En_n_add   - En_o_add))
  TermDen_remov = exp(-(En_n_remov - En_o_remov))

Ideal-gas terms (ratio of concentrations / particle counts):
  TermId_add   = conc_preferred / (CurrentConc_a * 2)
  TermId_remov = CurrentConc / (conc_preferred * NumType9)

Final acceptance probabilities (scaled by ConstRate = 1/100 for numerical stability):
  P_add   = ConstRate * TermDen_add   * TermId_add
  P_remov = ConstRate * TermDen_remov * TermId_remov

These are per-atom variables passed to the Nucleation and Death reactions as
the stochastic acceptance probability.
"""


# =============================================================================
# SECTION 6: BOND/REACT FIX (GCE + ASSEMBLY REACTIONS)
# =============================================================================
"""
The fix bond/react command now registers Nucleation and Death as the FIRST
two reactions, followed by all Assembly reactions unchanged:

  fix freact_creation all bond/react stabilization no reset_mol_ids no &
      react Nucleation all ${NeveryFast} ${RminF} ${RmaxFB} \
            nucleation_pre nucleation_post .../map_Nucleation.txt \
            prob 1 ${lseed1} modify_create nuc yes &
      react Death all ${NeveryFast} ${RminF} ${RmaxFB} \
            death_pre death_post .../map_Death.txt \
            prob 1 ${lseed1} &
      react Rc1_2N_01 ... prob v_probF ... &
      react BRc1_2N_01 ... prob v_probFB ... &
      ... (all Assembly reactions as before)

Key differences from Assembly:
  - Nucleation uses `modify_create nuc yes` to place the new monomer at an
    overlap-free site (determined by the neighbour list cutoff; OvaR is passed
    to the script but not actually used in collagen.in).
  - Both Nucleation and Death use NeveryFast (= 1): they are attempted every
    timestep, but acceptance is gated by P_add / P_remov.
  - The per-timestep acceptance (prob 1 for the react command, then filtered
    by the per-atom variable P_add/P_remov) means the effective rate is
    controlled entirely by the GCE variables above.

After the fix is defined, reaction counters are collected:
  variable NumNuc   equal f_freact_creation[1]  # cumulative insertions
  variable NumDeath equal f_freact_creation[2]  # cumulative deletions
"""


# =============================================================================
# SECTION 7: OUTPUT (additional thermo columns)
# =============================================================================
"""
The thermo.dat file has the same columns as Assembly plus GCE diagnostics:

  step etot ke peBond peAngle temp press stressX stressY stressZ
  AvBondForce NumMons BoxLength BoxVol CurrentConc ConcPref
  P_add P_remov NreactNuc NreactDeath T_den_add T_den_remov T_id_add T_id_remov NumTot

New columns:
  P_add        -- insertion acceptance probability at this step
  P_remov      -- deletion acceptance probability at this step
  NreactNuc    -- cumulative number of successful insertions
  NreactDeath  -- cumulative number of successful deletions
  T_den_add    -- Boltzmann factor for insertion (exp term)
  T_den_remov  -- Boltzmann factor for deletion
  T_id_add     -- ideal-gas factor for insertion
  T_id_remov   -- ideal-gas factor for deletion
  NumTot       -- total protomer count (bound + free)
"""


# =============================================================================
# SECTION 8: RUN AND OUTPUT
# =============================================================================
"""
run         ${tbonds}           # Run for Time_Eq timesteps with active GCE
write_data  data.postmakebonds  # Save final configuration
write_restart restart.postmakebonds
"""


# =============================================================================
# SECTION 9: SLURM SUBMISSION SCRIPT (runscript.sh)
# =============================================================================
"""
Identical to Assembly. The Python script writes runscript.sh which:
  1. Creates output subdirectories: dumplin/, dumplin_bonds/, restart/
  2. Generates random seeds for the Langevin thermostat
  3. Runs LAMMPS in serial (mpirun -np 1) on a single node

To adapt for your cluster, replace the `module load` commands and the
path to the LAMMPS executable with your local equivalents.
"""
