"""
build_Assembly_annotated.py
===========================
THIS SCRIPT IS NON-RUNNABLE. It is a heavily commented walkthrough of
build_Assembly.py, intended to explain what the LAMMPS input deck
(collagen.in) does and why. For the actual simulation script, see
build_Assembly.py.

The real script is a Python program that writes two files into a run
directory: `collagen.in` (the LAMMPS input deck) and `runscript.sh`
(the SLURM submission script). The content below mirrors the structure
of collagen.in section by section.

Overview
--------
Each collagen IV protomer is represented as a two-bead rod:
  - Type 1: NC1 domain bead (end-to-end binding, max 1 partner)
  - Type 2-8: 7S domain beads (lateral binding, max 3 partners)
  - Type 9: NC1 free monomer (unbound protomer)

Protomers self-assemble from a lattice initial configuration via
stochastic bond formation and breaking, implemented using LAMMPS'
`fix bond/react` package. The simulation runs in NVT (Langevin
thermostat + NVE integrator).
"""

# =============================================================================
# SECTION 1: COMMAND-LINE PARAMETERS
# =============================================================================
# build_Assembly.py reads all parameters from sys.argv so that Assemble.sh
# can sweep over them in a loop. Key parameters:
#
#   workingdir  -- output folder for this run
#   tEq         -- number of NVT timesteps for the assembly run
#   Frames      -- number of trajectory snapshots to save
#   Nevery      -- how often (in timesteps) bond reactions are attempted
#   tstep       -- integration timestep (LJ units)
#   MakeDist    -- NC1 bond formation distance cutoff (RmaxF)
#   BreakDist   -- NC1 bond breaking distance cutoff (RminFB)
#   MakeDist_7S -- 7S bond formation distance cutoff
#   BreakDist_7S-- 7S bond breaking distance cutoff
#   MakeProb    -- bond formation probability per attempt
#   FactorMult  -- BreakProb = FactorMult * MakeProb
#   KangNC1     -- angular stiffness of NC1 bonds (harmonic, 180° rest angle)
#   seed        -- random seed (varied across replicates)


# =============================================================================
# SECTION 2: LAMMPS INITIALISATION
# =============================================================================
"""
units           lj              # Lennard-Jones reduced units
boundary        p p p           # Periodic in all three dimensions
atom_style      molecular       # Atoms carry bond/angle topology

bond_style      harmonic        # U_bond = k(r - r0)^2
angle_style     harmonic        # U_angle = k(theta - theta0)^2

# Temporary LJ pair style just to start — replaced below after read_data
pair_style      lj/cut 1.0

# Read the initial configuration (no bonds formed yet).
# extra/bond/per/atom etc. pre-allocates memory for bonds that will form
# during the simulation via bond/react.
read_data  data  extra/bond/per/atom 10  extra/special/per/atom 100  extra/angle/per/atom 20
reset_timestep 0
"""


# =============================================================================
# SECTION 3: FORCE FIELD
# =============================================================================
"""
# Bond types:
#   1: stiff backbone bond holding the two beads of each protomer together
#      (k=200, r0=3.0 — very stiff, keeps bead separation fixed)
#   2: NC1-NC1 binding bond  (k=6, r0=0.5)
#   3: 7S-7S binding bond    (k=6, r0=0.5)
#   4: extra stiff bond for nucleation template (k=1000, r0=2.0) - these are dummy/ghost atoms
bond_coeff  1   200.0  3.0
bond_coeff  2   6.0    0.5
bond_coeff  3   6.0    0.5
bond_coeff  4   1000.0 2.0

# Angle types:
#   1: NC1 angle — stiffness KangNC1 (paper value: 4.0), rest angle 180°
#      This enforces approximate linearity of NC1-NC1 bonds
#   2: 7S angle — softer (k=1.0), rest angle 155°
#   3: 7S angle orthogonal — (k=1.0), rest angle 60° (for closed 7S trimers)
angle_coeff  1  KangNC1  180.0
angle_coeff  2  1.0      155.0
angle_coeff  3  1.0       60.0

# Pair style: WCA (purely repulsive) + zero potential.
# WCA = LJ truncated at r_min (cutoff = 2^(1/6) ≈ 1.122), giving only
# steric repulsion between all beads.
# The zero potential with a larger cutoff ensures all pairs within that
# range are included in the neighbour list, which bond/react needs to
# find candidate bonding partners.
pair_style   hybrid/overlay  lj/cut 1.0  zero 1.0
pair_coeff   * *  lj/cut  0.0  0.0  1.122462   # WCA: epsilon=0, sigma=0 except...
pair_coeff   * *  zero    6.0                   # ...zero potential out to 6.0
"""


# =============================================================================
# SECTION 4: MOLECULE TEMPLATES
# =============================================================================
"""
Each bond/react reaction needs three files per topology change:
  - XX_pre.txt  : the molecule graph BEFORE the reaction
  - XX_post.txt : the molecule graph AFTER the reaction
  - XX_map.txt  : atom mapping between pre and post states

The templates are organised by the bond topology they describe.
Multiple numbered variants (01, 02, ...) cover all distinct ways
that topology can appear given the current bond valencies of the
participating beads (i.e. how many other bonds each bead already has).

NC1 templates (end-to-end binding):
  NC1_2N/  -- NC1 bead with 2 existing NC1 neighbours gaining one more
  NC1_3N/  -- NC1 bead with 3 existing NC1 neighbours
  NC1_4N/ through NC1_8N/  -- higher valency states
  NC1_loop/ -- NC1 bead forming a loop back to itself

Example (abbreviated — the real script loads all variants):
  molecule  Nc1_2N_01_pre   ../../templates/NC1_2N/01_pre.txt
  molecule  Nc1_2N_01_post  ../../templates/NC1_2N/01_post.txt
  # ... (11 variants for NC1_2N alone)

7S templates (lateral binding, up to 3 bonds per bead):
  MonomerToDimer/          -- first 7S bond forms (monomer → dimer)
  DimerToLineTrimer/       -- second 7S bond forms (dimer → linear trimer)
  LineTrimerToTriangleTrimer/   -- linear trimer closes into triangle
  TriangleTrimerToOpenTetramer/ -- triangle gains a fourth protomer
  OpenTetramerToRhomboidTetramer/
  RhomboidTetramerToSquareTetramer/

Nucleation templates:
  Nucleation/pre_Nucleation.txt / post_Nucleation.txt
  Nucleation/pre_Death.txt      / post_Death.txt
  These handle creation/removal of type-9 free monomers.
"""


# =============================================================================
# SECTION 5: COMPUTES AND DIAGNOSTICS
# =============================================================================
"""
# Bond and angle property computes (for dump output):
compute  1  all  property/local  batom1 batom2 btype
compute  4  all  property/local  atype aatom1 aatom2 aatom3
compute  2  all  pe  bond
compute  3  all  pe  angle

# Per-atom stress tensor, summed to give total stress in each direction:
compute  perAtomStress  all  stress/atom  NULL
compute  totalStressX   all  reduce  sum  c_perAtomStress[1]
compute  totalStressY   all  reduce  sum  c_perAtomStress[2]
compute  totalStressZ   all  reduce  sum  c_perAtomStress[3]

# Average bond force magnitude:
compute  forceBond     all  bond/local  force
compute  forceBondAll  all  reduce  ave  c_forceBond

# Monomer count: type-9 atoms are free (unbound) protomers.
# This is tracked throughout the run to monitor equilibration.
group    type9      dynamic all  var vMaskType9  every 1
variable NumType9   equal   count(type9)
"""


# =============================================================================
# SECTION 6: NVT INTEGRATOR
# =============================================================================
"""
# Langevin thermostat + NVE gives an NVT ensemble.
# damp = 0.1 (friction coefficient in LJ units)
fix  fix_langevin  all  langevin  1.0  1.0  0.1  ${lseed}
fix  fNVE          all  nve

# Initialise velocities from Maxwell-Boltzmann at T=1.0, then rescale:
velocity  all  create  1.0  ${vseed}  dist gaussian  mom yes  rot yes
run 0
velocity  all  scale  1.0
"""


# =============================================================================
# SECTION 7: OUTPUT
# =============================================================================
"""
# Thermo file: written every `thermodump` steps.
# Columns: step, total energy, KE, bond PE, angle PE, temperature,
#          pressure, stress X/Y/Z, average bond force,
#          number of free monomers, box length X, box volume,
#          current monomer concentration, preferred concentration.
# (Note: preferred concentration is computed but GCE insertion/deletion
#  is not active in this stage — monomer count is monitored only.)
fix  fix_print  all  print  ${thermodump}  "..."  file  thermo.dat

# Trajectory dump: particle positions every `trajdump` steps.
dump  2  all  custom  ${trajdump}  dumplin/dump.${simname}.*.lammpstrj  type id xu yu zu

# Bond topology dump: bond index, atom IDs, bond type.
dump  3  all  local   ${trajdump}  dumplin_bonds/bonds.${simname}.*  index c_1[1] c_1[2] c_1[3]

# Restart files saved every trestart steps:
restart  ${trestart}  restart/${simname}.*.restart
"""


# =============================================================================
# SECTION 8: BOND/REACT FIX (ASSEMBLY)
# =============================================================================
"""
This is the core of the simulation. A single `fix bond/react` command
registers ALL allowed bond formation and breaking reactions simultaneously.

Structure of each reaction line:
  react  <name>  all  <Nevery>  <Rmin>  <Rmax>  <pre>  <post>  <map>  prob  <p>  <seed>

For bond FORMATION (forward reactions, prefixed 'R' or 'Mono', 'Dimer' etc.):
  - Rmin = 0 (no minimum distance)
  - Rmax = RmaxF (NC1) or RmaxF_7S (7S) — formation distance cutoff
  - prob = MakeProb

For bond BREAKING (reverse reactions, prefixed 'B'):
  - pre and post are SWAPPED relative to the forward reaction
  - Rmin = BreakDist (NC1) or BreakDist_7S (7S) — bond breaks if stretched beyond this
  - Rmax = 30 (large; effectively no upper limit)
  - prob = BreakProb = FactorMult * MakeProb

NC1 bond reactions (abbreviated — 11 variants for 2N, 10 for 3N, etc.):
  react Rc1_2N_01  all ${NeverySlow} 0 ${RmaxF}    Nc1_2N_01_pre  Nc1_2N_01_post  .../NC1_2N/01_map.txt   prob v_probF    ${lseed1}
  react BRc1_2N_01 all ${NeverySlow} ${RminFB} 30  Nc1_2N_01_post Nc1_2N_01_pre   .../NC1_2N/01_map.txt   prob v_probFB   ${lseed14}
  # ... (many more NC1 variants)

7S bond reactions (abbreviated):
  react MonoToDimer_01          all ${NeverySlow} 0 ${RmaxF_7S}  ...  prob v_probF_7S   ${lseed2}
  react DimerToLineTrimer_01    all ${NeverySlow} 0 ${RmaxF_7S}  ...  prob v_probF_7S   ${lseed3}
  react LineTrimerToTriangleTrimer_01  all ${NeveryFast} ...  # triangle closure is fast
  react TriangleTrimerToOpenTetramer_01 all ${NeverySlow} ...
  react OpenTetramerToRhomboidTetramer_01  all ${NeveryFast} ...
  react RhomboidTetramerToSquareTetramer_01 all ${NeveryFast} ...
  # ... and corresponding B (breaking) variants for each

Note on NeveryFast vs NeverySlow:
  - NeverySlow = Nevery (500) — used for bond formation/breaking events
    that require two beads to diffuse within range
  - NeveryFast = 1 — used for topological rearrangements (e.g. linear
    trimer → triangle trimer) that are purely geometric and should
    happen as soon as the geometry allows
"""


# =============================================================================
# SECTION 9: RUN AND OUTPUT
# =============================================================================
"""
run         ${tbonds}           # Run for tEq timesteps
write_data  data.postmakebonds  # Save final configuration with bond topology
write_restart restart.postmakebonds
"""


# =============================================================================
# SECTION 10: SLURM SUBMISSION SCRIPT (runscript.sh)
# =============================================================================
"""
The Python script also writes runscript.sh, which:
  1. Creates output subdirectories: dumplin/, dumplin_bonds/, restart/
  2. Generates random seeds for the Langevin thermostat and velocities
  3. Runs LAMMPS in serial (mpirun -np 1) on a single node

To adapt for your cluster, replace the `module load` commands and the
path to the LAMMPS executable with your local equivalents.
"""
