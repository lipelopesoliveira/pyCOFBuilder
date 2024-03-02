QuantumESPRESSO input (pwscf)
=========================================

QuantumESPRESSO is a suite of programs for electronic structure calculations based on density-functional theory (DFT).
The input file for QuantumESPRESSO, commonly named `input.in`, is crucial for specifying the details of the simulation,
including system information, calculation parameters, and method settings. Below is a simplified example of a
QuantumESPRESSO input file for a molecular structure:

.. code-block:: bash

   &CONTROL
   calculation = 'scf'
   restart_mode = 'from_scratch'
   pseudo_dir = '/path/to/pseudopotentials/'
   outdir = './tmp'
   prefix = 'molecule'
   disk_io = 'low'
   verbosity = 'high'
   /

   &SYSTEM
   ibrav = 1
   celldm(1) = 10.0
   nat = 3
   ntyp = 2
   ecutwfc = 30.0
   occupations = 'smearing'
   smearing = 'mv'
   degauss = 0.02
   /

   &ELECTRONS
   conv_thr = 1.0e-8
   mixing_beta = 0.7
   /

   ATOMIC_SPECIES
   O 16.00  O.pbe-mt_fhi.UPF
   H 1.008  H.pbe-mt_fhi.UPF

   ATOMIC_POSITIONS crystal
   O 0.0000 0.0000 0.0000
   H 0.7570 0.5860 0.0000
   H -0.7570 0.5860 0.0000

   K_POINTS automatic
   4 4 4 0 0 0


In this example:

- The `&CONTROL` block specifies global control parameters, such as the type of calculation (`'scf'` for self-consistent field), file I/O settings, and verbosity.
- The `&SYSTEM` block contains information about the crystal structure (`ibrav`, `celldm`), the number of atoms and types, energy cutoff for wavefunctions (`ecutwfc`), and parameters related to electronic occupations and smearing.
- The `&ELECTRONS` block sets convergence criteria for the electronic self-consistency loop.
- The `ATOMIC_SPECIES` block defines the atomic species, masses, and pseudopotential files.
- The `ATOMIC_POSITIONS` block specifies the atomic positions.
- The `K_POINTS` block determines the k-point sampling scheme.

This input file is for a self-consistent field (SCF) calculation of a simple molecular structure with oxygen and hydrogen atoms. 
Note that the specific details, such as pseudopotential files and k-point sampling, depend on the system being studied and
the desired level of accuracy. Always refer to the QuantumESPRESSO documentation for accurate and up-to-date information regarding input file parameters.
