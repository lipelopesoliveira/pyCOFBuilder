# -*- coding: utf-8 -*-
"""
Created on Sat Oct 16 12:31:19 2021
@author: Felipe Lopes de Oliveira
"""

import os
import numpy as np
import tools as Tools
from tqdm import tqdm
from textwrap import dedent

def create_opt_input(database_pah, out_path, prefix, opt_n='opt_1', kpoints=False, smearing=False):
    
    try:
        os.mkdir(os.path.join(out_path, prefix))
    except: None

    try:
        os.mkdir(os.path.join(out_path, prefix, opt_n.upper()))
    except: None
    
    if opt_n.lower() == 'opt_1':
        max_opt = 20
    if opt_n.lower() == 'opt_2':
        max_opt = 500


    cof_json = Tools.read_json(database_pah, prefix)
    cof_json['system'][opt_n] = 'Runing'
    Tools.write_json(database_pah, prefix, cof_json)

    atom_labels = cof_json['geometry']['atom_labels']
    a, b, c, alpha, beta, gamma  = cof_json['geometry']['cell_parameters']
    cell_matrix = cof_json['geometry']['cell_matrix']
    atom_pos = cof_json['geometry']['atom_pos']


    kx, ky, kz = Tools.get_kgrid(cell_matrix, distance=0.3)

    opt_file = dedent(f"""&GLOBAL
  PROJECT {prefix}           !Project name. Output files will use this name
  RUN_TYPE CELL_OPT          !Calculation Type : MD (Mol. Dyn), GEO_OPT (Geom. Opt.), Energy (Ener. Calc.), CELL_OPT (Geom. and Cell Opt), VIBRATIONAL_ANALYSIS            
  IOLEVEL  MEDIUM            !Control the amount of information writed
  !PRINTLEVEL HIGH           !Control de amount of information printed   
&END GLOBAL

&FORCE_EVAL
  &DFT
    CHARGE 0
    MULTIPLICITY 1
    &QS
      METHOD xTB
      &xTB
         DO_EWALD T
         &PARAMETER
             DISPERSION_PARAMETER_FILE dftd3.dat
         &END PARAMETER
      &END XTB
    &END QS
    &POISSON
      &EWALD
       ALPHA 1.0
       EWALD_TYPE SPME
       GMAX 25
      &END EWALD
    &END POISSON
    """)
    if kpoints is True:
        opt_file += dedent(f"""
    &KPOINTS
       SCHEME  MONKHORST-PACK  {kx} {ky} {kz}
       SYMMETRY ON
       EPS_GEO 1.e-7
       FULL_GRID ON
       VERBOSE T
       PARALLEL_GROUP_SIZE   2
    &END KPOINTS""")
    opt_file += dedent(f"""
    &PRINT
         &HIRSHFELD OFF
         &END HIRSHFELD
         &LOWDIN OFF
         &END LOWDIN
         &MULLIKEN OFF
         &END MULLIKEN
      &END PRINT
    &SCF
      SCF_GUESS MOPAC
      MAX_SCF  150
      EPS_SCF 1.0E-7          !Accuracy of the SCF procedure typically 1.0E-6 - 1.0E-7
      &OUTER_SCF             !Repeat the inner SCF cycle 20 times
        MAX_SCF 20
        EPS_SCF 1.0E-7       !Musst match the above
      &END
      """)
    if kpoints is True:
        opt_file += dedent(f"""
      &MIXING
          METHOD KERKER_MIXING
          ALPHA   0.3
      &END""")
    if kpoints is False:
        opt_file += dedent(f"""
      &MIXING
          METHOD DIRECT_P_MIXING
          ALPHA   0.3
      &END""")
    if smearing is True:
        opt_file += dedent(f"""
      &SMEAR
         METHOD FERMI_DIRAC
         ELECTRONIC_TEMPERATURE  250
      &END
      ADDED_MOS 25""")
    opt_file += dedent(f"""
      
    &END SCF
  &END DFT
  STRESS_TENSOR ANALYTICAL 
  &SUBSYS    #System description
   &CELL
     ABC [angstrom] {a} {b} {c}   !Unit cells that are orthorhombic are more efficient with CP2K
     ALPHA_BETA_GAMMA [deg] {alpha} {beta} {gamma}        !Specify the angles between the vectors A, B and C when using the ABC keyword""")
    if opt_n == 'opt_1':
        opt_file += dedent(f"""
   SYMMETRY HEXAGONAL                    !Imposes an initial cell symmetry""")
    opt_file += dedent(f"""
   &END CELL
   &COORD
   """)

    for i in range(len(atom_labels)):
        opt_file += dedent(f"""     {atom_labels[i]}    {atom_pos[i][0]:>14.10f}  {atom_pos[i][1]:>14.10f}  {atom_pos[i][2]:>14.10f}
        """)
      
    opt_file += dedent(f"""&END
   &END SUBSYS
&END FORCE_EVAL
&MOTION
   &CELL_OPT
      &LBFGS
         TRUST_RADIUS [angstrom] 0.25
      &END LBFGS""")
    if opt_n == 'opt_1':
        opt_file += dedent(f"""
      KEEP_ANGLES .TRUE.""")
    opt_file += dedent(f"""
      MAX_DR [bohr] 0.030
      MAX_FORCE [bohr^-1*hartree] 0.0010
      MAX_ITER {max_opt}
      OPTIMIZER LBFGS
      RMS_DR [bohr] 0.015
      RMS_FORCE [bohr^-1*hartree] 0.0007
   &END CELL_OPT
   &PRINT
      &CELL
        &EACH
            CELL_OPT 1
            GEO_OPT 1
            MD 1
        &END EACH
      &END CELL
      &FORCES 
        &EACH
            CELL_OPT 1
            GEO_OPT 1
            MD 1
        &END EACH
      &END FORCES
      &RESTART
         BACKUP_COPIES 0
         &EACH
            CELL_OPT 1
            GEO_OPT 1
            MD 1
         &END EACH
      &END RESTART
      &RESTART_HISTORY OFF
      &END RESTART_HISTORY
      &STRESS
        &EACH
            CELL_OPT 1
            GEO_OPT 1
            MD 1
        &END EACH
      &END STRESS
      &TRAJECTORY
         &EACH
            CELL_OPT 1
            GEO_OPT 1
            MD 1
         &END EACH
         FORMAT XYZ
      &END TRAJECTORY
      &VELOCITIES
        &EACH
            CELL_OPT 1
            GEO_OPT 1
            MD 1
        &END EACH
      &END VELOCITIES
   &END PRINT
&END MOTION
!&EXT_RESTART
!  RESTART_FILE_NAME {prefix}-1.restart
!&END
""")

    file = open(os.path.join(out_path, prefix, opt_n.upper(), prefix + '.cell_opt.in'), 'w')
    
    for i in opt_file:
        file.write(i)
    file.close()

def create_bomd_input(database_pah, out_path, prefix, kpoints=False, smearing=False):

    try:
        os.mkdir(os.path.join(out_path, prefix))
    except: None

    try:
        os.mkdir(os.path.join(out_path, prefix, 'BOMD'))
    except: None

    cof_json = Tools.read_json(database_pah, prefix)
    cof_json['system']['bomd'] = 'Runing'
    Tools.write_json(database_pah, prefix, cof_json)

    atom_labels = cof_json['geometry']['atom_labels']
    a, b, c, alpha, beta, gamma  = cof_json['geometry']['cell_parameters']
    cell_matrix = cof_json['geometry']['cell_matrix']
    atom_pos = cof_json['geometry']['atom_pos']

    kx, ky, kz = Tools.get_kgrid(cell_matrix, distance=0.3)

    opt_file = dedent(f"""&GLOBAL
  PROJECT {prefix}           !Project name. Output files will use this name
  RUN_TYPE MD          !Calculation Type : MD (molecular dynamics), GEO_OPT (Geometry Optimization), Energy (Energy Calculation), CELL_OPT (Geometry and Cell Optimization), VIBRATIONAL_ANALYSIS            
  IOLEVEL  MEDIUM            !Control the amount of information writed
  !PRINTLEVEL HIGH           !Control de amount of information printed   
&END GLOBAL

&FORCE_EVAL
  &DFT
    CHARGE 0
    MULTIPLICITY 1
    &QS
      METHOD xTB
      &xTB
         DO_EWALD T
         &PARAMETER
             DISPERSION_PARAMETER_FILE dftd3.dat
         &END PARAMETER
      &END XTB
    &END QS
    &POISSON
      PERIODIC XYZ ! the default, gas phase systems should have 'NONE' and a wavelet solver
      &EWALD
       ALPHA 1.0
       EWALD_TYPE SPME
       GMAX 25
      &END EWALD
    &END POISSON
    """)
    if kpoints is True:
        opt_file += dedent(f"""
    &KPOINTS
       SCHEME  MONKHORST-PACK  {kx} {ky} {kz}
       SYMMETRY ON
       EPS_GEO 1.e-7
       FULL_GRID ON
       VERBOSE T
       PARALLEL_GROUP_SIZE   2
    &END KPOINTS""")
    opt_file += dedent(f"""
    &PRINT
         &HIRSHFELD OFF
         &END HIRSHFELD
         &LOWDIN OFF
         &END LOWDIN
         &MULLIKEN OFF
         &END MULLIKEN
      &END PRINT
    &SCF
      SCF_GUESS MOPAC
      MAX_SCF  150
      EPS_SCF 1.0E-7          !Accuracy of the SCF procedure typically 1.0E-6 - 1.0E-7
      &OUTER_SCF             !Repeat the inner SCF cycle 20 times
        MAX_SCF 20
        EPS_SCF 1.0E-7       !Musst match the above
      &END""")
    if kpoints is True:
        opt_file += dedent(f"""
      &MIXING
          METHOD KERKER_MIXING
          ALPHA   0.3
      &END""")
    if kpoints is False:
        opt_file += dedent(f"""
      &MIXING
          METHOD DIRECT_P_MIXING
          ALPHA   0.3
      &END""")
    if smearing is True:
        opt_file += dedent(f"""
      &SMEAR
         METHOD FERMI_DIRAC
         ELECTRONIC_TEMPERATURE  250
      &END
      ADDED_MOS 25""")
    opt_file += dedent(f"""
    &END SCF
  &END DFT
  STRESS_TENSOR ANALYTICAL 
  &SUBSYS    #System description
   &CELL
     ABC [angstrom] {a} {b} {c}   !Unit cells that are orthorhombic are more efficient with CP2K
     ALPHA_BETA_GAMMA [deg] {alpha} {beta} {gamma}        !Specify the angles between the vectors A, B and C when using the ABC keyword
   &END CELL
   &COORD
   """)

    for i in range(len(atom_labels)):
        opt_file += dedent(f"""     {atom_labels[i]}    {atom_pos[i][0]:>14.10f}  {atom_pos[i][1]:>14.10f}  {atom_pos[i][2]:>14.10f}
        """)
      
    opt_file += dedent(f"""&END
   &END SUBSYS
&END FORCE_EVAL
&MOTION
   &MD
      &BAROSTAT
        PRESSURE [bar] 1.0
        TIMECON [fs] 500
      &END BAROSTAT
      ENSEMBLE NPT_F
      STEPS 100
      TEMPERATURE 400
      &THERMOSTAT
        &CSVR
           TIMECON [fs] 0.1
        &END CSVR
        TYPE CSVR
      &END THERMOSTAT
      TIMESTEP [fs]  0.5
   &END MD
   &PRINT
      &CELL
        &EACH
            CELL_OPT 1
            GEO_OPT 1
            MD 1
        &END EACH
      &END CELL
      &FORCES 
        &EACH
            CELL_OPT 1
            GEO_OPT 1
            MD 1
        &END EACH
      &END FORCES
      &RESTART
         BACKUP_COPIES 0
         &EACH
            CELL_OPT 1
            GEO_OPT 1
            MD 1
         &END EACH
      &END RESTART
      &RESTART_HISTORY OFF
      &END RESTART_HISTORY
      &STRESS
        &EACH
            CELL_OPT 1
            GEO_OPT 1
            MD 1
        &END EACH
      &END STRESS
      &TRAJECTORY
         &EACH
            CELL_OPT 1
            GEO_OPT 1
            MD 1
         &END EACH
         FORMAT XYZ
      &END TRAJECTORY
      &VELOCITIES
        &EACH
            CELL_OPT 1
            GEO_OPT 1
            MD 1
        &END EACH
      &END VELOCITIES
   &END PRINT
&END MOTION
!&EXT_RESTART
!  RESTART_FILE_NAME {prefix}-1.restart
!&END
""")
    
    file = open(os.path.join(out_path, prefix, 'BOMD', prefix + '.bomd.in'), 'w')
    
    for i in opt_file:
        file.write(i)
    file.close()


def create_vib_input(database_path, out_path, prefix, kpoints=False, smearing=False):

    try:
        os.mkdir(os.path.join(out_path, prefix))
    except: None

    try:
        os.mkdir(os.path.join(out_path, prefix, 'PHONON'))
    except: None

    cof_json = Tools.read_json(database_path, prefix)
    cof_json['system']['phonon'] = 'Runing'
    Tools.write_json(database_path, prefix, cof_json)

    atom_labels = cof_json['geometry']['atom_labels']
    a, b, c, alpha, beta, gamma  = cof_json['geometry']['cell_parameters']
    cell_matrix = cof_json['geometry']['cell_matrix']
    atom_pos = cof_json['geometry']['atom_pos']

    kx, ky, kz = Tools.get_kgrid(cell_matrix, distance=0.3)

    opt_file = dedent(f"""&GLOBAL
  PROJECT {prefix}                !Project name. Output files will use this name
  RUN_TYPE VIBRATIONAL_ANALYSIS   !Calculation Type : MD (molecular dynamics), GEO_OPT (Geometry Optimization), Energy (Energy Calculation), CELL_OPT (Geometry and Cell Optimization), VIBRATIONAL_ANALYSIS            
  IOLEVEL  MEDIUM                 !Control the amount of information writed
  !PRINTLEVEL HIGH                !Control de amount of information printed   
&END GLOBAL

&VIBRATIONAL_ANALYSIS
    FULLY_PERIODIC
    INTENSITIES
    NPROC_REP 2
&END

&FORCE_EVAL
  &DFT
    CHARGE 0
    MULTIPLICITY 1
    &QS
      METHOD xTB
      &xTB
         DO_EWALD T
         &PARAMETER
             DISPERSION_PARAMETER_FILE dftd3.dat
         &END PARAMETER
      &END XTB
    &END QS
    &POISSON
      PERIODIC XYZ ! the default, gas phase systems should have 'NONE' and a wavelet solver
      &EWALD
       ALPHA 1.0
       EWALD_TYPE SPME
       GMAX 25
      &END EWALD
    &END POISSON
    """)
    if kpoints is True:
        opt_file += dedent(f"""
    &KPOINTS
       SCHEME  MONKHORST-PACK  {kx} {ky} {kz}
       SYMMETRY ON
       EPS_GEO 1.e-7
       FULL_GRID ON
       VERBOSE T
       PARALLEL_GROUP_SIZE   2
    &END KPOINTS""")
    opt_file += dedent(f"""
    &PRINT
         &MOMENTS
           PERIODIC TRUE
         &END
         &HIRSHFELD OFF
         &END HIRSHFELD
         &LOWDIN OFF
         &END LOWDIN
         &MULLIKEN OFF
         &END MULLIKEN
      &END PRINT
    &SCF
      SCF_GUESS MOPAC
      MAX_SCF  200
      EPS_SCF 1.0E-7          !Accuracy of the SCF procedure typically 1.0E-6 - 1.0E-7
      &OUTER_SCF             !Repeat the inner SCF cycle 20 times
        MAX_SCF 20
        EPS_SCF 1.0E-7       !Musst match the above
      &END""")
    if kpoints is True:
        opt_file += dedent(f"""
      &MIXING
          METHOD KERKER_MIXING
          ALPHA   0.3
      &END""")
    if kpoints is False:
        opt_file += dedent(f"""
      &MIXING
          METHOD DIRECT_P_MIXING
          ALPHA   0.3
      &END""")
    if smearing is True:
        opt_file += dedent(f"""
      &SMEAR
         METHOD FERMI_DIRAC
         ELECTRONIC_TEMPERATURE  250
      &END
      ADDED_MOS 25""")
    opt_file += dedent(f"""
    &END SCF
  &END DFT
  STRESS_TENSOR ANALYTICAL 
  &SUBSYS    #System description
   &CELL
     ABC [angstrom] {a} {b} {c}   !Unit cells that are orthorhombic are more efficient with CP2K
     ALPHA_BETA_GAMMA [deg] {alpha} {beta} {gamma}        !Specify the angles between the vectors A, B and C when using the ABC keyword
   &END CELL
   &COORD
   """)

    for i in range(len(atom_labels)):
        opt_file += dedent(f"""     {atom_labels[i]}    {atom_pos[i][0]:>14.10f}  {atom_pos[i][1]:>14.10f}  {atom_pos[i][2]:>14.10f}
        """)
      
    opt_file += dedent(f"""&END
   &END SUBSYS
&END FORCE_EVAL
!&EXT_RESTART
!  RESTART_FILE_NAME {prefix}-1.restart
!&END
""")
    
    file = open(os.path.join(out_path, prefix, 'PHONON', prefix + '.ph.in'), 'w')
    
    for i in opt_file:
        file.write(i)
    file.close()

