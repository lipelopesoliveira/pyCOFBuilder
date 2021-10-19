import os
import numpy as np
import tools as Tools
from tqdm import tqdm


def get_available_runs():
    f = os.popen('qstat -f').readlines()
    user = os.popen('echo $USER').readlines()[0].rstrip('\n')
    jobid_pos = [i for i in range(len(f)) if 'Job Id:' in f[i]]

    ids = []
    job_names = []
    state = []

    for i in jobid_pos:
        if user in f[i + 2]:
            ids += [f[i].lstrip('Job Id:')[:-11]]
            job_names += [f[i+1].rstrip('\n').lstrip('    Job_Name = ')]
            if 'job_state' in f[i+3]:
                state += ['Q']
            if 'resources_used' in f[i+3]:
                state += ['R']
    rodando = 0
    fila = 0
    for i in state:
        if i == 'R':
            rodando +=1
        if i == 'Q':
            fila +=1

    n_run = 20 - fila

    return n_run, rodando, fila

def suc_finiched(file):
    temp = open(file).readlines()
    success = False
    for i in temp:
        if 'PROGRAM ENDED AT ' in i:
            success = True
    return success

def used_smearing(file):
    temp = open(file).readlines()
    smearing = True
    for i in temp:
        if '!&SMEAR' in i:
            smearing = False
    return smearing

def read_cp2k_result(path, prefix):
    pos_file = open(os.path.join(path, prefix + '-pos-1.xyz'), 'r').readlines()
    pos_file = [i.rstrip('\n') for i in pos_file]
    
    n_atoms = int(pos_file[0])
    
    energies = [float(i.split()[-1]) for i in pos_file if 'i' in i]
    
    splited_pos = []
    
    for i in range(len(energies)):
        splited_pos += [pos_file[i*(n_atoms + 2): (i+1)*(n_atoms + 2)]]
    
    for i in range(len(splited_pos)):
        splited_pos[i] = splited_pos[i][2:]
    
    for i in range(len(splited_pos)):
        for j in range(len(splited_pos[i])):
            splited_pos[i][j] = splited_pos[i][j].split()
            
    last_pos = np.transpose(splited_pos[-1])
    
    atom_list = last_pos[0]
    
    atom_pos = np.transpose(last_pos[1:])
    
    floated_pos = []
    for i in range(len(atom_pos)):
        floated_pos += [[float(j) for j in atom_pos[i]]]
    
    cell_file = open(os.path.join(path, prefix + '-1.cell'), 'r').readlines()[1:]
    cell_file = [i.rstrip('\n').split()[2:-1] for i in cell_file]
    
    for i in range(len(cell_file)):
        cell_file[i] = np.array([float(j) for j in cell_file[i]]).reshape(3,3)
        
    return energies, cell_file[-1], atom_list, floated_pos


#####################################################################


def cread_opt_input(database_pah, out_path, prefix, opt_n='opt_1', kpoints=False, smearing=False):
    
    try:
        os.mkdir(os.path.join(out_path, prefix))
    except: None

    try:
        os.mkdir(os.path.join(out_path, prefix, opt_n.upper()))
    except: None
    
    opt_2_label=''
    if opt_n.lower() == 'opt_1':
        max_opt = 20
    if opt_n.lower() == 'opt_2':
        opt_2_label='!'
        max_opt = 500


    cof_json = Tools.read_json(database_pah, prefix)
    cof_json['system'][opt_n] = 'Runing'
    Tools.write_json(database_pah, prefix, cof_json)

    atom_labels = cof_json['geometry']['atom_labels']
    a, b, c, alpha, beta, gamma  = cof_json['geometry']['cell_parameters']
    cell_matrix = cof_json['geometry']['cell_matrix']
    atom_pos = cof_json['geometry']['atom_pos']

    atom_types = list(set(atom_labels))

    kx, ky, kz = Tools.get_kgrid(cell_matrix, distance=0.3)

    k_string = '!'
    if kpoints is True:
        k_string = ''

    s_string = '!'
    OT_label = ''
    if smearing is True:
        s_string = ''
        OT_label = '!'
    
    opt_file = [ '&GLOBAL\n',
                 f'  PROJECT {prefix}   !Project name. Output files will use this name\n',
                 '  RUN_TYPE CELL_OPT !Calculation Type : MD (molecular dynamics), GEO_OPT (Geometry Optimization), Energy (Energy Calculation), CELL_OPT (Geometry and Cell Optimization)           \n',
                 '  !WALLTIME 1800    !Time limit in seconds for calculation. Finishes allowing restart   \n',
                 '  IOLEVEL  MEDIUM   !Control the amount of information writed\n',
                 '  !PRINTLEVEL HIGH  !Control de amount of information printed   \n',
                 '&END GLOBAL\n',
                 '\n',
                 '&FORCE_EVAL\n',
                 '  METHOD Quickstep   !Method to calculate force: Fist (Molecular Mechanics), QS or QUICKSTEP (Electronic structure methods, like DFT)\n',
                 '  &DFT\n',
                 '    ! basis sets and pseudopotential files can be found in cp2k/data\n',
                 '    BASIS_SET_FILE_NAME BASIS_MOLOPT\n',
                 '    BASIS_SET_FILE_NAME BASIS_MOLOPT_UCL\n',
                 '    POTENTIAL_FILE_NAME GTH_POTENTIALS            \n',
                 '\n',
                 '    ! Charge and multiplicity\n',
                 '    CHARGE 0\n',
                 '    MULTIPLICITY 1\n',
                 '\n',
                 '    &MGRID\n',
                 '       ! PW cutoff ... depends on the element (basis) too small cutoffs lead to the eggbox effect.\n',
                 '       ! Certain calculations (e.g. geometry optimization, vibrational frequencies, NPT and cell optimizations, need higher cutoffs) Squared value of QE\n',
                 '       CUTOFF [Ry] 400 \n',
                 '       NGRIDS 4\n',
                 '       REL_CUTOFF 50\n',
                 '    &END\n',
                 '\n',
                 '    &QS\n',
                 '       METHOD GPW            !Use the GPW method (i.e. pseudopotential based calculations with the Gaussian and Plane Waves scheme).\n',
                 '       EPS_DEFAULT 1.0E-10   !Threshold for numerics ~ roughly numerical accuracy of the total energy per electron sets reasonable values for all other thresholds.\n',
                 '       EXTRAPOLATION ASPC    !Used for MD, the method used to generate the initial guess.\n',
                 '    &END\n',
                f'    {k_string}&KPOINTS\n'
                f'    {k_string}    SCHEME  MONKHORST-PACK  {kx} {ky} {kz}\n'
                f'    {k_string}&END KPOINTS\n',
                 '    &POISSON\n',
                 "       PERIODIC XYZ ! the default, gas phase systems should have 'NONE' and a wavelet solver\n",
                 '    &END\n',
                 '\n',
                 '    &PRINT\n',
                 '         &HIRSHFELD OFF\n',
                 '         &END HIRSHFELD\n',
                 '         &LOWDIN OFF\n',
                 '         &END LOWDIN\n',
                 '         &MO_CUBES \n',
                 '            ADD_LAST SYMBOLIC\n',
                 '            &EACH\n',
                 '               CELL_OPT 0\n',
                 '               GEO_OPT 0\n',
                 '               MD 0\n',
                 '            &END EACH\n',
                 '            NHOMO 1\n',
                 '            NLUMO 1\n',
                 '            WRITE_CUBE F\n',
                 '         &END MO_CUBES\n',
                 '         &MULLIKEN OFF\n',
                 '         &END MULLIKEN\n',
                 '      &END PRINT\n',
                 '\n',
                 '    ! use the OT METHOD for robust and efficient SCF, suitable for all non-metallic systems.\n',
                 '    &SCF                      !Use the OT METHOD for robust and efficient SCF, suitable for all non-metallic systems.                              \n',
                 '      SCF_GUESS ATOMIC        !#Can be used to RESTART an interrupted calculation\n',
                 '      MAX_SCF 50\n',
                f'      {s_string}ADDED_MOS 25\n',
                f'      {s_string}&SMEAR\n',
                f'      {s_string}  METHOD FERMI_DIRAC\n',
                f'      {s_string}  ELECTRONIC_TEMPERATURE 300\n',
                f'      {s_string}&END SMEAR\n',
                 '      EPS_SCF 1.0E-7          !Accuracy of the SCF procedure typically 1.0E-6 - 1.0E-7\n',
                 '      MAX_ITER_LUMO 10000\n',
                 '      &MIXING\n',
                 '       ALPHA 0.4\n',
                 '       BETA 0.5\n',
                 '       METHOD BROYDEN_MIXING\n',
                 '      &END MIXING\n',
                f'     {OT_label}&OT\n',
                f'       {OT_label}PRECONDITIONER FULL_SINGLE_INVERSE   ! OR (FULL_ALL) An accurate preconditioner suitable also for larger systems\n',
                f'       {OT_label}MINIMIZER CG                         !The most robust choice (DIIS might sometimes be faster, but not as stable)\n',
                f'     {OT_label}&END OT\n',
                 '      &OUTER_SCF             !Repeat the inner SCF cycle 20 times\n',
                 '        MAX_SCF 50\n',
                 '        EPS_SCF 1.0E-7       !Musst match the above\n',
                 '      &END\n',
                 '      ! do not store the wfn during MD\n',
                 '      &PRINT\n',
                 '        &RESTART ON          !Set on if you want to generate the file NAME.restar\n',
                 '        &END\n',
                 '      &END\n',
                 '    &END SCF\n',
                 '\n',
                 '    &XC                     !Specify the exchange and correlation treatment\n',
                 '      &XC_FUNCTIONAL        !Use PBE functional for XC\n',
                 '         &PBE\n',
                 '            PARAMETRIZATION ORIG\n',
                 '         &END PBE\n',
                 '      &END XC_FUNCTIONAL\n',
                 '      &VDW_POTENTIAL        !Add a dispersion correction. \n',
                 '         POTENTIAL_TYPE PAIR_POTENTIAL \n',
                 "         &PAIR_POTENTIAL                      !Adding Grimme's D3 correction (by default without C9 terms)\n",
                 '            PARAMETER_FILE_NAME dftd3.dat   \n',
                 '            TYPE DFTD3(BJ)\n',
                 '            REFERENCE_FUNCTIONAL PBE\n',
                 '            R_CUTOFF [angstrom] 16\n',
                 '         &END\n',
                 '      &END VDW_POTENTIAL\n',
                 '    &END XC\n',
                 '  &END DFT\n',
                 '  STRESS_TENSOR ANALYTICAL \n',
                 '  &SUBSYS    #System description\n',
                 '    &CELL \n',
                 f'      ABC [angstrom] {a} {b} {c}   !Unit cells that are orthorhombic are more efficient with CP2K\n',
                 f'      ALPHA_BETA_GAMMA [deg] {alpha} {beta} {gamma}        !Specify the angles between the vectors A, B and C when using the ABC keyword\n',
                 f'      {opt_2_label}SYMMETRY HEXAGONAL                    !Imposes an initial cell symmetry\n',
                 '    &END CELL\n',
                 '\n',
                 '    #Atom coordinates can be in the &COORD section or provided as an external file on TOPOLOGY.\n',
                 '    &COORD\n']

    for i in range(len(atom_labels)):
        opt_file += [f'      {atom_labels[i]}    {atom_pos[i][0]:>14.10f}  {atom_pos[i][1]:>14.10f}  {atom_pos[i][2]:>14.10f}\n']
        
    opt_file += ['    &END COORD\n',
                 '\n',
                 "    ! MOLOPT basis sets are fairly costly, but in the 'DZVP-MOLOPT-SR-GTH' available for all elements\n",
                 '    ! their contracted nature makes them suitable for condensed and gas phase systems alike.\n']
    if 'H' in atom_types:
        opt_file += ['    &KIND H                              \n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q1\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q1             \n',
                 '    &END KIND\n']
    if 'B' in atom_types:
        opt_file += ['    &KIND B                              \n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q3\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q3             \n',
                 '    &END KIND\n']

    if 'Be' in atom_types:
        opt_file += ['    &KIND Be                              \n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q4\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q4             \n',
                 '    &END KIND\n']
    if 'C' in atom_types:
        opt_file += ['    &KIND C\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q4\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q4\n',
                 '    &END KIND\n']
    if 'N' in atom_types:   
        opt_file += ['    &KIND N\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q5\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q5\n',
                 '    &END KIND\n']
    if 'O' in atom_types:   
        opt_file += ['    &KIND O\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q6\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q6\n',
                 '    &END KIND  \n']
    if 'F' in atom_types:   
        opt_file += ['    &KIND F\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q7\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q7\n',
                 '    &END KIND  \n']
    if 'S' in atom_types:   
        opt_file += ['    &KIND S\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q6\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q6\n',
                 '    &END KIND  \n']
                
    opt_file += ['  &END SUBSYS\n',
                 '&END FORCE_EVAL\n',
                 '\n',
                 '&MOTION\n',
                 '   &CELL_OPT\n',
                 '      &LBFGS\n',
                 '         TRUST_RADIUS [angstrom] 0.25\n',
                 '      &END LBFGS\n',
                 '      KEEP_ANGLES .TRUE.\n',
                 '      MAX_DR [bohr] 0.030\n',
                 '      MAX_FORCE [bohr^-1*hartree] 0.0010\n',
                f'      MAX_ITER {max_opt}\n',
                 '      OPTIMIZER LBFGS\n',
                 '      RMS_DR [bohr] 0.015\n',
                 '      RMS_FORCE [bohr^-1*hartree] 0.0007\n',
                 '   &END CELL_OPT\n',
                 '   &PRINT\n',
                 '      &CELL\n',
                 '         &EACH\n',
                 '            CELL_OPT 1\n',
                 '            GEO_OPT 1\n',
                 '            MD 1\n',
                 '         &END EACH\n',
                 '      &END CELL\n',
                 '      &FORCES OFF\n',
                 '      &END FORCES\n',
                 '      &RESTART\n',
                 '         BACKUP_COPIES 0\n',
                 '         &EACH\n',
                 '            CELL_OPT 1\n',
                 '            GEO_OPT 1\n',
                 '            MD 1\n',
                 '         &END EACH\n',
                 '      &END RESTART\n',
                 '      &RESTART_HISTORY OFF\n',
                 '      &END RESTART_HISTORY\n',
                 '      &STRESS OFF\n',
                 '      &END STRESS\n',
                 '      &TRAJECTORY\n',
                 '         &EACH\n',
                 '            CELL_OPT 1\n',
                 '            GEO_OPT 1\n',
                 '            MD 1\n',
                 '         &END EACH\n',
                 '         FORMAT XYZ\n',
                 '      &END TRAJECTORY\n',
                 '      &VELOCITIES OFF\n',
                 '      &END VELOCITIES\n',
                 '   &END PRINT\n',
                 '&END MOTION\n',
                 '!&EXT_RESTART\n',
                 f'!  RESTART_FILE_NAME {prefix}-1.restart\n',
                 '!&END\n']

    file = open(os.path.join(out_path, prefix, opt_n.upper(), prefix + '.cell_opt.in'), 'w')
    
    for i in opt_file:
        file.write(i)
    file.close()

def cread_bomd_input(database_pah, out_path, prefix, kpoints=False, smearing=False):

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

    atom_types = list(set(atom_labels))

    kx, ky, kz = Tools.get_kgrid(cell_matrix, distance=0.3)
    k_string = '!'

    if kpoints is True:
        k_string = ''
    s_string = '!'
    OT_label = ''
    if smearing is True:
        s_string = ''
        OT_label = '!'
    
    opt_file = [ '&GLOBAL\n',
                 f'  PROJECT {prefix}   !Project name. Output files will use this name\n',
                 '  RUN_TYPE MD !Calculation Type : MD (molecular dynamics), GEO_OPT (Geometry Optimization), Energy (Energy Calculation), CELL_OPT (Geometry and Cell Optimization)           \n',
                 '  !WALLTIME 1800    !Time limit in seconds for calculation. Finishes allowing restart   \n',
                 '  IOLEVEL  MEDIUM   !Control the amount of information writed\n',
                 '  !PRINTLEVEL HIGH  !Control de amount of information printed   \n',
                 '&END GLOBAL\n',
                 '\n',
                 '&FORCE_EVAL\n',
                 '  METHOD Quickstep   !Method to calculate force: Fist (Molecular Mechanics), QS or QUICKSTEP (Electronic structure methods, like DFT)\n',
                 '  &DFT\n',
                 '    ! basis sets and pseudopotential files can be found in cp2k/data\n',
                 '    BASIS_SET_FILE_NAME BASIS_MOLOPT\n',
                 '    BASIS_SET_FILE_NAME BASIS_MOLOPT_UCL\n',
                 '    POTENTIAL_FILE_NAME GTH_POTENTIALS            \n',
                 '\n',
                 '    ! Charge and multiplicity\n',
                 '    CHARGE 0\n',
                 '    MULTIPLICITY 1\n',
                 '\n',
                 '    &MGRID\n',
                 '       ! PW cutoff ... depends on the element (basis) too small cutoffs lead to the eggbox effect.\n',
                 '       ! Certain calculations (e.g. geometry optimization, vibrational frequencies, NPT and cell optimizations, need higher cutoffs) Squared value of QE\n',
                 '       CUTOFF [Ry] 400 \n',
                 '       NGRIDS 4\n',
                 '       REL_CUTOFF 50\n',
                 '    &END\n',
                 '\n',
                 '    &QS\n',
                 '       METHOD GPW            !Use the GPW method (i.e. pseudopotential based calculations with the Gaussian and Plane Waves scheme).\n',
                 '       EPS_DEFAULT 1.0E-10   !Threshold for numerics ~ roughly numerical accuracy of the total energy per electron sets reasonable values for all other thresholds.\n',
                 '       EXTRAPOLATION ASPC    !Used for MD, the method used to generate the initial guess.\n',
                 '    &END\n',
                f'    {k_string}&KPOINTS\n'
                f'    {k_string}    SCHEME  MONKHORST-PACK  {kx} {ky} {kz}\n'
                f'    {k_string}&END KPOINTS\n',
                 '    &POISSON\n',
                 "       PERIODIC XYZ ! the default, gas phase systems should have 'NONE' and a wavelet solver\n",
                 '    &END\n',
                 '\n',
                 '    &PRINT\n',
                 '         &HIRSHFELD OFF\n',
                 '         &END HIRSHFELD\n',
                 '         &LOWDIN OFF\n',
                 '         &END LOWDIN\n',
                 '         &MO_CUBES\n',
                 '            ADD_LAST SYMBOLIC\n',
                 '            &EACH\n',
                 '               CELL_OPT 0\n',
                 '               GEO_OPT 0\n',
                 '               MD 0\n',
                 '            &END EACH\n',
                 '            NHOMO 1\n',
                 '            NLUMO 1\n',
                 '            WRITE_CUBE F\n',
                 '         &END MO_CUBES\n',
                 '         &MULLIKEN ON\n',
                 '            ADD_LAST SYMBOLIC\n',
                 '            &EACH\n',
                 '               CELL_OPT 0\n',
                 '               GEO_OPT 0\n',
                 '               MD 0\n',
                 '            &END EACH\n',
                 '         &END MULLIKEN\n',
                 '      &END PRINT\n',
                 '\n',
                 '    ! use the OT METHOD for robust and efficient SCF, suitable for all non-metallic systems.\n',
                 '    &SCF                      !Use the OT METHOD for robust and efficient SCF, suitable for all non-metallic systems.                              \n',
                 '      SCF_GUESS ATOMIC        !#Can be used to RESTART an interrupted calculation\n',
                 '      MAX_SCF 50\n',
                f'      {s_string}ADDED_MOS 25\n',
                f'      {s_string}&SMEAR\n',
                f'      {s_string}  METHOD FERMI_DIRAC\n',
                f'      {s_string}  ELECTRONIC_TEMPERATURE 300\n',
                f'      {s_string}&END SMEAR\n',
                 '      EPS_SCF 1.0E-7          !Accuracy of the SCF procedure typically 1.0E-6 - 1.0E-7\n',
                 '      MAX_ITER_LUMO 10000\n',
                 '      &MIXING\n',
                 '       ALPHA 0.4\n',
                 '       BETA 0.5\n',
                 '       METHOD BROYDEN_MIXING\n',
                 '      &END MIXING\n',
                f'     {OT_label}&OT\n',
                f'       {OT_label}PRECONDITIONER FULL_SINGLE_INVERSE   ! OR (FULL_ALL) An accurate preconditioner suitable also for larger systems\n',
                f'       {OT_label}MINIMIZER CG                         !The most robust choice (DIIS might sometimes be faster, but not as stable)\n',
                f'     {OT_label}&END OT\n',
                 '      &OUTER_SCF             !Repeat the inner SCF cycle 20 times\n',
                 '        MAX_SCF 50\n',
                 '        EPS_SCF 1.0E-7       !Musst match the above\n',
                 '      &END\n',
                 '      ! do not store the wfn during MD\n',
                 '      &PRINT\n',
                 '        &RESTART ON          !Set on if you want to generate the file NAME.restar\n',
                 '        &END\n',
                 '      &END\n',
                 '    &END SCF\n',
                 '\n',
                 '    &XC                     !Specify the exchange and correlation treatment\n',
                 '      &XC_FUNCTIONAL        !Use PBE functional for XC\n',
                 '         &PBE\n',
                 '            PARAMETRIZATION ORIG\n',
                 '         &END PBE\n',
                 '      &END XC_FUNCTIONAL\n',
                 '      &VDW_POTENTIAL        !Add a dispersion correction. \n',
                 '         POTENTIAL_TYPE PAIR_POTENTIAL \n',
                 "         &PAIR_POTENTIAL                      !Adding Grimme's D3 correction (by default without C9 terms)\n",
                 '            PARAMETER_FILE_NAME dftd3.dat   \n',
                 '            TYPE DFTD3(BJ)\n',
                 '            REFERENCE_FUNCTIONAL PBE\n',
                 '            R_CUTOFF [angstrom] 16\n',
                 '         &END\n',
                 '      &END VDW_POTENTIAL\n',
                 '    &END XC\n',
                 '  &END DFT\n',
                 '  STRESS_TENSOR ANALYTICAL \n',
                 '  &SUBSYS    #System description\n',
                 '    &CELL \n',
                 f'      ABC [angstrom] {a} {b} {c}   !Unit cells that are orthorhombic are more efficient with CP2K\n',
                 f'      ALPHA_BETA_GAMMA [deg] {alpha} {beta} {gamma}        !Specify the angles between the vectors A, B and C when using the ABC keyword\n',
                 '      SYMMETRY HEXAGONAL                    !Imposes an initial cell symmetry\n',
                 '    &END CELL\n',
                 '\n',
                 '    #Atom coordinates can be in the &COORD section or provided as an external file on TOPOLOGY.\n',
                 '    &COORD\n']
    for i in range(len(atom_labels)):
        opt_file += [f'      {atom_labels[i]}    {atom_pos[i][0]:>14.10f}  {atom_pos[i][1]:>14.10f}  {atom_pos[i][2]:>14.10f}\n']
        
    opt_file += ['    &END COORD\n',
                 '\n',
                 "    ! MOLOPT basis sets are fairly costly, but in the 'DZVP-MOLOPT-SR-GTH' available for all elements\n",
                 '    ! their contracted nature makes them suitable for condensed and gas phase systems alike.\n']
    if 'H' in atom_types:
        opt_file += ['    &KIND H                              \n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q1\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q1             \n',
                 '    &END KIND\n']
    if 'B' in atom_types:
        opt_file += ['    &KIND B                              \n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q3\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q3             \n',
                 '    &END KIND\n']

    if 'Be' in atom_types:
        opt_file += ['    &KIND Be                              \n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q4\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q4             \n',
                 '    &END KIND\n']
    if 'C' in atom_types:
        opt_file += ['    &KIND C\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q4\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q4\n',
                 '    &END KIND\n']
    if 'N' in atom_types:   
        opt_file += ['    &KIND N\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q5\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q5\n',
                 '    &END KIND\n']
    if 'O' in atom_types:   
        opt_file += ['    &KIND O\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q6\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q6\n',
                 '    &END KIND  \n']
    if 'F' in atom_types:   
        opt_file += ['    &KIND F\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q7\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q7\n',
                 '    &END KIND  \n']
    if 'S' in atom_types:   
        opt_file += ['    &KIND S\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q6\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q6\n',
                 '    &END KIND  \n']
                
    opt_file += ['  &END SUBSYS\n',
                     '&END FORCE_EVAL\n',
                     '\n',
                     '&MOTION\n',
                     '   &MD\n',
                     '      &BAROSTAT\n',
                     '         PRESSURE [bar] 1.0\n',
                     '         TIMECON [fs] 500\n',
                     '      &END BAROSTAT\n',
                     '      ENSEMBLE NPT_F\n',
                     '      STEPS 100\n',
                     '      TEMPERATURE 400\n',
                     '      &THERMOSTAT\n',
                     '         &CSVR\n',
                     '            TIMECON [fs] 0.1\n',
                     '         &END CSVR\n',
                     '         TYPE CSVR\n',
                     '      &END THERMOSTAT\n',
                     '      TIMESTEP [fs]  0.5\n',
                     '   &END MD\n',
                     '   &PRINT\n',
                     '      &CELL\n',
                     '         &EACH\n',
                     '            CELL_OPT 1\n',
                     '            GEO_OPT 1\n',
                     '            MD 1\n',
                     '         &END EACH\n',
                     '      &END CELL\n',
                     '      &FORCES OFF\n',
                     '      &END FORCES\n',
                     '      &RESTART\n',
                     '         BACKUP_COPIES 0\n',
                     '         &EACH\n',
                     '            CELL_OPT 1\n',
                     '            GEO_OPT 1\n',
                     '            MD 1\n',
                     '         &END EACH\n',
                     '      &END RESTART\n',
                     '      &RESTART_HISTORY OFF\n',
                     '      &END RESTART_HISTORY\n',
                     '      &STRESS OFF\n',
                     '      &END STRESS\n',
                     '      &TRAJECTORY\n',
                     '         &EACH\n',
                     '            CELL_OPT 1\n',
                     '            GEO_OPT 1\n',
                     '            MD 1\n',
                     '         &END EACH\n',
                     '         FORMAT XYZ\n',
                     '      &END TRAJECTORY\n',
                     '      &VELOCITIES OFF\n',
                     '      &END VELOCITIES\n',
                     '   &END PRINT\n',
                     '&END MOTION\n',
                     '!&EXT_RESTART\n',
                     f'!  RESTART_FILE_NAME {prefix}-1.restart\n',
                     '!&END\n']
    
    file = open(os.path.join(out_path, prefix, 'BOMD', prefix + '.bomd.in'), 'w')
    
    for i in opt_file:
        file.write(i)
    file.close()
    

def cread_charges_input(database_pah, out_path, prefix, kpoints=False, smearing=False, pbe0=True):
    
    try:
        os.mkdir(os.path.join(out_path, prefix))
    except: None

    try:
        os.mkdir(os.path.join(out_path, prefix, 'CHARGES'))
    except: None

    cof_json = Tools.read_json(database_pah, prefix)
    cof_json['system']['DDEC_charges'] = 'Runing'
    Tools.write_json(database_pah, prefix, cof_json)

    atom_labels = cof_json['geometry']['atom_labels']
    a, b, c, alpha, beta, gamma  = cof_json['geometry']['cell_parameters']
    cell_matrix = cof_json['geometry']['cell_matrix']
    atom_pos = cof_json['geometry']['atom_pos']

    atom_types = list(set(atom_labels))

    kx, ky, kz = Tools.get_kgrid(cell_matrix, distance=0.3)

    k_string = '!'
    if kpoints is True:
        k_string = ''
    
    pbe0_string = '!'
    if pbe0 is True:
        pbe0_string = ''

    s_string = '!'
    OT_label = ''
    if smearing is True:
        s_string = ''
        OT_label = '!'
    
    opt_file = [ '&GLOBAL\n',
                 f'  PROJECT {prefix}   !Project name. Output files will use this name\n',
                 '  RUN_TYPE Energy !Calculation Type : MD (molecular dynamics), GEO_OPT (Geometry Optimization), Energy (Energy Calculation), CELL_OPT (Geometry and Cell Optimization)           \n',
                 '  !WALLTIME 1800    !Time limit in seconds for calculation. Finishes allowing restart   \n',
                 '  IOLEVEL  MEDIUM   !Control the amount of information writed\n',
                 '  !PRINTLEVEL HIGH  !Control de amount of information printed   \n',
                 '&END GLOBAL\n',
                 '\n',
                 '&FORCE_EVAL\n',
                 '  METHOD Quickstep   !Method to calculate force: Fist (Molecular Mechanics), QS or QUICKSTEP (Electronic structure methods, like DFT)\n',
                 '  &DFT\n',
                 '    ! basis sets and pseudopotential files can be found in cp2k/data\n',
                 '    BASIS_SET_FILE_NAME BASIS_MOLOPT\n',
                 '    BASIS_SET_FILE_NAME BASIS_MOLOPT_UCL\n',
                f'   {pbe0_string} BASIS_SET_FILE_NAME BASIS_ADMM\n',
                 '    POTENTIAL_FILE_NAME GTH_POTENTIALS            \n',
                 '\n',
                 '    ! Charge and multiplicity\n',
                 '    CHARGE 0\n',
                 '    MULTIPLICITY 1\n',
                f'   {pbe0_string} &AUXILIARY_DENSITY_MATRIX_METHOD\n',
                f'   {pbe0_string}    METHOD BASIS_PROJECTION\n',
                f'   {pbe0_string}    ADMM_PURIFICATION_METHOD MO_DIAG\n',
                f'   {pbe0_string} &END AUXILIARY_DENSITY_MATRIX_METHOD\n',
                 '\n',
                 '    &MGRID\n',
                 '       ! PW cutoff ... depends on the element (basis) too small cutoffs lead to the eggbox effect.\n',
                 '       ! Certain calculations (e.g. geometry optimization, vibrational frequencies, NPT and cell optimizations, need higher cutoffs) Squared value of QE\n',
                 '       CUTOFF [Ry] 600 \n',
                 '       NGRIDS 4\n',
                 '       REL_CUTOFF 50\n',
                 '    &END\n',
                 '\n',
                 '    &QS\n',
                 '       METHOD GPW            !Use the GPW method (i.e. pseudopotential based calculations with the Gaussian and Plane Waves scheme).\n',
                 '       EPS_DEFAULT 1.0E-10   !Threshold for numerics ~ roughly numerical accuracy of the total energy per electron sets reasonable values for all other thresholds.\n',
                 '       EXTRAPOLATION ASPC    !Used for MD, the method used to generate the initial guess.\n',
                f'      {pbe0_string} EPS_PGF_ORB 1.0E-9\n',
                 '    &END\n',
                f'    {k_string}&KPOINTS\n'
                f'    {k_string}    SCHEME  MONKHORST-PACK  {kx} {ky} {kz}\n'
                f'    {k_string}&END KPOINTS\n',
                 '    &POISSON\n',
                 "       PERIODIC XYZ ! the default, gas phase systems should have 'NONE' and a wavelet solver\n",
                 '    &END\n',
                 '\n',
                 '    &PRINT\n',
                 '         &HIRSHFELD ON\n',
                 '            ADD_LAST SYMBOLIC\n',
                 '            &EACH\n',
                 '               CELL_OPT 0\n',
                 '               GEO_OPT 0\n',
                 '               MD 0\n',
                 '            &END EACH\n',
                 '         &END HIRSHFELD\n',
                 '         &LOWDIN ON\n',
                 '            ADD_LAST SYMBOLIC\n',
                 '            &EACH\n',
                 '               CELL_OPT 0\n',
                 '               GEO_OPT 0\n',
                 '               MD 0\n',
                 '            &END EACH\n',
                 '         &END LOWDIN\n',
                 '         &MO_CUBES\n',
                 '            ADD_LAST SYMBOLIC\n',
                 '            &EACH\n',
                 '               CELL_OPT 0\n',
                 '               GEO_OPT 0\n',
                 '               MD 0\n',
                 '            &END EACH\n',
                 '            NHOMO 1\n',
                 '            NLUMO 1\n',
                 '            WRITE_CUBE F\n',
                 '         &END MO_CUBES\n',
                 '         &MULLIKEN ON\n',
                 '            ADD_LAST SYMBOLIC\n',
                 '            &EACH\n',
                 '               CELL_OPT 0\n',
                 '               GEO_OPT 0\n',
                 '               MD 0\n',
                 '            &END EACH\n',
                 '         &END MULLIKEN\n',
                 '         &E_DENSITY_CUBE ON\n',
                 '            ADD_LAST SYMBOLIC\n',
                 '            FILENAME valence_density\n',
                 '           STRIDE 1\n',
                 '            &EACH\n',
                 '               CELL_OPT 0\n',
                 '               GEO_OPT 0\n',
                 '               MD 0\n',
                 '            &END EACH\n',
                 '         &END E_DENSITY_CUBE\n',
                 '      &END PRINT\n',
                 '\n',
                 '    ! use the OT METHOD for robust and efficient SCF, suitable for all non-metallic systems.\n',
                 '    &SCF                      !Use the OT METHOD for robust and efficient SCF, suitable for all non-metallic systems.                              \n',
                 '      SCF_GUESS ATOMIC        !#Can be used to RESTART an interrupted calculation\n',
                 '      MAX_SCF 50\n',
                f'      {s_string}ADDED_MOS 25\n',
                f'      {s_string}&SMEAR\n',
                f'      {s_string}  METHOD FERMI_DIRAC\n',
                f'      {s_string}  ELECTRONIC_TEMPERATURE 300\n',
                f'      {s_string}&END SMEAR\n',
                 '      EPS_SCF 1.0E-7          !Accuracy of the SCF procedure typically 1.0E-6 - 1.0E-7\n',
                 '      MAX_ITER_LUMO 10000\n',
                 '      &MIXING\n',
                 '       ALPHA 0.4\n',
                 '       BETA 0.5\n',
                 '       METHOD BROYDEN_MIXING\n',
                 '      &END MIXING\n',
                f'     {OT_label}&OT\n',
                f'       {OT_label}PRECONDITIONER FULL_SINGLE_INVERSE   ! OR (FULL_ALL) An accurate preconditioner suitable also for larger systems\n',
                f'       {OT_label}MINIMIZER CG                         !The most robust choice (DIIS might sometimes be faster, but not as stable)\n',
                f'     {OT_label}&END OT\n',
                 '      &OUTER_SCF             !Repeat the inner SCF cycle 20 times\n',
                 '        MAX_SCF 50\n',
                 '        EPS_SCF 1.0E-7       !Musst match the above\n',
                 '      &END\n',
                 '      ! do not store the wfn during MD\n',
                 '      &PRINT\n',
                 '        &RESTART ON          !Set on if you want to generate the file NAME.restar\n',
                 '        &END\n',
                 '      &END\n',
                 '    &END SCF\n',
                 '\n',
                 '    &XC                     !Specify the exchange and correlation treatment\n',
                 '      &XC_FUNCTIONAL        !Use PBE functional for XC\n',
                 '         &PBE\n',
                 '            PARAMETRIZATION ORIG\n',
                f'          {pbe0_string} SCALE_X 0.75\n',
                f'          {pbe0_string} SCALE_C 1.0\n',
                 '         &END PBE\n',
                f'       {pbe0_string} &PBE_HOLE_T_C_LR\n',
                f'       {pbe0_string}    CUTOFF_RADIUS 2.0\n',
                f'       {pbe0_string}    SCALE_X 0.25\n',
                f'       {pbe0_string} &END PBE_HOLE_T_C_LR\n',
                 '      &END XC_FUNCTIONAL\n',
                 '      &VDW_POTENTIAL        !Add a dispersion correction. \n',
                 '         POTENTIAL_TYPE PAIR_POTENTIAL \n',
                 "         &PAIR_POTENTIAL                      !Adding Grimme's D3 correction (by default without C9 terms)\n",
                 '            PARAMETER_FILE_NAME dftd3.dat   \n',
                 '            TYPE DFTD3(BJ)\n',
                 '            REFERENCE_FUNCTIONAL PBE\n',
                 '            R_CUTOFF [angstrom] 16\n',
                 '         &END\n',
                 '      &END VDW_POTENTIAL\n',
                f'   {pbe0_string} &HF\n',
                f'        {pbe0_string} &SCREENING\n',
                f'        {pbe0_string} EPS_SCHWARZ 1.0E-6\n',
                f'        {pbe0_string} SCREEN_ON_INITIAL_P FALSE\n',
                f'        {pbe0_string} &END SCREENING\n',
                f'        {pbe0_string} &INTERACTION_POTENTIAL\n',
                f'        {pbe0_string} POTENTIAL_TYPE TRUNCATED\n',
                f'        {pbe0_string} CUTOFF_RADIUS 2.0\n',
                f'        {pbe0_string} T_C_G_DATA ./t_c_g.dat\n',
                f'        {pbe0_string} &END INTERACTION_POTENTIAL\n',
                f'        {pbe0_string} &MEMORY\n',
                f'        {pbe0_string} MAX_MEMORY  2000\n',
                f'        {pbe0_string} EPS_STORAGE_SCALING 0.1\n',
                f'        {pbe0_string} &END MEMORY\n',
                f'        {pbe0_string} FRACTION 0.25\n',
                f'    {pbe0_string} &END HF\n',
                 '    &END XC\n',
                 '  &END DFT\n',
                 '  STRESS_TENSOR ANALYTICAL \n',
                 '  &SUBSYS    #System description\n',
                 '    &CELL \n',
                 f'      ABC [angstrom] {a} {b} {c}   !Unit cells that are orthorhombic are more efficient with CP2K\n',
                 f'      ALPHA_BETA_GAMMA [deg] {alpha} {beta} {gamma}        !Specify the angles between the vectors A, B and C when using the ABC keyword\n',
                 '      SYMMETRY HEXAGONAL                    !Imposes an initial cell symmetry\n',
                 '    &END CELL\n',
                 '\n',
                 '    #Atom coordinates can be in the &COORD section or provided as an external file on TOPOLOGY.\n',
                 '    &COORD\n']
    for i in range(len(atom_labels)):
        opt_file += [f'      {atom_labels[i]}    {atom_pos[i][0]:>14.10f}  {atom_pos[i][1]:>14.10f}  {atom_pos[i][2]:>14.10f}\n']
        
    opt_file += ['    &END COORD\n',
                 '\n',
                 "    ! MOLOPT basis sets are fairly costly, but in the 'DZVP-MOLOPT-SR-GTH' available for all elements\n",
                 '    ! their contracted nature makes them suitable for condensed and gas phase systems alike.\n']
    if 'H' in atom_types:
        opt_file += ['    &KIND H                              \n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q1\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q1             \n',
                f'     {pbe0_string} BASIS_SET AUX_FIT cFIT3       \n',
                 '    &END KIND\n']
    if 'B' in atom_types:
        opt_file += ['    &KIND B                              \n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q3\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q3             \n',
                f'     {pbe0_string} BASIS_SET AUX_FIT cFIT3       \n',
                 '    &END KIND\n']

    if 'Be' in atom_types:
        opt_file += ['    &KIND Be                              \n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q4\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q4             \n',
                f'     {pbe0_string} BASIS_SET AUX_FIT cFIT3       \n',
                 '    &END KIND\n']
    if 'C' in atom_types:
        opt_file += ['    &KIND C\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q4\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q4\n',
                f'     {pbe0_string} BASIS_SET AUX_FIT cFIT3       \n',
                 '    &END KIND\n']
    if 'N' in atom_types:   
        opt_file += ['    &KIND N\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q5\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q5\n',
                f'     {pbe0_string} BASIS_SET AUX_FIT cFIT3       \n',
                 '    &END KIND\n']
    if 'O' in atom_types:   
        opt_file += ['    &KIND O\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q6\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q6\n',
                f'     {pbe0_string} BASIS_SET AUX_FIT cFIT3       \n',
                 '    &END KIND  \n']
    if 'F' in atom_types:   
        opt_file += ['    &KIND F\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q7\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q7\n',
                f'     {pbe0_string} BASIS_SET AUX_FIT cFIT3       \n',
                 '    &END KIND  \n']
    if 'S' in atom_types:   
        opt_file += ['    &KIND S\n',
                 '      BASIS_SET DZVP-MOLOPT-SR-GTH-q6\n',
                 '      MAGNETIZATION 0\n',
                 '      POTENTIAL GTH-PBE-q6\n',
                f'     {pbe0_string} BASIS_SET AUX_FIT cFIT3       \n',
                 '    &END KIND  \n']
                
    opt_file += ['  &END SUBSYS\n',
                 '&END FORCE_EVAL\n',
                 '!&EXT_RESTART\n',
                 f'!  RESTART_FILE_NAME {prefix}-1.restart\n',
                 '!&END\n']

    file = open(os.path.join(out_path, prefix, 'CHARGES', prefix + '.scf.in'), 'w')
    
    for i in opt_file:
        file.write(i)
    file.close()


def cread_sh_script(path, prefix, opt_n='OPT_1', system_type='PBS'):

    if opt_n.upper() == 'OPT_1':
        sufix = 'cell_opt'
    if opt_n.upper() == 'BOMD':
        sufix = 'bomd'
    if opt_n.upper() == 'OPT_2':
        sufix = 'cell_opt'
    if opt_n.upper() == 'CHARGES':
        sufix = 'scf'

    if system_type == 'PBS':
        opt_file = ['#!/bin/sh\n',
                    '#PBS -l walltime=36:00:00\n',
                    '#PBS -l select=1:ncpus=48:mpiprocs=48\n',
                    f'#PBS -N {prefix}_{opt_n}\n',
                    '#PBS -m bea \n',
                    '#PBS -j oe\n',
                    '#PBS -V\n',
                    '\n',
                    'cd ${PBS_O_WORKDIR}\n']

    if system_type == 'SLURM':
        opt_file = ['#!/bin/sh\n',
                    '#PBS -l walltime=36:00:00\n',
                    '#PBS -l select=1:ncpus=48:mpiprocs=48\n',
                    f'#PBS -N {prefix}_{opt_n}\n',
                    '#PBS -m bea \n',
                    '#PBS -j oe\n',
                    '#PBS -V\n',
                    '\n',
                    'cd ${PBS_O_WORKDIR}\n']

    opt_file += ['\n',
                'module load intel/2019.4\n',
                'module load gcc/7.4.0\n',
                'module load anaconda3/2020.11\n',
                '\n',
                'source /home/users/lipelopes/bin/cp2k-8.2/tools/toolchain/install/setup\n',
                '\n',
                'export I_MPI_PIN_DOMAIN=auto \n',
                'export I_MPI_PIN_ORDER=bunch \n',
                'export OMP_PLACES=threads \n',
                'export OMP_PROC_BIND=SPREAD \n',
                'export OMP_NUM_THREADS=4\n',
                '\n',
                f'mpirun -np 12 /home/users/lipelopes/bin/cp2k-8.2/exe/local/cp2k.psmp -i {prefix}.{sufix}.in -o {prefix}.{sufix}.out\n',
                'rm *.wfn* \n']

    if opt_n == 'CHARGES':
        opt_file += ['\n',
                     f'mv {prefix}-valence_density-ELECTRON_DENSITY-1_0.cube valence_density.cube \n',
                     '\n',
                     'ATOMIC_D_PATH="/home/users/lipelopes/atomic_densities/"\n',
                     'DDCE_P="/home/users/lipelopes/bin/Chargemol_09_26_2017_linux_parallel"\n',
                     '\n',
                     'cat>job_control.txt<<EOF\n',
                     '<net charge>\n',
                     '0.0 <-- specifies the net charge of the unit cell (defaults to 0.0 if nothing specified)\n',
                     '</net charge>\n',
                     '\n',
                     '<periodicity along A, B, and C vectors>\n',
                     '.true.\n',
                     '.true.\n',
                     '.true.\n',
                     '</periodicity along A, B, and C vectors>\n',
                     '\n',
                     '<atomic densities directory complete path>\n',
                     '$ATOMIC_D_PATH\n',
                     '</atomic densities directory complete path>\n',
                     '\n',
                     '<charge type>\n',
                     'DDEC6 <-- specifies the charge type (DDEC3 or DDEC6)\n',
                     '</charge type>\n',
                     '\n',
                     '<compute BOs>\n',
                     '.true. <-- specifies whether to compute bond orders or not\n',
                     '</compute BOs>\n',
                     '\n',
                     '<number of core electrons>\n',
                     '</number of core electrons>\n',
                     'EOF\n',
                     '\n',
                     '$DDCE_P job_control.txt \n',
                     'rm *.cube\n']

    file = open(os.path.join(path, prefix, opt_n, prefix + '.sh'), 'w')

    for i in opt_file:
        file.write(i)
    file.close()


#############################################################################

COF_DATABASE_PATH = '/home/users/lipelopes/COF_Database/'
OUT_PATH = '/scratch/31081a/lipelopes/Teste_Database/'

COF_DATABASE_PATH = '/home/users/lipelopes/Size_matter_Database/'
OUT_PATH = '/scratch/31081a/lipelopes/Size_matter_Database_Calc/'

try:
    os.mkdir(OUT_PATH)
except:
    None

COFs_list = [i.rstrip('.json') for i in os.listdir(COF_DATABASE_PATH) if '.json' in i]

OPT_1_list = []
BOMD_list = []
OPT_2_list = []
CHARGES_LIST = []

for i in COFs_list:
    cof_json = Tools.read_json(COF_DATABASE_PATH, i)
    opt_1, bomd, opt_2, charges = False, False, False, False

    if 'opt_1' in cof_json['system'].keys():
        opt_1 = cof_json['system']["opt_1"]
    if 'bomd' in cof_json['system'].keys():
        bomd = cof_json['system']['bomd']
    if 'opt_2' in cof_json['system'].keys():
        opt_2 = cof_json['system']['opt_2']
    if 'DDEC_charges' in cof_json['system'].keys():
        charges = cof_json['system']['DDEC_charges']

    if opt_1 is False:
        OPT_1_list += [i]
    if opt_1 is True and bomd is False:
        BOMD_list += [i]
    if opt_1 is True and bomd is True and opt_2 is False:
        OPT_2_list += [i]
    if opt_1 is True and bomd is True and opt_2 is True and charges == False:
        CHARGES_LIST += [i]

job_type = 'opt_full' # opt1, bomd, opt2, opt_full

n_run, rodando, fila = get_available_runs()

max_submits = n_run
submit = 0

print('OPT_1', len(OPT_1_list))
print('BOMD', len(BOMD_list))
print('OPT_2', len(OPT_2_list))
print('CHARGES', len(CHARGES_LIST))

if len(OPT_1_list) > 0:
    for cof in OPT_1_list:
        if submit < max_submits:
            cread_opt_input(COF_DATABASE_PATH, OUT_PATH, cof, opt_n='opt_1', kpoints=False, smearing=False)
            cread_sh_script(OUT_PATH, cof, opt_n='OPT_1')
            os.chdir(os.path.join(OUT_PATH, cof, 'OPT_1'))
            os.system(f'qsub {cof}.sh')
    
            submit += 1
            os.chdir(OUT_PATH)

if len(BOMD_list) > 0:
    for cof in BOMD_list:
        if submit < max_submits:
            cread_bomd_input(COF_DATABASE_PATH, OUT_PATH, cof, kpoints=False, smearing=False)
            cread_sh_script(OUT_PATH, cof, opt_n='BOMD')
            os.chdir(os.path.join(OUT_PATH, cof, 'BOMD'))
            os.system(f'qsub {cof}.sh')
            submit += 1
            os.chdir(OUT_PATH)

if len(OPT_2_list) > 0:
    for cof in OPT_2_list:
        if submit < max_submits:
            cread_opt_input(COF_DATABASE_PATH, OUT_PATH, cof, opt_n='opt_2', kpoints=False, smearing=False)
            cread_sh_script(OUT_PATH, cof, opt_n='OPT_2')
            os.chdir(os.path.join(OUT_PATH, cof, 'OPT_2'))
            os.system(f'qsub {cof}.sh')
            submit += 1
            os.chdir(OUT_PATH)

if len(CHARGES_LIST) > 0:
    for cof in CHARGES_LIST:
        if submit < max_submits:
            cread_charges_input(COF_DATABASE_PATH, OUT_PATH, cof, kpoints=False, smearing=False)
            cread_sh_script(OUT_PATH, cof, opt_n='CHARGES')
            os.chdir(os.path.join(OUT_PATH, cof, 'CHARGES'))
            os.system(f'qsub {cof}.sh')
            submit += 1
            os.chdir(OUT_PATH)
