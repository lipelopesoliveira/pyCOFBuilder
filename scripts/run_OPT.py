# -*- coding: utf-8 -*-
"""
Created on Sat Oct 16 13:49:50 2021
@author: Felipe Lopes de Oliveira
"""

import os
import numpy as np
import tools as Tools
from tqdm import tqdm
from textwrap import dedent
from create_xTB_input import *
from get_results import *
import time

def create_sh_script(CMD_SOURCE, CMD_CP2K, name, opt_type):
    ssh_script = dedent(f"""
    {CMD_SOURCE}
    {CMD_CP2K} -i {name}.{opt_type}.in -o {name}.{opt_type}.out""")

    file = open('run.sh', 'w')
    
    for i in ssh_script:
        file.write(i)
    file.close()

def concluded(out_file):

    if os.path.exists(out_file):
        concluded = False
        for line in open(out_file, 'r').readlines():
            if 'PROGRAM ENDED AT ' in line:
                concluded = True
        return concluded
    else:
        return False

###############          DEFINE THE PATH FOR DATABASE AND OUTPUT FILES          ###############

COF_DATABASE_PATH = '/home/felipe/Simulations/Teste_Optmization/Data_Base'
OUT_PATH = '/home/felipe/Simulations/Teste_Optmization/Data_Out'

###############               DEFINE THE COMMANDS FOR EXECUTABLES              ###############

CMD_SOURCE = 'source /home/felipe/Programas/cp2k/tools/toolchain/install/setup'
CMD_CP2K = 'mpirun -np 2 /home/felipe/Programas/cp2k/exe/local/cp2k.popt'
EQEQ_PATH = '/home/felipe/Programas/EQeq/'
ZEO_CMD = '/home/felipe/Programas/zeo++-0.3/network'

K_POINTS = False

###############            LOOP OVER THE STRUCTURES AND OPTMIZE ALL            ###############

COF_list = [i for i in os.listdir(COF_DATABASE_PATH) if '.json' in i]

ROOT_PATH = os.getcwd()

for cof in tqdm(COF_list):
    COF_NAME = cof.rstrip('.json')
    COF_PATH = os.path.join(OUT_PATH, COF_NAME)

    COF_JSON = Tools.read_json(COF_DATABASE_PATH, COF_NAME)

    if 'opt_1' not in list(COF_JSON['system'].keys()):

        start_time = time.time()

        create_opt_input(COF_DATABASE_PATH, OUT_PATH, COF_NAME, opt_n='opt_1', kpoints=K_POINTS, smearing=False)
        os.chdir(os.path.join(COF_PATH, 'OPT_1'))
        create_sh_script(CMD_SOURCE, CMD_CP2K, COF_NAME, 'cell_opt')
        os.system('bash run.sh >/dev/null')

        if concluded(COF_NAME + '.cell_opt.out'):
            os.system('rm *.restart *.wfn* ')
            os.chdir(ROOT_PATH)
            add_opt_geo(OUT_PATH, COF_DATABASE_PATH, COF_NAME, save_opt_history=False)

            opt1_time = time.time()

            create_bomd_input(COF_DATABASE_PATH, OUT_PATH, COF_NAME, kpoints=K_POINTS, smearing=False)
            os.chdir(os.path.join(COF_PATH, 'BOMD'))
            create_sh_script(CMD_SOURCE, CMD_CP2K, COF_NAME, 'bomd')
            os.system('bash run.sh >/dev/null')
        else:
            print(f'ERROR ON OPT_1 of {COF_NAME}')
            continue

        if concluded(COF_NAME + '.bomd.out'):

            os.system('rm *.restart *.wfn* ')
            os.chdir(ROOT_PATH)
            add_opt_geo(OUT_PATH, COF_DATABASE_PATH, COF_NAME, save_opt_history=False)

            bomd_time = time.time()

            create_opt_input(COF_DATABASE_PATH, OUT_PATH, COF_NAME, opt_n='opt_2', kpoints=K_POINTS, smearing=True)
            os.chdir(os.path.join(COF_PATH, 'OPT_2'))
            create_sh_script(CMD_SOURCE, CMD_CP2K, COF_NAME, 'cell_opt')
            os.system('bash run.sh >/dev/null')
        else:
            print(f'ERROR ON BOMD of {COF_NAME}')
            continue

        if concluded(COF_NAME + '.cell_opt.out'):
            os.system('rm *.restart *.wfn* ')
            os.chdir(ROOT_PATH)
            add_opt_geo(OUT_PATH, COF_DATABASE_PATH, COF_NAME, save_opt_history=False)

            opt2_time = time.time()

            create_vib_input(COF_DATABASE_PATH, OUT_PATH, COF_NAME, kpoints=K_POINTS, smearing=True)
            os.chdir(os.path.join(COF_PATH, 'PHONON'))
            create_sh_script(CMD_SOURCE, CMD_CP2K, COF_NAME, 'ph')
            os.system('bash run.sh >/dev/null')

            add_PHONON_results(os.path.join(COF_PATH, 'PHONON'), COF_DATABASE_PATH, COF_NAME)
            phonon_time = time.time()

            try:
                os.mkdir(os.path.join(COF_PATH, 'EQeq'))
            except:
                None

            COF_JSON = Tools.read_json(COF_DATABASE_PATH, COF_NAME)

            Tools.save_cif(os.path.join(COF_PATH, 'EQeq'), 
                        COF_NAME, 
                        COF_JSON['geometry']['cell_matrix'], 
                        COF_JSON['geometry']['atom_labels'], 
                        COF_JSON['geometry']['atom_pos'], 
                        partial_charges=False, 
                        frac_coords=False)

            os.chdir(os.path.join(COF_PATH, 'EQeq'))
            
            if os.path.exists('EQeq.log') is False: 
                os.system(f'cp {EQEQ_PATH}data/ionizationdata.dat .')
                os.system(f'cp {EQEQ_PATH}data/chargecenters.dat .')
                os.system(f'{EQEQ_PATH}eqeq {COF_NAME}.cif chargePrecision=6 2>EQeq.log')

            add_EQeq_data(OUT_PATH, COF_DATABASE_PATH, COF_NAME)

            eqeq_time = time.time()

            try:
                os.mkdir(os.path.join(COF_PATH, 'ZEO'))
            except:
                None

            Tools.save_cif(os.path.join(COF_PATH, 'ZEO'), 
                        COF_NAME, COF_JSON['geometry']['cell_matrix'], 
                        COF_JSON['geometry']['atom_labels'], 
                        COF_JSON['geometry']['atom_pos'], 
                        partial_charges=False, 
                        frac_coords=False)

            os.chdir(os.path.join(COF_PATH, 'ZEO'))
            
            #Calculate the largest included sphere, free sphere and included sphere along free sphere path
            #os.system(f'{ZEO_CMD} -ha -res large_sphere_includ.res {cof_name}.cif >/dev/null')

            #Chanel identification and dimensionality - 1.4 is the vdW radius of He
            if os.path.exists(f'{COF_NAME}.chan') is False: 
                os.system(f'{ZEO_CMD} -ha -chan 1.4 {COF_NAME}.cif >/dev/null')

            #Acessible surface area - 1.82 is the vdW radius of N2
            if os.path.exists(f'{COF_NAME}.sa') is False: 
                os.system(f'{ZEO_CMD} -ha -sa 1.82 1.82 2000 {COF_NAME}.cif >/dev/null')

            #Acessible volume - 1.4 is the vdW radius of He
            if os.path.exists(f'{COF_NAME}.vol') is False: 
                os.system(f'{ZEO_CMD} -ha -vol 1.4 1.4 2000 {COF_NAME}.cif >/dev/null')

            #Pore size distribution
            #if os.path.exists(f'{cof_name}.psd_histo') is False: 
            #    os.system(f'{ZEO_CMD} -ha -psd 1.82 1.82 5000 {cof_name}.cif &>/dev/null')

            add_ZEO_data(OUT_PATH, COF_DATABASE_PATH, COF_NAME)

            zepp_time = time.time()

            os.chdir(OUT_PATH)

            COF_JSON = Tools.read_json(COF_DATABASE_PATH, COF_NAME)

            times = {'opt_1': opt1_time - start_time,
                    'bomd': bomd_time - opt1_time,
                    'opt_2': opt2_time - bomd_time,
                    'phonon': phonon_time - opt2_time,
                    'eqeq': eqeq_time - phonon_time,
                    'zeo++': zepp_time - eqeq_time}
            COF_JSON['system']['execution_times_seconds'] = times

            Tools.write_json(COF_DATABASE_PATH, COF_NAME, COF_JSON)
        else:
            print(f'ERROR ON OPT_2 of {COF_NAME}')