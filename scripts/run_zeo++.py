# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020
@author: lipel
"""
import tools as Tools
import os
from tqdm import tqdm

ZEO_CMD = '/home/felipelopes/zeo++-0.3/network'

COF_DATABASE_PATH = '/home/felipelopes/Size_M/concluido'

COF_DATABASE_PATH = '/home/felipelopes/concluido/results/'
COF_PATH = '/home/felipelopes/concluido'

########################
ZEO_CMD = '/home/felipe/Programas/zeo++-0.3/network'
COF_DATABASE_PATH = '/home/felipe/Simulations/Size_Matter/'
COF_PATH = '/home/felipe/Simulations/Size_EQeq'
#######################

COF_list = [i for i in os.listdir(COF_DATABASE_PATH) if '.json' in i]

ROOT_PATH = os.getcwd()

for cof in tqdm(COF_list):
    cof_name = cof.rstrip('.json')
    cof_path = os.path.join(COF_PATH, cof_name)
    cof_data = Tools.read_json(COF_DATABASE_PATH, cof_name)

    if 'geo_opt' in cof_data['system'].keys() and "textural" not in cof_data['system'].keys():
        try:
            os.mkdir(cof_path)
        except:
            None
        try:
            os.mkdir(os.path.join(cof_path, 'ZEO'))
        except:
            None

        Tools.save_cif(os.path.join(cof_path, 'ZEO'), 
                       cof_name, cof_data['geometry']['cell_matrix'], 
                       cof_data['geometry']['atom_labels'], 
                       cof_data['geometry']['atom_pos'], 
                       partial_charges=False, 
                       frac_coords=False)

        os.chdir(os.path.join(cof_path, 'ZEO'))
        
        #Calculate the largest included sphere, free sphere and included sphere along free sphere path
        #os.system(f'{ZEO_CMD} -ha -res large_sphere_includ.res {cof_name}.cif >/dev/null')

        #Chanel identification and dimensionality - 1.4 is the vdW radius of He
        if os.path.exists(f'{cof_name}.chan') is False: 
            os.system(f'{ZEO_CMD} -ha -chan 1.4 {cof_name}.cif >/dev/null')

        #Acessible surface area - 1.82 is the vdW radius of N2
        if os.path.exists(f'{cof_name}.sa') is False: 
            os.system(f'{ZEO_CMD} -ha -sa 1.82 1.82 2000 {cof_name}.cif >/dev/null')

        #Acessible volume - 1.4 is the vdW radius of He
        if os.path.exists(f'{cof_name}.vol') is False: 
            os.system(f'{ZEO_CMD} -ha -vol 1.4 1.4 2000 {cof_name}.cif >/dev/null')

        #Pore size distribution
        #if os.path.exists(f'{cof_name}.psd_histo') is False: 
        #    os.system(f'{ZEO_CMD} -ha -psd 1.82 1.82 5000 {cof_name}.cif &>/dev/null')

        os.chdir(ROOT_PATH)

    
