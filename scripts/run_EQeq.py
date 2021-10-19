# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020
@author: lipel
"""
import tools as Tools
import os
from tqdm import tqdm

EQEQ_PATH = '/home/felipelopes/eqeq/'

COF_DATABASE_PATH = '/home/felipelopes/Size_M/Size_matter/concluido'
OUT_PATH = '/home/felipelopes/Size_M/Size_matter/'


COF_DATABASE_PATH = '/home/felipelopes/concluido/results/'

OUT_PATH = '/home/felipelopes/concluido'


EQEQ_PATH = '/home/felipe/Programas/EQeq/'

COF_DATABASE_PATH = '/home/felipe/Simulations/Size_Matter/'

OUT_PATH = '/home/felipe/Simulations/Size_EQeq'

COF_list = [i for i in os.listdir(COF_DATABASE_PATH) if '.json' in i]

for cof in tqdm(COF_list):
    cof_name = cof.rstrip('.json')
    cof_path = os.path.join(OUT_PATH, cof_name)
    cof_data = Tools.read_json(COF_DATABASE_PATH, cof_name)

    if 'geo_opt' in cof_data['system'].keys() and cof_data['system']['geo_opt'] is True:
        try:
            os.mkdir(cof_path)
        except:
            None
        try:
            os.mkdir(os.path.join(cof_path, 'EQeq'))
        except:
            None

        Tools.save_cif(os.path.join(cof_path, 'EQeq'), 
                       cof_name, 
                       cof_data['geometry']['cell_matrix'], 
                       cof_data['geometry']['atom_labels'], 
                       cof_data['geometry']['atom_pos'], 
                       partial_charges=False, 
                       frac_coords=False)

        os.chdir(os.path.join(cof_path, 'EQeq'))
        
        if os.path.exists('EQeq.log') is False: 
            os.system(f'cp {EQEQ_PATH}data/ionizationdata.dat .')
            os.system(f'cp {EQEQ_PATH}data/chargecenters.dat .')
            os.system(f'{EQEQ_PATH}eqeq {cof_name}.cif chargePrecision=6 >/dev/null | tee EQeq.log')

        os.chdir(OUT_PATH)

    
