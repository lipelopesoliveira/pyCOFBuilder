# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020
@author: lipel
"""
import simplejson
import os
import tools as Tools
import numpy as np
from tqdm import tqdm
import traceback

from get_CP2K_data import *
from get_RASPA_data import *
from get_zeo_data import *

########################### Add Data to JSON ##########################

def add_opt_geo(source_path, databank_path, name, save_opt_history=False, symetrize=False):
    
    COF_json = Tools.read_json(databank_path, name)

    if cp2k_concluded(os.path.join(source_path, name, 'OPT_1', name + '.cell_opt.out')) is not False:
        COF_json['optimization']['opt1'] = get_CP2k_results(source_path, name, 'OPT_1', save_opt_history)
        COF_json['system']['opt_1'] = True

        cp2k_results = get_cp2k_opt_results(source_path, name, 'OPT_1', symetrize)

        if cp2k_results is not False:
            cell, opt_labels, opt_pos, fermi, gap, energy_list, density, composition, symm_info = cp2k_results

            cell_par = Tools.cell_to_cellpar(np.array(cell)).tolist()
            cell_par =  [round(i, 10) for i in cell_par]

            COF_json['geometry']['cell_matrix'] = cell
            COF_json['geometry']['cell_parameters'] = cell_par
            COF_json['geometry']['atom_labels'] = opt_labels
            COF_json['geometry']['atom_pos'] = opt_pos
            COF_json['geometry']['energy'] = energy_list[-1]
            COF_json['geometry']['band_gap'] = gap
            COF_json['geometry']['fermi_level'] = fermi
            COF_json['geometry']['density'] = density
            COF_json['geometry']['composition'] = composition

            if symm_info is not None:
                lattice_type, space_group, space_group_n, symm_op, hall = symm_info

                COF_json['geometry']['lattice_type'] = lattice_type
                COF_json['geometry']['space_group'] = space_group
                COF_json['geometry']['space_group_number'] = space_group_n
                COF_json['geometry']['symm_op'] = symm_op
                COF_json['geometry']['hall'] = hall
    
    if cp2k_concluded(os.path.join(source_path, name, 'BOMD', name + '.bomd.out')) is not False:
        COF_json['optimization']['bomd'] = get_CP2k_results(source_path, name, 'BOMD', save_opt_history)
        COF_json['system']['bomd'] = True

        cp2k_results = get_cp2k_opt_results(source_path, name, 'BOMD', symetrize)

        if cp2k_results is not False:
            cell, opt_labels, opt_pos, fermi, gap, energy_list, density, composition, symm_info = cp2k_results

            cell_par = Tools.cell_to_cellpar(np.array(cell)).tolist()
            cell_par =  [round(i, 10) for i in cell_par]

            COF_json['geometry']['cell_matrix'] = cell
            COF_json['geometry']['cell_parameters'] = cell_par
            COF_json['geometry']['atom_labels'] = opt_labels
            COF_json['geometry']['atom_pos'] = opt_pos
            COF_json['geometry']['energy'] = energy_list[-1]
            COF_json['geometry']['band_gap'] = gap
            COF_json['geometry']['fermi_level'] = fermi
            COF_json['geometry']['density'] = density
            COF_json['geometry']['composition'] = composition

            if symm_info is not None:
                lattice_type, space_group, space_group_n, symm_op, hall = symm_info

                COF_json['geometry']['lattice_type'] = lattice_type
                COF_json['geometry']['space_group'] = space_group
                COF_json['geometry']['space_group_number'] = space_group_n
                COF_json['geometry']['symm_op'] = symm_op
                COF_json['geometry']['hall'] = hall

    if cp2k_concluded(os.path.join(source_path, name, 'OPT_2', name + '.cell_opt.out')) is not False:
        COF_json['optimization']['opt2'] = get_CP2k_results(source_path, name, 'OPT_2', save_opt_history)
        COF_json['system']['opt_2'] = True
        COF_json['system']['geo_opt'] = True

        cp2k_results = get_cp2k_opt_results(source_path, name, 'OPT_2', symetrize)

        if cp2k_results is not False:
            cell, opt_labels, opt_pos, fermi, gap, energy_list, density, composition, symm_info = cp2k_results

            cell_par = Tools.cell_to_cellpar(np.array(cell)).tolist()
            cell_par =  [round(i, 10) for i in cell_par]

            COF_json['geometry']['cell_matrix'] = cell
            COF_json['geometry']['cell_parameters'] = cell_par
            COF_json['geometry']['atom_labels'] = opt_labels
            COF_json['geometry']['atom_pos'] = opt_pos
            COF_json['geometry']['energy'] = energy_list[-1]
            COF_json['geometry']['band_gap'] = gap
            COF_json['geometry']['fermi_level'] = fermi
            COF_json['geometry']['density'] = density
            COF_json['geometry']['composition'] = composition

            if symm_info is not None:
                lattice_type, space_group, space_group_n, symm_op, hall = symm_info

                COF_json['geometry']['lattice_type'] = lattice_type
                COF_json['geometry']['space_group'] = space_group
                COF_json['geometry']['space_group_number'] = space_group_n
                COF_json['geometry']['symm_op'] = symm_op
                COF_json['geometry']['hall'] = hall

    Tools.write_json(databank_path, name, COF_json)

def add_DDEC_charges(source_path, databank_path, cof_name):

    if DDEC_charges_concluded(source_path, cof_name):
        COF_json = Tools.read_json(databank_path, cof_name)

        charges, bond_orders = read_DDEC_charges(source_path, cof_name)
        
        COF_json['system']['DDEC_charges'] = True
        COF_json['geometry']['DDEC_charges'] = charges
        COF_json['geometry']['bond_orders'] = bond_orders
        
        Tools.write_json(databank_path, cof_name, COF_json)

def add_PHONON_results(source_path, databank_path, cof_name):

    if os.path.exists(os.path.join(source_path, cof_name + '-VIBRATIONS-1.mol')):

        freq, intensity = read_VIB_data(source_path, cof_name)

        COF_json = Tools.read_json(databank_path, cof_name)
        
        COF_json['system']['phonon'] = True
        COF_json['spectrum']['FTIR'] = {'frequency_cm-1':freq, 'IR_intensity':intensity}
        
        Tools.write_json(databank_path, cof_name, COF_json)

def add_charges(source_path, databank_path, cof_name, pbe0=False):

    if DDEC_charges_concluded(source_path, cof_name):
        COF_json = Tools.read_json(databank_path, cof_name)

        lowdin_charges, hirshfeld_charges, mulliken_charges = read_charges(source_path, cof_name)

        COF_json['system']['lowdin_charges'] = True
        COF_json['system']['hirshfeld_charges'] = True
        COF_json['system']['mulliken_charges'] = True
        COF_json['geometry']['lowdin_charges'] = lowdin_charges
        COF_json['geometry']['hirshfeld_charges'] = hirshfeld_charges
        COF_json['geometry']['mulliken_charges'] = mulliken_charges

        if pbe0:
            fermi, gap, t_energy = get_gap(source_path, cof_name, 'CHARGES')
            COF_json['geometry']['band_gap_PBE0'] = gap
            COF_json['geometry']['fermi_level_PBE0'] = fermi
            COF_json['geometry']['energy_PBE0'] = t_energy[0]

        
        Tools.write_json(databank_path, cof_name, COF_json)
    
def add_RASPA_data(source_path, databank_path, cof_name, gas='CO2', temp='298'):
    ADS_dir = os.path.join(source_path, cof_name, 'RASPA')
    
    if os.path.exists(ADS_dir):
        ads_dict = read_raspa_single_ads(ADS_dir, gas, temp)
        if ads_dict is None:
            print(cof_name, 'ERRO!!')
        else:
            COF_json = Tools.read_json(databank_path, cof_name)
            if gas in COF_json['adsorption'].keys():
                COF_json['adsorption'][gas].update({round(ads_dict['temperature'], 2):ads_dict})
            else:
                COF_json['adsorption'][gas] = {ads_dict['temperature']:ads_dict}

            Tools.write_json(databank_path, cof_name, COF_json)


def add_ZEO_data(source_path, databank_path, cof_name):

    ZEO_PATH = os.path.join(source_path, cof_name, 'ZEO')

    if os.path.exists(ZEO_PATH):

        cof_data = Tools.read_json(databank_path, cof_name)
  
        zeo_data = get_zeo_results(ZEO_PATH, cof_name)
    
        if len(zeo_data.keys()) > 1:
            cof_data['system']['textural'] = True
            cof_data['textural'].update(zeo_data)

            Tools.write_json(databank_path, cof_name, cof_data)

def add_EQeq_data(source_path, databank_path, cof_name):

    EQeq_PATH = os.path.join(source_path, cof_name, 'EQeq')

    if os.path.exists(EQeq_PATH):

        cof_data = Tools.read_json(databank_path, cof_name)
  
        EQeq_data = open(os.path.join(EQeq_PATH, cof_name + '.cif_EQeq_ewald_0.00_-2.00.json'), 'r').readlines()[0]
        EQeq_data = [float(i) for i in EQeq_data.lstrip('[').rstrip(']\n').split(',')]

        cof_data['system']['eqeq_charges'] = True
        cof_data['geometry']['eqeq_charges'] = EQeq_data

        Tools.write_json(databank_path, cof_name, cof_data)

if __name__ == '__main__':
    #############################################################################################

    COF_DATABASE_PATH = '/home/felipelopes/concluido/COF_Database'
    OUT_PATH = '/home/felipelopes/concluido/'

    ######################################
    COF_DATABASE_PATH = '/home/felipe/Simulations/Size_Matter/'
    OUT_PATH = '/home/felipe/Simulations/Size_EQeq'
    ######################################

    if os.path.exists(COF_DATABASE_PATH) is False:
        print(COF_DATABASE_PATH, 'created succesfully!')
        os.mkdir(COF_DATABASE_PATH)

    COFs = [i for i in os.listdir(OUT_PATH) if os.path.isdir(os.path.join(OUT_PATH, i))]

    folders_to_remove = ['original', '__pycache__', 'rodando', 'concluido', '.vscode', 'results', COF_DATABASE_PATH.split('/')[-1]]

    for r in folders_to_remove:
        try:
            COFs.remove(r)
        except:
            None

    for cof in tqdm(COFs):
        try:
            if os.path.exists(os.path.join(COF_DATABASE_PATH, f'{cof}.json')) is False:
                Tools.write_json(COF_DATABASE_PATH, cof, Tools.create_COF_json(cof))
            add_opt_geo(OUT_PATH, COF_DATABASE_PATH, cof, save_opt_history=False)
            add_DDEC_charges(OUT_PATH, COF_DATABASE_PATH, cof)
            add_charges(OUT_PATH, COF_DATABASE_PATH, cof, pbe0=True)
            for g in ['CO2']:
                add_RASPA_data(OUT_PATH, COF_DATABASE_PATH, cof, g)
            add_ZEO_data(OUT_PATH, COF_DATABASE_PATH, cof)
            add_EQeq_data(OUT_PATH, COF_DATABASE_PATH, cof)
        except Exception:
            print('Erro:', cof)
            traceback.print_exc()