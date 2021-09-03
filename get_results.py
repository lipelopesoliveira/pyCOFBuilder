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

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

COF_PATH= os.getcwd()

COFs = [i for i in os.listdir(COF_PATH) if os.path.isdir(i)]

folders_to_remove = ['original', '__pycache__', 'rodando', 'concluido', '.vscode']

for r in folders_to_remove:
    try:
        COFs.remove(r)
    except:
        None
    
########################### JSON related ##########################  


def save_json(path, name, COF_json):

    if os.path.exists(os.path.join(path, 'concluido')) is not True:
        os.mkdir(os.path.join(path, 'concluido'))
    
    save_path = os.path.join(path, 'concluido', name + '.json')
    
    with open(save_path, 'w', encoding='utf-8') as f:
        simplejson.dump(COF_json, f, ensure_ascii=False, separators=(',', ':'), indent=2, ignore_nan=True)
        
def read_json(path, cof_name):
    
    cof_path = os.path.join(path, 'concluido', cof_name + '.json')
    
    with open(cof_path, 'r') as r:
        json_object = simplejson.loads(r.read())
    
    return json_object


########################### CP2K related ##########################

def opt_concluded(cof_name, opt_type='OPT_1'):
    
    out_path = os.path.join(COF_PATH, cof_name, opt_type, cof_name + '.cell_opt.out')
    
    if os.path.exists(out_path):
    
        opt_out = open(out_path, 'r').readlines()
        
        concluded = False
        for line in opt_out:
            if 'PROGRAM ENDED AT' in line:
                concluded = True
        return concluded
    else:
        return False

def cp2k_concluded(out_file):

    if os.path.exists(out_file):
        
        opt_out = open(out_file, 'r').readlines()
        
        concluded = False
        for line in opt_out:
            if 'PROGRAM ENDED AT ' in line:
                concluded = True
        return concluded
    else:
        return False

def charges_concluded(cof_name):
    
    out_path = os.path.join(COF_PATH, cof_name, 'CHARGES', 'valence_cube_DDEC_analysis.output')
    
    if os.path.exists(out_path):
    
        opt_out = open(out_path, 'r').readlines()
        
        concluded = False
        for line in opt_out:
            if 'Finished chargemol' in line:
                concluded = True
        return concluded
    else:
        return False

def get_optimized_geometry(cof_name, opt_type='OPT_1'):
    
    xyz_path = os.path.join(COF_PATH, cof_name, opt_type, cof_name + '-pos-1.xyz')
    xyz_opt = open(xyz_path, 'r').readlines()
    
    n_atoms = int(xyz_opt[0].rstrip('\n'))
    opt = [i.split() for i in xyz_opt[-n_atoms:]]
    
    opt_labels = []
    opt_pos = []
    
    for i in opt:
        opt_labels += [i[0]]
        opt_pos += [[float(j) for j in i[1:]]]
        
    return opt_labels, opt_pos

def get_optimized_cell_parameters(cof_name, opt_type='OPT_1'):
    
    cell_path = os.path.join(COF_PATH, cof_name, opt_type, cof_name + '-1.cell')
    cell_opt = open(cell_path, 'r').readlines()
    cell = [float(i) for i in cell_opt[-1].split()[2:-1]]
    cell = [cell[:3], cell[3:6], cell[6:]]
    return cell 

def get_all_structures_CP2k(name, save_opt_history=False):
    
    pos_1_file = os.path.join(COF_PATH, name, 'OPT_1', name + '-pos-1.xyz')
    pos_2_file = os.path.join(COF_PATH, name, 'OPT_2', name + '-pos-1.xyz')
    pos_bomd_file = os.path.join(COF_PATH, name, 'BOMD', name + '-pos-1.xyz')

    cell_1_file = os.path.join(COF_PATH, name, 'OPT_1', name + '-1.cell')
    cell_bomd_file = os.path.join(COF_PATH, name, 'BOMD', name + '-1.cell')
    cell_2_file = os.path.join(COF_PATH, name, 'OPT_2', name + '-1.cell')

    opt_1 = {}
    bomd = {}
    opt_2 = {}

    if cp2k_concluded(os.path.join(COF_PATH, name, 'OPT_1', name + '.cell_opt.out')):

        if save_opt_history:

            tmp = open(pos_1_file).readlines()
            n_atoms = int(tmp[0])
            n_steps = [i for i in tmp if 'i' in i]

            for i in range(len(n_steps)):
                start = i*(n_atoms + 2)+2
                end =i*(n_atoms + 2) + 2 + n_atoms
                step = tmp[start -1].split()[2].rstrip(',')
                opt_1[step] = {
                                'cell_matrix':[],
                                'atom_label':[],
                                'atom_pos':[],
                                'energy':[]}

                tmp_pos = tmp[start: end]
                energy = float(tmp[start -1].split()[-1])
                opt_1[step]['energy'] = energy

                tmp_atom_labels, tmp_atom_pos = [], []
                for j in range(len(tmp_pos)):
                    tmp_atom_labels += [tmp_pos[j].split()[0]]
                    tmp_atom_pos += [[float(k) for k in tmp_pos[j].split()[1:]]]

                opt_1[step]['atom_label'] = tmp_atom_labels
                opt_1[step]['atom_pos'] = tmp_atom_pos

            tmp_cell = open(cell_1_file).readlines()[1:]
            for t in range(len(tmp_cell)):
                step = tmp_cell[t].split()[0]
                cell = np.reshape(np.array([float(i) for i in tmp_cell[t].split()[2:-1]]), (3,3)).tolist()
                opt_1[step]['cell_matrix'] = cell
        
        if save_opt_history is False:
            tmp = open(pos_1_file).readlines()
            opt_1['energy'] = [float(i.split()[-1]) for i in tmp if 'i' in i]

    if cp2k_concluded(os.path.join(COF_PATH, name, 'BOMD', name + '.bomd.out')):

        if save_opt_history:

            tmp = open(pos_bomd_file).readlines()
            n_atoms = int(tmp[0])
            n_steps = [i for i in tmp if 'i' in i]

            for i in range(len(n_steps)):

                start = i*(n_atoms + 2)+2
                end =i*(n_atoms + 2) + 2 + n_atoms

                step = tmp[start -1].split()[2].rstrip(',')
                bomd[step] = {
                                'cell_matrix':[],
                                'atom_label':[],
                                'atom_pos':[],
                                'energy':[]}

                tmp_pos = tmp[start: end]

                energy = float(tmp[start -1].split()[-1])
                bomd[step]['energy'] = energy

                tmp_atom_labels, tmp_atom_pos = [], []
                for j in range(len(tmp_pos)):
                    tmp_atom_labels += [tmp_pos[j].split()[0]]
                    tmp_atom_pos += [[float(k) for k in tmp_pos[j].split()[1:]]]

                bomd[step]['atom_label'] = tmp_atom_labels
                bomd[step]['atom_pos'] = tmp_atom_pos

            tmp_cell = open(cell_bomd_file).readlines()[1:]
            for t in range(len(tmp_cell)):
                step = tmp_cell[t].split()[0]
                cell = np.reshape(np.array([float(i) for i in tmp_cell[t].split()[2:-1]]), (3,3)).tolist()
                bomd[step]['cell_matrix'] = cell
        
        if save_opt_history is False:
            tmp = open(pos_bomd_file).readlines()
            bomd['energy'] = [float(i.split()[-1]) for i in tmp if 'i' in i]
    
    if cp2k_concluded(os.path.join(COF_PATH, name, 'OPT_2', name + '.cell_opt.out')):

        if save_opt_history:

            tmp = open(pos_2_file).readlines()
            n_atoms = int(tmp[0])

            n_steps = [i for i in tmp if 'i' in i]

            for i in range(len(n_steps)):

                start = i*(n_atoms + 2)+2
                end =i*(n_atoms + 2) + 2 + n_atoms

                step = tmp[start -1].split()[2].rstrip(',')
                opt_2[step] = {
                                'cell_matrix':[],
                                'atom_label':[],
                                'atom_pos':[],
                                'energy':[]}

                tmp_pos = tmp[start: end]

                energy = float(tmp[start -1].split()[-1])
                opt_2[step]['energy'] = energy

                tmp_atom_labels, tmp_atom_pos = [], []
                for j in range(len(tmp_pos)):
                    tmp_atom_labels += [tmp_pos[j].split()[0]]
                    tmp_atom_pos += [[float(k) for k in tmp_pos[j].split()[1:]]]

                opt_2[step]['atom_label'] = tmp_atom_labels
                opt_2[step]['atom_pos'] = tmp_atom_pos

            tmp_cell = open(cell_2_file).readlines()[1:]
            for t in range(len(tmp_cell)):
                step = tmp_cell[t].split()[0]
                cell = np.reshape(np.array([float(i) for i in tmp_cell[t].split()[2:-1]]), (3,3)).tolist()
                opt_2[step]['cell_matrix'] = cell

        if save_opt_history is False:
            tmp = open(pos_2_file).readlines()
            opt_2['energy'] = [float(i.split()[-1]) for i in tmp if 'i' in i]

    return opt_1, bomd, opt_2
    
def get_charges(cof_name):
    
    charges_path = os.path.join(COF_PATH, cof_name, 'CHARGES', 'DDEC6_even_tempered_net_atomic_charges.xyz')
    
    if os.path.exists(charges_path):
        tmp = open(charges_path, 'r').readlines()
        n_atoms = int(tmp[0].rstrip('\n'))
        
        atom_list = tmp[2: 2+n_atoms]
        charges = []
        
        for i in atom_list:
            charges += [float(i.split()[-1])]
            
        return charges
    else:
        return False

def get_bond_order(cof_name):
    
    charges_path = os.path.join(COF_PATH, cof_name, 'CHARGES', 'DDEC6_even_tempered_bond_orders.xyz')
    
    if os.path.exists(charges_path):
        tmp = open(charges_path, 'r').readlines()
        n_atoms = int(tmp[0].rstrip('\n'))
        
        atom_list = tmp[2: 2+n_atoms]
        bo = []
        
        for i in atom_list:
            bo += [float(i.split()[-1])]
            
        return bo
    else:
        return False
        
def get_gap(cof_name, opt_type='OPT_1'):
    
    out_path = os.path.join(COF_PATH, cof_name, opt_type, cof_name + '.cell_opt.out')
    fermi = 0
    gap = 0
    t_energy = []
    
    if cp2k_concluded(out_path) == True:
        opt_out = open(out_path, 'r').readlines()
        for i in opt_out:
            if 'Fermi Energy [eV]' in i:
                fermi += float(i.split()[-1])
            if 'HOMO - LUMO gap [eV] :' in i:
                gap += float(i.split()[-1])
            if 'ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: ' in i:
                t_energy += [float(i.split()[-1])]
        
    return fermi, gap, t_energy
        
def get_cp2k_opt_results(cof_name, opt_type='OPT_1', symetrize=False):

    out_path = os.path.join(COF_PATH, cof_name, opt_type, cof_name + '.cell_opt.out')
    
    if cp2k_concluded(out_path) == True:
        opt_labels, opt_pos = get_optimized_geometry(cof_name, opt_type)
        cell = get_optimized_cell_parameters(cof_name, opt_type)
        fermi, gap, energy_list = get_gap(cof_name, opt_type)

        # Creates the COF as a pymatgen structure and get density and composition
        struct = Structure(cell, opt_labels, opt_pos, coords_are_cartesian=True)
        density = float(struct.density)
        composition = str(struct.formula)

        symm_info = None

        if symetrize is True:
            symm = SpacegroupAnalyzer(struct, symprec=0.3, angle_tolerance=5.0)
            symm_structure = symm.get_primitive_standard_structure()
            dict_structure = symm_structure.as_dict()

            cell = dict_structure['lattice']['matrix']

            opt_labels = [i['label'] for i in dict_structure['sites']]
            opt_pos = [i['xyz'] for i in dict_structure['sites']]
            composition = str(symm_structure.formula)
            density = float(symm_structure.density)

            # Get the simmetry information of the generated structure
            lattice_type = symm.get_lattice_type()
            space_group = symm.get_space_group_symbol()
            space_group_n = symm.get_space_group_number()

            symm_op = len(symm.get_point_group_operations())
            hall = symm.get_hall()

            symm_info = [lattice_type, space_group, space_group_n, symm_op, hall]
        
        return cell, opt_labels, opt_pos, fermi, gap, energy_list, density, composition, symm_info
    else:
        return False

def read_charges(cof_name):
    
    if charges_concluded(cof_name) == True:
        charges = get_charges(cof_name)
        bond_orders = get_bond_order(cof_name)
        
        return charges, bond_orders
    else:
        return False
    

########################### RASPA related ##########################

def read_raspa_single_ads(path, gas):
    raspa_out_path = os.path.join(path, gas, 'Output', 'System_0')
    if os.path.exists(raspa_out_path):
        out_files = os.listdir(raspa_out_path)
        ads_dict = {'gas': gas,
                    'temperature': None,
                    'pressure': [],
                    'mol/unit_cell':[],
                    'sd_mol/unit_cell':[],
                    'mol/kg':[],
                    'sd_mol/kg':[],
                    'mg/g':[],
                    'sd_mg/g':[],
                    'cm3/g':[],
                    'sd_cm3/g':[],
                    'cm3/cm3':[],
                    'sd_cm3/cm3':[],
                    'enthalpy':[],
                    'sd_enthalpy':[],
                    'heat_capacity_J/mol/K':[],
                    'sd_heat_capacity_J/mol/K':[],
                    'heat_capacity_cal/mol/K':[],
                    'sd_heat_capacity_cal/mol/K':[]
                    }

        for file in out_files:
            tmp = open(os.path.join(raspa_out_path, file), 'r').readlines()
            if 'Simulation finished' in tmp[-3]:
                for line in tmp:
                    if 'External temperature' in line: 
                        ads_dict['temperature'] = float(line.split()[-2])
                    if 'External Pressure' in line: 
                        ads_dict['pressure'] += [float(line.split()[-2])]
                    if '[J/mol/K] +/-  ' in line:
                        ads_dict['heat_capacity_J/mol/K'] += [float(line.split()[1])]
                        ads_dict['sd_heat_capacity_J/mol/K'] += [float(line.split()[-2])]
                    if '[cal/mol/K] +/- ' in line:
                        ads_dict['heat_capacity_cal/mol/K'] += [float(line.split()[1])]
                        ads_dict['sd_heat_capacity_cal/mol/K'] += [float(line.split()[-2])]
                    if '[KJ/MOL]' in line:
                        ads_dict['enthalpy'] += [float(line.split()[0])]
                        ads_dict['sd_enthalpy'] += [float(line.split()[-2])]
                    if 'Average loading absolute [molecules/unit cell]' in line:
                        ads_dict['mol/unit_cell'] += [float(line.split()[-4])]
                        ads_dict['sd_mol/unit_cell'] += [float(line.split()[-2])]
                    if 'Average loading absolute [mol/kg framework]' in line:
                        ads_dict['mol/kg'] += [float(line.split()[-4])]
                        ads_dict['sd_mol/kg'] += [float(line.split()[-2])]
                    if 'Average loading absolute [milligram/gram framework]' in line:
                        ads_dict['mg/g'] += [float(line.split()[-4])]
                        ads_dict['sd_mg/g'] += [float(line.split()[-2])]
                    if 'Average loading absolute [cm^3 (STP)/gr framework]' in line:
                        ads_dict['cm3/g'] += [float(line.split()[-4])]
                        ads_dict['sd_cm3/g'] += [float(line.split()[-2])]
                    if 'Average loading absolute [cm^3 (STP)/cm^3 framework]' in line:
                        ads_dict['cm3/cm3'] += [float(line.split()[-4])]
                        ads_dict['sd_cm3/cm3'] += [float(line.split()[-2])]

        return ads_dict
    else:
        print(raspa_out_path, 'do not exist! Skipping this structure.')
                

########################### JSON creation ##########################
    
def create_COF_json(name):
    
    if os.path.exists(name + '.json') is not True:

        system_info = 'Informations about the system such as name, if it is optimized and other relevant information.'
        geometry_info = 'Informations about the geometry: cell parameters, cell matrix, atomic positions, partial charges, bond orders, simmetry information'
        optimization_info = 'Information about the optimization process such as level of calculations, optimization schema and optimization steps.'
        adsorption_info = 'Information about the adsorption simulation experiments on RASPA2'
        textural_info = 'Information about the textural calculations of the structure such as specific area, pore volume, void fraction.'
        spectrum_info = 'Information about spectra simulation like DRX, FTIR, ssNMR, UV-VIS, Band dispersion, Phonon dispersion...'
        experimental_info = 'Experimental data DRX, FTIR, ssNMR, UV-VIS...'
        
        COF_json = {'system':{'description':system_info,
                               'name':name},
                    'geometry':{'description':geometry_info},
                    'optimization':{'description':optimization_info},
                    'adsorption':{'description':adsorption_info},
                    'textural':{'description':textural_info},
                    'spectrum':{'description':spectrum_info},
                    'experimental':{'description':experimental_info}
                    }
        
        save_json(COF_PATH, name, COF_json)

def add_opt_geo(name, save_opt_history=False):
    
    opt1, bomd, opt2 = get_all_structures_CP2k(name, save_opt_history)
    COF_json = read_json(COF_PATH, name)

    opt_type = False
    
    if len(opt1.keys()) >= 1:
        opt_type = ['OPT_1']
        COF_json['optimization']['opt1'] = opt1
    
    if len(bomd.keys()) >= 1:
        opt_type += ['BOMD']
        COF_json['optimization']['bomd'] = bomd
        
    if len(opt2.keys()) >= 1:
        opt_type += ['OPT_2']
        COF_json['optimization']['opt2'] = opt2

    if opt_type is not False:
    
        cp2k_results = get_cp2k_opt_results(name, opt_type=opt_type[-1])

        if cp2k_results is not False:
            cell, opt_labels, opt_pos, fermi, gap, energy_list, density, composition, symm_info = cp2k_results

            cell_par = Tools.cell_to_cellpar(np.array(cell)).tolist()
            cell_par =  [round(i, 10) for i in cell_par]

            opt_file = os.path.join(COF_PATH, name, opt_type[-1], name + '.cell_opt.out')

            COF_json['system']['geo_opt'] = cp2k_concluded(opt_file)
            COF_json['system']['opt_type'] = opt_type
            COF_json['system']['charges'] = charges_concluded(name)
            
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
            
            save_json(COF_PATH, name, COF_json)

def add_charges(cof_name):

    if charges_concluded(cof_name):
        COF_json = read_json(COF_PATH, cof_name)

        charges, bond_orders = read_charges(cof_name)
        
        COF_json['geometry']['DDEC_charges'] = charges
        COF_json['geometry']['bond_orders'] = bond_orders
        
        save_json(COF_PATH, cof_name, COF_json)
    
def add_RASPA_data(cof_name, gas='CO2'):
    ADS_dir = os.path.join(COF_PATH, cof_name,'ADS')
    
    if os.path.exists(ADS_dir):
        ads_dir = read_raspa_single_ads(ADS_dir, gas)

        COF_json = read_json(COF_PATH, cof_name)
        if gas in COF_json['adsorption'].keys():
            COF_json['adsorption'][gas].update({round(ads_dir['temperature'], 2):ads_dir})
        else:
            COF_json['adsorption'][gas] = {round(ads_dir['temperature'], 2):ads_dir}

        save_json(COF_PATH, cof_name, COF_json)

        
for cof in tqdm(COFs):
    try:
        create_COF_json(cof)
        add_opt_geo(cof, save_opt_history=False)
        add_charges(cof)
        for g in ['CO2', 'methane', 'CO', 'H2']:
            add_RASPA_data(cof, g)
    except Exception:
        print('Erro:', cof)
        traceback.print_exc()
