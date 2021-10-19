# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020
@author: lipel
"""
import os
import tools as Tools
import numpy as np

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    

########################### CP2K related ##########################

def cp2k_concluded(out_file):

    if os.path.exists(out_file):

        for line in open(out_file, 'r').readlines():
            if 'PROGRAM ENDED AT ' in line:
                return True
    else:
        return False

def DDEC_charges_concluded(path, cof_name):
    
    out_path = os.path.join(path, cof_name, 'CHARGES', 'valence_cube_DDEC_analysis.output')
    
    if os.path.exists(out_path):

        for line in open(out_path, 'r').readlines():
            if 'Finished chargemol' in line:
                return True
    else:
        return False

def get_optimized_geometry(path, cof_name, opt_type='OPT_1'):
    
    xyz_path = os.path.join(path, cof_name, opt_type, cof_name + '-pos-1.xyz')
    xyz_opt = open(xyz_path, 'r').readlines()
    
    n_atoms = int(xyz_opt[0].rstrip('\n'))
    opt = [i.split() for i in xyz_opt[-n_atoms:]]
    
    opt_labels = []
    opt_pos = []
    
    for i in opt:
        opt_labels += [i[0]]
        opt_pos += [[float(j) for j in i[1:]]]
        
    return opt_labels, opt_pos

def get_optimized_cell_parameters(path, cof_name, opt_type='OPT_1'):
    
    cell_path = os.path.join(path, cof_name, opt_type, cof_name + '-1.cell')
    cell_opt = open(cell_path, 'r').readlines()
    cell = [float(i) for i in cell_opt[-1].split()[2:-1]]
    cell = [cell[:3], cell[3:6], cell[6:]]
    return cell 

def get_CP2k_results(path, name, opt_type, save_opt_history=False):

    if opt_type.upper() == 'OPT_1':
        pos_file = os.path.join(path, name, 'OPT_1', name + '-pos-1.xyz')
        cell_file = os.path.join(path, name, 'OPT_1', name + '-1.cell')

    if opt_type.upper() == 'BOMD': 
        pos_file = os.path.join(path, name, 'BOMD', name + '-pos-1.xyz') 
        cell_file = os.path.join(path, name, 'BOMD', name + '-1.cell')

    if opt_type.upper() == 'OPT_2':
        pos_file = os.path.join(path, name, 'OPT_2', name + '-pos-1.xyz')
        cell_file = os.path.join(path, name, 'OPT_2', name + '-1.cell')

    results = {}
    pos_file_exist = os.path.exists(pos_file)
    cell_file_exist = os.path.exists(cell_file)
    if pos_file_exist is True and cell_file_exist:
        if save_opt_history:

            tmp = open(pos_file).readlines()
            n_atoms = int(tmp[0])
            n_steps = [i for i in tmp if 'i' in i]

            for i in range(len(n_steps)):
                start = i*(n_atoms + 2)+2
                end =i*(n_atoms + 2) + 2 + n_atoms
                step = tmp[start -1].split()[2].rstrip(',')
                results[step] = {
                                'cell_matrix':[],
                                'atom_label':[],
                                'atom_pos':[],
                                'energy':[]}

                tmp_pos = tmp[start: end]
                energy = float(tmp[start -1].split()[-1])
                results[step]['energy'] = energy

                tmp_atom_labels, tmp_atom_pos = [], []
                for j in range(len(tmp_pos)):
                    tmp_atom_labels += [tmp_pos[j].split()[0]]
                    tmp_atom_pos += [[float(k) for k in tmp_pos[j].split()[1:]]]

                results[step]['atom_label'] = tmp_atom_labels
                results[step]['atom_pos'] = tmp_atom_pos

            tmp_cell = open(cell_file).readlines()[1:]
            for t in range(len(tmp_cell)):
                step = tmp_cell[t].split()[0]
                cell = np.reshape(np.array([float(i) for i in tmp_cell[t].split()[2:-1]]), (3,3)).tolist()
                results[step]['cell_matrix'] = cell
        
        if save_opt_history is False:
            tmp = open(pos_file).readlines()
            results['energy'] = [float(i.split()[-1]) for i in tmp if 'i' in i]
        
        return results
    else:
        return False

def get_all_structures_CP2k(path, name, save_opt_history=False):
    
    initial_file = os.path.join(path, name, 'OPT_1', name + '.xyz')

    pos_1_file = os.path.join(path, name, 'OPT_1', name + '-pos-1.xyz')
    pos_2_file = os.path.join(path, name, 'OPT_2', name + '-pos-1.xyz')
    pos_bomd_file = os.path.join(path, name, 'BOMD', name + '-pos-1.xyz')

    cell_1_file = os.path.join(path, name, 'OPT_1', name + '-1.cell')
    cell_bomd_file = os.path.join(path, name, 'BOMD', name + '-1.cell')
    cell_2_file = os.path.join(path, name, 'OPT_2', name + '-1.cell')

    initial = {}
    opt_1 = {}
    bomd = {}
    opt_2 = {}

    if cp2k_concluded(os.path.join(path, name, 'OPT_1', name + '.cell_opt.out')) is False:
        try:
            tmp = open(initial_file).readlines()
            n_atoms = int(tmp[0])

            initial = { 'cell_parameters':[],
                        'cell_matrix':[],
                        'atom_labels':[],
                        'atom_pos':[],
                        'energy':[]}

            tmp_pos = tmp[2:]

            tmp_atom_labels, tmp_atom_pos = [], []
            for j in range(len(tmp_pos)):
                tmp_atom_labels += [tmp_pos[j].split()[0]]
                tmp_atom_pos += [[float(k) for k in tmp_pos[j].split()[1:]]]

            initial['atom_labels'] = tmp_atom_labels
            initial['atom_pos'] = tmp_atom_pos

            tmp_cell = [float(i) for i in tmp[1].split()]
            cell_matrix = Tools.cellpar_to_cell(tmp_cell)
            initial['cell_parameters'] = tmp_cell
            initial['cell_matrix'] = cell_matrix.tolist()
        except:
            None

    if cp2k_concluded(os.path.join(path, name, 'OPT_1', name + '.cell_opt.out')):

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

    if cp2k_concluded(os.path.join(path, name, 'BOMD', name + '.bomd.out')):

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
    
    if cp2k_concluded(os.path.join(path, name, 'OPT_2', name + '.cell_opt.out')):

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

    return initial, opt_1, bomd, opt_2
    
def get_DDEC_charges(path, cof_name):
    
    charges_path = os.path.join(path, cof_name, 'CHARGES', 'DDEC6_even_tempered_net_atomic_charges.xyz')
    
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

def get_lowdin_charges(path, cof_name):
    
    charges_path = os.path.join(path, cof_name, 'CHARGES', f'{cof_name}.scf.out')
    
    if os.path.exists(charges_path):
        tmp = open(charges_path, 'r').readlines()
        for i in tmp:
            if '- Atoms: ' in i:
                n_atoms = int(i.split()[-1])
        for i in range(len(tmp)):
            if ' LOWDIN POPULATION ANALYSIS' in tmp[i]:
                atom_list = tmp[i + 3: i+3+n_atoms]

        charges = []
        
        for i in atom_list:
            charges += [float(i.split()[-1])]
            
        return charges
    else:
        return False

def get_hirshfeld_charges(path, cof_name):
    
    charges_path = os.path.join(path, cof_name, 'CHARGES', f'{cof_name}.scf.out')
    
    if os.path.exists(charges_path):
        tmp = open(charges_path, 'r').readlines()
        for i in tmp:
            if '- Atoms: ' in i:
                n_atoms = int(i.split()[-1])
        for i in range(len(tmp)):
            if 'Hirshfeld Charges' in tmp[i]:
                atom_list = tmp[i + 3: i+3+n_atoms]

        charges = []
        
        for i in atom_list:
            charges += [float(i.split()[-1])]
            
        return charges
    else:
        return False

def get_mulliken_charges(path, cof_name):
    
    charges_path = os.path.join(path, cof_name, 'CHARGES', f'{cof_name}.scf.out')
    
    if os.path.exists(charges_path):
        tmp = open(charges_path, 'r').readlines()
        for i in tmp:
            if '- Atoms: ' in i:
                n_atoms = int(i.split()[-1])
        for i in range(len(tmp)):
            if ' Mulliken Population Analysis' in tmp[i]:
                atom_list = tmp[i + 3: i+3+n_atoms]

        charges = []
        
        for i in atom_list:
            charges += [float(i.split()[-1])]
            
        return charges
    else:
        return False

def get_bond_order(path, cof_name):
    
    charges_path = os.path.join(path, cof_name, 'CHARGES', 'DDEC6_even_tempered_bond_orders.xyz')
    
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
        
def get_gap(path, cof_name, opt_type='OPT_1'):
    
    if opt_type == 'OPT_1':
        file_sufix = '.cell_opt.out'
    if opt_type == 'BOMD':
        file_sufix = '.bomd.out'
    if opt_type == 'OPT_2':
        file_sufix = '.cell_opt.out'
    if opt_type == 'CHARGES':
        file_sufix = '.scf.out'

    out_path = os.path.join(path, cof_name, opt_type, cof_name + file_sufix)
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
            if 'ENERGY| Total FORCE_EVAL ( QS ) energy' in i:
                t_energy += [float(i.split()[-1])]

    return fermi, gap, t_energy
        
def get_cp2k_opt_results(path, cof_name, opt_type='OPT_1', symetrize=False):

    if opt_type == 'OPT_1':
        file_sufix = '.cell_opt.out'
    if opt_type == 'BOMD':
        file_sufix = '.bomd.out'
    if opt_type == 'OPT_2':
        file_sufix = '.cell_opt.out'
    if opt_type == 'CHARGES':
        file_sufix = '.scf.out'

    out_path = os.path.join(path, cof_name, opt_type, cof_name + file_sufix)
    
    if cp2k_concluded(out_path) == True:
        opt_labels, opt_pos = get_optimized_geometry(path, cof_name, opt_type)
        cell = get_optimized_cell_parameters(path, cof_name, opt_type)
        fermi, gap, energy_list = get_gap(path, cof_name, opt_type)

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

def read_DDEC_charges(path, cof_name):
    
    if DDEC_charges_concluded(path, cof_name) == True:
        charges = get_DDEC_charges(path, cof_name)
        bond_orders = get_bond_order(path, cof_name)
        
        return charges, bond_orders
    else:
        return False

def read_VIB_data(path, cof_name):
    
    tmp = open(os.path.join(path, cof_name + '-VIBRATIONS-1.mol'), 'r').readlines()

    freq_begin_pos = 0
    freq_end_pos = 0
    int_begin_pos = 0
    
    for i in range(len(tmp)):
        if ' [FREQ]' in tmp[i]:
            freq_begin_pos = i+ 1
        if ' [FR-COORD]' in tmp[i]:
            freq_end_pos = i
        if ' [INT]' in tmp[i]:
            int_begin_pos = i + 1
            
    freq = [float(i.rstrip()) for i in tmp[freq_begin_pos:freq_end_pos]]
    inten = [float(i.rstrip()) for i in tmp[int_begin_pos:]]
    
    return freq, inten

def read_charges(path, cof_name):

    if DDEC_charges_concluded(path, cof_name) == True:
        lowdin_charges = get_lowdin_charges(path, cof_name)
        hirshfeld_charges = get_hirshfeld_charges(path, cof_name)
        mulliken_charges = get_mulliken_charges(path, cof_name)
    
        return lowdin_charges, hirshfeld_charges, mulliken_charges
    else:
        return None
