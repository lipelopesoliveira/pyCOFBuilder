# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020

@author: lipel
"""

import os
import numpy as np

from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from scipy.spatial.transform import Rotation as R

import pycofbuilder.tools as Tools
from pycofbuilder.building_block import Building_Block

class Reticulum():

    def __init__(self, bb_lib='bb_lib', verbosity=False):

        _ROOT = os.path.abspath(os.path.dirname(__file__))

        self.verbosity = verbosity
        self.available_2D_topologies = ['HCB', 'HCB_A', 'SQL', 'SQL_A', 'KGM', 'KGM_A', 'HXL', 'HXL_A', 'KGD', 'KGD_A']
        self.available_3D_topologies = ['dia', 'bor', 'srs', 'pts', 'ctn', 'rra', 'fcc', 'lon', 'stp', 'acs', 'tbo', 'bcu', 'fjh', 'ceq']

        self.available_stacking = {'HCB': ['AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'HCB_A': ['AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'SQL':['AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'SQL_A':['AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'KGD': ['A', 'AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2']
                                   }
        self.available_topologies = self.available_2D_topologies + self.available_3D_topologies
        self.lib_bb = bb_lib
        self.main_path = os.path.join(_ROOT, 'data')
        self.lib_path = os.path.join(self.main_path, bb_lib)
        self.out_path = os.path.join(os.getcwd(), 'out')
        self.name = None
        self.topology = None
        self.dimention = None
        self.lattice = None
        self.lattice_sgs = None
        self.space_group = None
        self.space_group_n = None
        self.stacking = None
        self.mass = None
        self.composition = None
        self.charge = 0
        self.multiplicity = 1
        self.chirality = False
        self.atom_labels = []
        self.atom_pos = []
        self.lattice = [[], [], []]
        self.symm_tol = 0.2
        self.angle_tol = 0.2
        self.n_atoms = self.n_atoms()

    def n_atoms(self):
        ''' Returns the number of atoms in the unitary cell'''
        return len(self.atom_labels)

    def print_available_topologies(self, dimensionality='all'):

        if dimensionality == 'all' or dimensionality == '2D':
            print('Available 2D Topologies:')
            for i in self.available_2D_topologies:
                print(i.upper())

        if dimensionality == 'all' or dimensionality == '3D':
            print('Available 3D Topologies:')
            for i in self.available_3D_topologies:
                print(i.upper())

    def create_hcb_structure(self, name_bb_a, name_bb_b, stacking='AA', bond_atom='N', c_parameter_base=3.6, print_result=True, slab=10, shift_vector=[1.0,1.0,0], tilt_angle=5):
        '''Creates a COF with HCB network'''

        self.topology = 'hcb'
        self.dimension = 2

        bb_1 = Building_Block(name_bb_a, self.lib_bb)
        bb_2 = Building_Block(name_bb_b, self.lib_bb)

        self.name = f'{bb_1.name}-{bb_2.name}-HCB-{stacking}'
        self.charge = bb_1.charge + bb_2.charge
        self.chirality = bb_1.chirality or bb_2.chirality

        if self.verbosity is True:
            print('Starting the creation of a hcb net')

        if bb_1.connectivity != 3:
            print('Building block A must present connectivity 3 insted of', len(bb_1.connectivity))
            return None
        if bb_2.connectivity != 3:
            print('Building block B must present connectivity 3 insted of', len(bb_2.connectivity))
            return None

        # Calculate the cell parameter based on the size of the building blocks
        size_a = bb_1.size
        if self.verbosity:
            print('BB_A size:', size_a)
        size_b = bb_2.size
        if self.verbosity:
            print('BB_B size:', size_b)

        # Defines the cell parameter a
        a = np.cos(np.radians(30))*2*(size_a[0] + size_b[0])

        if self.verbosity:
            print('Calculated cell parameter a:', a)

        # Gets the maximum distace in the z axis to create the c parameter
        delta_a = abs(max(np.transpose(bb_1.atom_pos)[2])) + abs(min(np.transpose(bb_1.atom_pos)[2]))
        delta_b = abs(max(np.transpose(bb_2.atom_pos)[2])) + abs(min(np.transpose(bb_2.atom_pos)[2]))

        # Build the matrix of unitary hexagonal cell
        if stacking == 'A':
            c = slab
        else:
            c = c_parameter_base + max([delta_a, delta_b])

        # Define the cell lattice
        lattice = [[a, 0, 0], 
                   [-0.5*a, np.sqrt(3)/2*a, 0], 
                   [0, 0, c]]

        if self.verbosity is True:
            print('Unitary cell built:', lattice)

        # Adds the building block A in the origin of the unitary cell (A1 site)
        final_label = bb_1.atom_labels
        final_pos = bb_1.atom_pos

        # Rotates and add the building block A to site A2 of unitary cell
        r_pos_a_2 = np.dot(bb_2.atom_pos, R.from_euler('z', 180, degrees=True).as_matrix()) + np.array([0, np.sqrt(3)/3, 0])*a

        final_pos = np.vstack((final_pos, r_pos_a_2))
        final_label += bb_2.atom_labels

        # Changes the X atoms by the desirede bond_atom
        final_label, final_pos = Tools.change_X_atoms(final_label, final_pos, bond_atom)

        # Creates the COF as a pymatgen structure
        struct = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)

        # Remove duplicate atoms and translate the structure to the center of the cell
        struct.sort(reverse=True)
        struct.merge_sites(tol=.5, mode='delete')
        struct.translate_sites(range(len(struct.as_dict()['sites'])), [0, 0, 0.5], frac_coords=True, to_unit_cell=True)

        # Simetrizes the structure using pymatgen
        symm = SpacegroupAnalyzer(struct, symprec=self.symm_tol, angle_tolerance=self.angle_tol)
        struct_symm_prim = symm.get_primitive_standard_structure()

        if stacking not in self.available_stacking['HCB']:
            raise Exception(f"""{stacking} is not in the available stack list for HCB net.
    Available options are: {self.available_stacking['HCB']}""")

        # Create A stacking, a 2D isolated sheet with slab
        if stacking == 'A':
            self.stacking = 'A'
            self.symm_structure = struct_symm_prim

        # Create AA staking. By default one sheet per unitary cell is used
        if stacking == 'AA':
            self.stacking = 'AA'
            self.symm_structure = struct_symm_prim

        # Create AB1 staking. 
        if stacking == 'AB1':
            self.stacking = 'AB1'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal
            B = ion_conv_crystal

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_1 = Structure(lattice, AB_label + AB_label, AB, coords_are_cartesian=True)

            B_list = np.arange(len(AB_label)) + len(AB_label)

            AB_1.translate_sites(B_list, [2/3, 1/3, 0.5], frac_coords=True, to_unit_cell=True)
            AB_1_symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = AB_1_symm.get_refined_structure()

        # Create AB2 stacking
        if stacking == 'AB2':
            self.stacking = 'AB2'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal
            B = ion_conv_crystal

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_2 = Structure(lattice, AB_label + AB_label, AB, coords_are_cartesian=True)

            B_list = np.arange(len(AB_label)) + len(AB_label)
            AB_2.translate_sites(B_list, [1/2, 0, 0.5], frac_coords=True, to_unit_cell=True)
            AB_2_symm = SpacegroupAnalyzer(AB_2, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = AB_2_symm.get_refined_structure()
            
        # Create AAl stacking. Hexagonal cell with two sheets per cell shifited by the shift_vector in angstroms. 
        if stacking == 'AAl':
            self.stacking = 'AAl'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 2)

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal * np.array([1, 1, 0.5])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(*Tools.cell_to_cellpar(cell))
            shift_vector = np.dot(r, np.array(shift_vector))

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal * np.array([1, 1, 1.5]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            self.symm_structure = Structure(lattice, AB_label+AB_label, AB, coords_are_cartesian=False)

        # Create AA tilted stacking. Tilted Hexagonal cell by tilt_angle with two sheets per cell shifited by the shift_vector in angstroms.
        if stacking == 'AAt':
            self.stacking = 'AAt'

            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz'] for i in struct_symm_prim.as_dict()['sites']])
            cell = struct_symm_prim.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta= cell['beta'] - tilt_angle
            gamma = cell['gamma']

            new_cell = Tools.cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + np.array([0, 0, 0.5*c_cell]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(new_cell)
            self.symm_structure = Structure(lattice, AB_label+AB_label, AB, coords_are_cartesian=True)

        if stacking == 'ABC1':
            self.stacking = 'ABC1'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(ion_conv_crystal*(1, 1, 1) + (2/3, 1/3, 0))
            C = Tools.translate_inside(ion_conv_crystal*(1, 1, 1/3) + (4/3, 2/3, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label+ABC_label, ABC, coords_are_cartesian=False)
            ABC_f_symm = SpacegroupAnalyzer(ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_primitive_standard_structure()

        if stacking == 'ABC2':
            self.stacking = 'ABC2'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim['sites']])
            cell = np.array(struct_symm_prim['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(ion_conv_crystal*(1, 1, 1) + (1/3, 0, 0))
            C = Tools.translate_inside(ion_conv_crystal*(1, 1, 1/3) + (2/3, 0, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label+ABC_label, ABC, coords_are_cartesian=False)
            ABC_f_symm = SpacegroupAnalyzer(ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_primitive_standard_structure()

        dict_structure = self.symm_structure.as_dict()

        self.lattice = dict_structure['lattice']['matrix']

        self.atom_labels = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.n_atoms = len(self.symm_structure)
        self.composition = self.symm_structure.formula

        if self.verbosity is True:
            print(self.symm_structure)

        # Get the simmetry information of the generated structure
        lattice_type = symm.get_lattice_type()
        self.lattice_type = lattice_type
        self.space_group = symm.get_space_group_symbol()
        self.space_group_n = symm.get_space_group_number()

        symm_op = symm.get_point_group_operations()
        self.hall = symm.get_hall()

        if print_result == True:
            Tools.print_result(self.name, str(self.lattice_type), str(self.hall[0:2]), str(self.space_group), str(self.space_group_n), len(symm_op))

        return [self.name, str(self.lattice_type), str(self.hall[0:2]), str(self.space_group), str(self.space_group_n), len(symm_op)]

    def create_hcb_a_structure(self, name_bb_a, name_bb_b, stacking='AA', bond_atom='N', c_parameter_base=3.6, print_result=True, slab=10, shift_vector=[1,1,0], tilt_angle=5):

        self.topology = 'hcb-a'
        self.dimension = 2

        bb_triangular = Building_Block(name_bb_a, self.lib_bb, verbosity=self.verbosity)

        bb_linear = Building_Block(name_bb_b, self.lib_bb, verbosity=self.verbosity)

        self.charge = bb_linear.charge + bb_triangular.charge
        self.chirality = bb_linear.chirality or bb_triangular.chirality

        self.name = f'{bb_triangular.name}-{bb_linear.name}-HCB_A-{stacking}'

        Tools.print_comand(f'Starting the creation of {self.name}', self.verbosity, ['debug', 'high'])

        if bb_triangular.connectivity != 3:
            print('Building block A must present connectivity 3 insted of', bb_triangular.connectivity)
            return None
        if bb_linear.connectivity != 2:
            print('Building block B must present connectivity 2 insted of', bb_linear.connectivity)
            return None

        # Calcula o parâmetro de célula com base no tamanho dos blocos de construção
        size_a = bb_triangular.size
        size_b = bb_linear.size

        Tools.print_comand(f'BB_A size: {size_a}', self.verbosity, ['debug'])
        Tools.print_comand(f'BB_B size: {size_b}', self.verbosity, ['debug'])

        a = 2*np.cos(np.radians(30))*2*(size_a[0] + size_b[0])

        Tools.print_comand(f'Calculated cell parameter a: {a}', self.verbosity, ['debug'])

        # Mede o valor do delta c dos blocos de construção
        delta_a = abs(max(np.transpose(bb_triangular.atom_pos)[2])) + abs(min(np.transpose(bb_triangular.atom_pos)[2]))
        delta_b = abs(max(np.transpose(bb_linear.atom_pos)[2])) + abs(min(np.transpose(bb_linear.atom_pos)[2]))

        # Build the matrix of unitary hexagonal cell
        if stacking == 'A':
            c = slab
        else:
            c = c_parameter_base + max([delta_a, delta_b])

        # Define the cell lattice
        lattice = [[a, 0, 0], 
                   [-0.5*a, np.sqrt(3)/2*a, 0], 
                   [0, 0, c]]

        if self.verbosity is True:
            print('Unitary cell built:', lattice)

        # Adiciona o bloco A na origem da célula unitária (sítio A1)
        final_label = bb_triangular.atom_labels
        final_pos = bb_triangular.atom_pos

        # Rotaciona o bloco A e adiciona no sítio A2 da célula unitária
        r_pos_a_1 = np.dot(bb_triangular.atom_pos, R.from_euler('z', 180, degrees=True).as_matrix()) + np.array([0, np.sqrt(3)/3, 0])*a

        final_pos = np.vstack((final_pos, r_pos_a_1))
        final_label += bb_triangular.atom_labels

        # Adiciona o bloco B no sítio B1 da célula unitária
        final_pos = np.vstack((final_pos, bb_linear.atom_pos + np.array([0, np.sqrt(3)/6, 0])*a))
        final_label += bb_linear.atom_labels

        # Rotaciona o bloco B e o adiciona no sítio B2 da célula unitária
        r_pos_b_1 = np.dot(bb_linear.atom_pos, R.from_euler('z', 120, degrees=True).as_matrix()) + np.array([-1/4, 5*np.sqrt(3)/12, 0])*a

        final_pos = np.vstack((final_pos, r_pos_b_1))
        final_label += bb_linear.atom_labels

        # Rotaciona o bloco B e o adiciona no sítio B3 da célula unitária
        r_pos_b_2 = np.dot(bb_linear.atom_pos, R.from_euler('z', 240, degrees=True).as_matrix()) + np.array([1/4, 5*np.sqrt(3)/12, 0])*a

        final_pos = np.vstack((final_pos, r_pos_b_2))
        final_label += bb_linear.atom_labels

        # Changes the X atoms by the desirede bond_atom
        final_label, final_pos = Tools.change_X_atoms(final_label, final_pos, bond_atom)

        # Cria a estrutura como entidade do pymatgen
        struct = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)

        # Remove os átomos duplicados
        struct.sort(reverse=True)
        struct.merge_sites(tol=.5, mode='delete')
        struct.translate_sites(range(len(struct.as_dict()['sites'])), [0, 0, 0.5], frac_coords=True, to_unit_cell=True)

        # Simetriza a estrutura
        try:
            symm = SpacegroupAnalyzer(struct, symprec=0.3, angle_tolerance=3.0)
            struct_symm_prim = symm.get_refined_structure()
        except Exception:
            return None

        if stacking not in self.available_stacking['HCB_A']:
            raise Exception(f"""{stacking} is not in the available stack list for HCB-A net.
    Available options are: {self.available_stacking['HCB_A']}""")

        # Create A stacking. The slab is defined by the c_cell parameter
        if stacking == 'A':
            self.stacking = 'A'
            self.symm_structure = struct_symm_prim

        # Create AA staking. By default one sheet per unitary cell is used
        if stacking == 'AA':
            self.stacking = 'AA'
            self.symm_structure = struct_symm_prim

        if stacking == 'AB1':
            self.stacking = 'AB1'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal
            B = ion_conv_crystal

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_1 = Structure(lattice, AB_label + AB_label, AB, coords_are_cartesian=True)

            B_list = np.arange(len(AB_label)) + len(AB_label)

            AB_1.translate_sites(B_list, [2/3, 1/3, 0.5], frac_coords=True, to_unit_cell=True)
            AB_1_symm = SpacegroupAnalyzer(AB_1, symprec=0.05, angle_tolerance=.5)

            self.symm_structure = AB_1_symm.get_refined_structure()

        if stacking == 'AB2':
            self.stacking = 'AB2'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal
            B = ion_conv_crystal

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_2 = Structure(lattice, AB_label + AB_label, AB, coords_are_cartesian=True)

            B_list = np.arange(len(AB_label)) + len(AB_label)
            AB_2.translate_sites(B_list, [1/2, 0, 0.5], frac_coords=True, to_unit_cell=True)
            AB_2_symm = SpacegroupAnalyzer(AB_2, symprec=0.05, angle_tolerance=0.5)

            self.symm_structure = AB_2_symm.get_refined_structure()

        # Create AAl stacking. Hexagonal cell with two sheets per cell shifited by the shift_vector in angstroms. 
        if stacking == 'AAl':
            self.stacking = 'AAl'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 2)

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal * np.array([1, 1, 0.5])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(*Tools.cell_to_cellpar(cell))
            shift_vector = np.dot(r, np.array(shift_vector))

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal * np.array([1, 1, 1.5]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            self.symm_structure = Structure(lattice, AB_label+AB_label, AB, coords_are_cartesian=False)

        # Create AA tilted stacking. Tilted Hexagonal cell with two sheets per cell tilted by tilt_angle. 
        if stacking == 'AAt':
            self.stacking = 'AAt'

            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz'] for i in struct_symm_prim.as_dict()['sites']])
            cell = struct_symm_prim.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta= cell['beta'] - tilt_angle
            gamma = cell['gamma']

            new_cell = Tools.cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + np.array([0, 0, 0.5*c_cell]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(new_cell)
            self.symm_structure = Structure(lattice, AB_label+AB_label, AB, coords_are_cartesian=True)

        if stacking == 'ABC1':
            self.stacking = 'ABC1'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(ion_conv_crystal*(1, 1, 1) + (2/3, 1/3, 0))
            C = Tools.translate_inside(ion_conv_crystal*(1, 1, 1/3) + (4/3, 2/3, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label+ABC_label, ABC, coords_are_cartesian=False)
            ABC_f_symm = SpacegroupAnalyzer(ABC_f, symprec=0.2, angle_tolerance=2.0)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        if stacking == 'ABC2':
            self.stacking = 'ABC2'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(ion_conv_crystal*(1, 1, 1) + (1/3, 0, 0))
            C = Tools.translate_inside(ion_conv_crystal*(1, 1, 1/3) + (2/3, 0, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label+ABC_label, ABC, coords_are_cartesian=False)
            ABC_f_symm = SpacegroupAnalyzer(ABC_f, symprec=0.2, angle_tolerance=2.0)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        dict_structure = self.symm_structure.as_dict()

        self.lattice = dict_structure['lattice']['matrix']

        self.atom_labels = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.n_atoms = len(self.symm_structure)
        self.composition = self.symm_structure.formula

        Tools.print_comand(self.symm_structure, self.verbosity, ['debug'])

        # Get the simmetry information of the generated structure
        self.lattice_type = symm.get_lattice_type()
        self.space_group = symm.get_space_group_symbol()
        self.space_group_n = symm.get_space_group_number()

        symm_op = symm.get_point_group_operations()
        self.hall = symm.get_hall()

        if print_result == True:
            Tools.print_result(self.name, str(self.lattice_type), str(self.hall[0:2]), str(self.space_group), str(self.space_group_n), len(symm_op))

        return [self.name, str(self.lattice_type), str(self.hall[0:2]), str(self.space_group), str(self.space_group_n), len(symm_op)]

    def create_sql_structure(self, name_a, name_b, stacking='AA', bond_atom='N', c_cell=3.6, print_result=True):

        self.topology = 'sql'
        self.dimension = 2

        bb_1 = Building_Block(name_a, self.lib_bb, verbosity=self.verbosity)

        bb_2 = Building_Block(name_b, self.lib_bb, verbosity=self.verbosity)

        self.charge = bb_1.charge + bb_2.charge
        self.chirality = bb_1.chirality or bb_2.chirality

        self.name = f'{bb_1.name}-{bb_2.name}-SQL-{stacking}'

        if self.verbosity is True:
            print('Starting the creation of a sql net')

        if bb_1.connectivity != 4:
            print(f'Building block A ({name_a}) must present connectivity 4 insted of', bb_1.connectivity)
            return None
        if bb_2.connectivity != 4:
            print(f'Building block B ({name_b}) must present connectivity 4 insted of', bb_2.connectivity)
            return None

        # Calculates the cell parameter based on the building blocks size
        size_a = bb_1.size
        if self.verbosity:
            print('BB_A size:', size_a)
        size_b = bb_2.size
        if self.verbosity:
            print('BB_B size:', size_b)

        a = 2*(np.average(size_a) + np.average(size_b))/np.sqrt(2)

        if self.verbosity:
            print('Calculated cell parameter a:', a)

        # Mede o valor do delta c dos blocos de construção
        delta_a = abs(max(np.transpose(bb_1.atom_pos)[2])) + abs(min(np.transpose(bb_1.atom_pos)[2]))
        delta_b = abs(max(np.transpose(bb_2.atom_pos)[2])) + abs(min(np.transpose(bb_2.atom_pos)[2]))

        # Constrói a matriz da célula unitária tetragonal
        lattice = [[a, 0, 0], [0, a, 0], [0, 0, c_cell + max([delta_a, delta_b])]]

        if self.verbosity is True:
            print('Unitary cell built:')
            print('a =', lattice[0])
            print('b =', lattice[1])
            print('c =', lattice[2])

        # Adiciona o bloco 1 na origem da célula unitária (sítio A1)
        final_pos = np.dot(bb_1.atom_pos, R.from_euler('z', 45, degrees=True).as_matrix())
        final_label = bb_1.atom_labels

        # Adiciona o bloco 2 no sítio B1 da célula unitária
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', 45, degrees=True).as_matrix()) + np.array([0.5, 0.5, 0])*a))
        final_label += bb_2.atom_labels

        # Changes the X atoms by the desirede bond_atom
        final_label, final_pos = Tools.change_X_atoms(final_label, final_pos, bond_atom)

        # Cria a estrutura como entidade do pymatgen
        struct = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)

        # Remove os átomos duplicados
        struct.sort(reverse=True)
        struct.merge_sites(tol=.5, mode='delete')
        struct.translate_sites(range(len(struct.as_dict()['sites'])), [0, 0, 0.5], frac_coords=True, to_unit_cell=True)

        # Simetriza a estrutura
        symm = SpacegroupAnalyzer(struct, symprec=1, angle_tolerance=5.0)
        struct_symm_prim = symm.get_refined_structure()

        if stacking not in self.available_stacking['SQL']:
            raise Exception(f"""{stacking} is not in the available stack list for SQL net.
    Available options are: {self.available_stacking['SQL']}""")

        if stacking == 'AA':
            self.stacking = 'AA'
            self.symm_structure = struct_symm_prim

        if stacking == 'AB1':
            self.stacking = 'AB1'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal*(1, 1, 0.5)
            B = Tools.translate_inside(ion_conv_crystal*(1, 1, 1.5) + (1/4, 1/4, 0))

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_1 = Structure(lattice, AB_label + AB_label, AB, coords_are_cartesian=False)
            AB_1_symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = AB_1_symm.get_refined_structure()

        if stacking == 'AB2':
            self.stacking = 'AB2'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal*(1, 1, 0.5)
            B = Tools.translate_inside(ion_conv_crystal*(1, 1, 1.5) + (1/2, 0, 0))

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_2 = Structure(lattice, AB_label + AB_label, AB, coords_are_cartesian=False)
            AB_2_symm = SpacegroupAnalyzer(AB_2, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = AB_2_symm.get_refined_structure()

        # Create AAl stacking. Tetragonal cell with two sheets per cell shifited by the shift_vector in angstroms. 
        if stacking == 'AAl':
            self.stacking = 'AAl'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 2)

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal * np.array([1, 1, 0.5])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(*Tools.cell_to_cellpar(cell))
            shift_vector = np.dot(r, np.array(shift_vector))

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal * np.array([1, 1, 1.5]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            self.symm_structure = Structure(lattice, AB_label+AB_label, AB, coords_are_cartesian=False)

        # Create AA tilted stacking. Tilted tetragonal cell with two sheets per cell tilted by tilt_angle. 
        if stacking == 'AAt':
            self.stacking = 'AAt'

            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz'] for i in struct_symm_prim.as_dict()['sites']])
            cell = struct_symm_prim.as_dict()['lattice']

            # Shift the cell by the tilt angle
            print(cell)
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta= cell['beta'] - tilt_angle
            gamma = cell['gamma']

            new_cell = Tools.cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + np.array([0, 0, 0.5*c_cell]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(new_cell)
            self.symm_structure = Structure(lattice, AB_label+AB_label, AB, coords_are_cartesian=True)

        if stacking == 'ABC1':
            self.stacking = 'ABC1'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(ion_conv_crystal*(1, 1, 1) + (2/3, 1/3, 0))
            C = Tools.translate_inside(ion_conv_crystal*(1, 1, 1/3) + (4/3, 2/3, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label+ABC_label, ABC, coords_are_cartesian=False)
            ABC_f_symm = SpacegroupAnalyzer(ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        if stacking == 'ABC2':
            self.stacking = 'ABC2'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(ion_conv_crystal*(1, 1, 1) + (1/3, 0, 0))
            C = Tools.translate_inside(ion_conv_crystal*(1, 1, 1/3) + (2/3, 0, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label+ABC_label, ABC, coords_are_cartesian=False)
            ABC_f_symm = SpacegroupAnalyzer(ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        dict_structure = self.symm_structure.as_dict()

        self.lattice = dict_structure['lattice']['matrix']

        self.atom_labels = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.n_atoms = len(self.symm_structure)
        self.composition = self.symm_structure.formula

        if self.verbosity is True:
            print(self.symm_structure)

        # Get the simmetry information of the generated structure
        self.lattice_type = symm.get_lattice_type()
        self.space_group = symm.get_space_group_symbol()
        self.space_group_n = symm.get_space_group_number()

        symm_op = symm.get_point_group_operations()
        self.hall = symm.get_hall()

        if print_result == True:
            Tools.print_result(self.name, str(self.lattice_type), str(self.hall[0:2]), str(self.space_group), str(self.space_group_n), len(symm_op))

        return [self.name, str(self.lattice_type), str(self.hall[0:2]), str(self.space_group), str(self.space_group_n), len(symm_op)]

    def create_sql_a_structure(self, name_a, name_b, stacking='AA', bond_atom='N', c_cell=3.6, print_result=True):

        self.topology = 'SQL_A'
        self.dimension = 2

        bb_1 = Building_Block(name_a, self.lib_bb, verbosity=self.verbosity)

        bb_2 = Building_Block(name_b, self.lib_bb, verbosity=self.verbosity)

        self.charge = bb_1.charge + bb_2.charge
        self.chirality = bb_1.chirality or bb_2.chirality

        self.name = f'{bb_1.name}-{bb_2.name}-SQL_A-{stacking}'

        if self.verbosity is True:
            print('Starting the creation of a sql net')

        if bb_1.connectivity != 4:
            print(f'Building block A ({name_a}) must present connectivity 4 insted of', bb_1.connectivity)
            return None
        if bb_2.connectivity != 2:
            print(f'Building block B ({name_b}) must present connectivity 2 insted of', bb_2.connectivity)
            return None

        # Calcula o parâmetro de célula com base no tamanho dos blocos de construção
        size_a = bb_1.size
        if self.verbosity:
            print('BB_A size:', size_a)
        size_b = bb_2.size
        if self.verbosity:
            print('BB_B size:', size_b)

        a = 4*(np.average(size_a) + np.average(size_b))/np.sqrt(2)

        if self.verbosity:
            print('Calculated cell parameter a:', a)

        # Mede o valor do delta c dos blocos de construção
        delta_a = abs(max(np.transpose(bb_1.atom_pos)[2])) + abs(min(np.transpose(bb_1.atom_pos)[2]))
        delta_b = abs(max(np.transpose(bb_2.atom_pos)[2])) + abs(min(np.transpose(bb_2.atom_pos)[2]))

        # Constrói a matriz da célula unitária tetragonal
        lattice = [[a, 0, 0], [0, a, 0], [0, 0, c_cell + max([delta_a, delta_b])]]

        if self.verbosity is True:
            print('Unitary cell built:')
            print('a =', lattice[0])
            print('b =', lattice[1])
            print('c =', lattice[2])

        # Adiciona o bloco 1 na origem da célula unitária (sítio A1)
        final_pos = np.dot(bb_1.atom_pos, R.from_euler('z', 45, degrees=True).as_matrix())
        final_label = bb_1.atom_labels

        # Adiciona o bloco 1 no centro da célula unitária (sítio A2)
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler('z', -45, degrees=True).as_matrix()) + np.array([0.5, 0.5, 0])*a))
        final_label += bb_1.atom_labels

        # Adiciona o bloco 2 na origem da célula unitária (sítio B1)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', 45, degrees=True).as_matrix()) + np.array([0.25, 0.25, 0])*a))
        final_label += bb_2.atom_labels

        # Adiciona o bloco 2 na origem da célula unitária (sítio B2)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', 135, degrees=True).as_matrix()) + np.array([0.75, 0.25, 0])*a))
        final_label = final_label + bb_2.atom_labels

        # Adiciona o bloco 2 na origem da célula unitária (sítio B3)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', 225, degrees=True).as_matrix()) + np.array([0.75, 0.75, 0])*a))
        final_label += bb_2.atom_labels

        # Adiciona o bloco 2 na origem da célula unitária (sítio B4)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', 315, degrees=True).as_matrix()) + np.array([0.25, 0.75, 0])*a))
        final_label += bb_2.atom_labels

        # Changes the X atoms by the desired bond_atom
        final_label, final_pos = Tools.change_X_atoms(final_label, final_pos, bond_atom)

        # Cria a estrutura como entidade do pymatgen
        struct = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)

        # Remove os átomos duplicados
        struct.sort(reverse=True)
        struct.merge_sites(tol=.7, mode='delete')
        struct.translate_sites(range(len(struct.as_dict()['sites'])), [0, 0, 0.5], frac_coords=True, to_unit_cell=True)

        # Simetriza a estrutura
        symm = SpacegroupAnalyzer(struct, symprec=self.symm_tol, angle_tolerance=self.angle_tol)
        struct_symm_prim = symm.get_refined_structure()

        if stacking not in self.available_stacking['SQL_A']:
            raise Exception(f"""{stacking} is not in the available stack list for SQL-A net.
    Available options are: {self.available_stacking['SQL_A']}""")

        if stacking == 'AA':
            self.stacking = 'AA'
            self.symm_structure = struct_symm_prim

        if stacking == 'AB1':
            self.stacking = 'AB1'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal*(1, 1, 0.5)
            B = Tools.translate_inside(ion_conv_crystal*(1, 1, 1.5) + (1/4, 1/4, 0))

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_1 = Structure(lattice, AB_label + AB_label, AB, coords_are_cartesian=False)
            AB_1_symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = AB_1_symm.get_refined_structure()

        if stacking == 'AB2':
            self.stacking = 'AB2'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal*(1, 1, 0.5)
            B = Tools.translate_inside(ion_conv_crystal*(1, 1, 1.5) + (1/2, 0, 0))

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_2 = Structure(lattice, AB_label + AB_label, AB, coords_are_cartesian=False)
            AB_2_symm = SpacegroupAnalyzer(AB_2, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = AB_2_symm.get_refined_structure()

        # Create AAl stacking. Tetragonal cell with two sheets per cell shifited by the shift_vector in angstroms. 
        if stacking == 'AAl':
            self.stacking = 'AAl'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 2)

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal * np.array([1, 1, 0.5])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(*Tools.cell_to_cellpar(cell))
            shift_vector = np.dot(r, np.array(shift_vector))

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal * np.array([1, 1, 1.5]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            self.symm_structure = Structure(lattice, AB_label+AB_label, AB, coords_are_cartesian=False)

        # Create AA tilted stacking. Tilted tetragonal cell with two sheets per cell tilted by tilt_angle. 
        if stacking == 'AAt':
            self.stacking = 'AAt'

            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz'] for i in struct_symm_prim.as_dict()['sites']])
            cell = struct_symm_prim.as_dict()['lattice']

            # Shift the cell by the tilt angle
            print(cell)
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta= cell['beta'] - tilt_angle
            gamma = cell['gamma']

            new_cell = Tools.cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + np.array([0, 0, 0.5*c_cell]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(new_cell)
            self.symm_structure = Structure(lattice, AB_label+AB_label, AB, coords_are_cartesian=True)

        if stacking == 'ABC1':
            self.stacking = 'ABC1'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(ion_conv_crystal*(1, 1, 1) + (2/3, 1/3, 0))
            C = Tools.translate_inside(ion_conv_crystal*(1, 1, 1/3) + (4/3, 2/3, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label+ABC_label, ABC, coords_are_cartesian=False)
            ABC_f_symm = SpacegroupAnalyzer(ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        if stacking == 'ABC2':
            self.stacking = 'ABC2'
            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(ion_conv_crystal*(1, 1, 1) + (1/3, 0, 0))
            C = Tools.translate_inside(ion_conv_crystal*(1, 1, 1/3) + (2/3, 0, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label+ABC_label, ABC, coords_are_cartesian=False)
            ABC_f_symm = SpacegroupAnalyzer(ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        dict_structure = self.symm_structure.as_dict()

        self.lattice = dict_structure['lattice']['matrix']

        self.atom_labels = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.n_atoms = len(self.symm_structure)
        self.composition = self.symm_structure.formula

        if self.verbosity is True:
            print(self.symm_structure)

        # Get the simmetry information of the generated structure
        self.lattice_type = symm.get_lattice_type()
        self.space_group = symm.get_space_group_symbol()
        self.space_group_n = symm.get_space_group_number()

        symm_op = symm.get_point_group_operations()
        self.hall = symm.get_hall()

        if print_result == True:
            Tools.print_result(self.name, str(self.lattice_type), str(self.hall[0:2]), str(self.space_group), str(self.space_group_n), len(symm_op))

        return [self.name, str(self.lattice_type), str(self.hall[0:2]), str(self.space_group), str(self.space_group_n), len(symm_op)]

    def create_kgd_structure(self,
                             name_bb_a,
                             name_bb_b,
                             stacking='AA',
                             bond_atom='N',
                             c_parameter_base=3.6,
                             print_result=True,
                             slab=10,
                             shift_vector=[1, 1, 0],
                             tilt_angle=5):

        '''Creates a COF with KGD network'''

        self.topology = 'KGD'
        self.dimension = 2

        bb_1 = Building_Block(name_bb_a, self.lib_bb)
        bb_2 = Building_Block(name_bb_b, self.lib_bb)

        self.name = f'{bb_1.name}-{bb_2.name}-KGD-{stacking}'
        self.charge = bb_1.charge + bb_2.charge
        self.chirality = bb_1.chirality or bb_2.chirality

        if self.verbosity is True:
            print('Starting the creation of a KGM net')

        if bb_1.connectivity != 6:
            print('Building block A must present connectivity 6 insted of', len(bb_1.connectivity))
            return None
        if bb_2.connectivity != 3:
            print('Building block B must present connectivity 3 insted of', len(bb_2.connectivity))
            return None

        # Calculate the cell parameter based on the size of the building blocks
        size_a = max(bb_1.size)
        if self.verbosity:
            print('BB_A size:', size_a)

        size_b = max(bb_2.size)
        if self.verbosity:
            print('BB_B size:', size_b)

        # Defines the cell parameter a
        a = 2*np.cos(np.radians(30))*(size_a + size_b) # np.cos(np.radians(30))*2*(size_a + size_b)

        if self.verbosity:
            print('Calculated cell parameter a:', a)

        # Gets the maximum distace in the z axis to create the c parameter
        delta_a = abs(max(np.transpose(bb_1.atom_pos)[2])) + abs(min(np.transpose(bb_1.atom_pos)[2]))
        delta_b = abs(max(np.transpose(bb_2.atom_pos)[2])) + abs(min(np.transpose(bb_2.atom_pos)[2]))

        # Build the matrix of unitary hexagonal cell
        if stacking == 'A':
            c = max([slab, c_parameter_base + max([delta_a, delta_b])])
        else:
            c = c_parameter_base + max([delta_a, delta_b])

        # Define the cell lattice
        lattice = [[a, 0, 0], 
                   [-0.5*a, np.sqrt(3)/2*a, 0], 
                   [0, 0, c]]

        if self.verbosity is True:
            print('Unitary cell built:', lattice)

        # Add tbe building block 1 (C6) on the origin of the unitary cell (A1 site)
        final_pos = bb_1.atom_pos
        final_label = list(bb_1.atom_labels)

        # Add tbe building block 2 (C3) on [0.5, 0.0, 0.0] of the unitary cell (B1 site)
        #final_pos = np.dot(bb_2.atom_pos, R.from_euler('z', -180, degrees=True).as_matrix()) + np.array([0, 0.57735027, 0])*a
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', -180, degrees=True).as_matrix()) + np.array([0, 0.57735027, 0])*a))
        final_label += list(bb_2.atom_labels)

        # Add tbe building block 2 (C3) on [0.0, 0.5, 0.0] of the unitary cell (B2 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', 0, degrees=True).as_matrix()) + np.array([0.5, 0.28867513, 0])*a))
        final_label += list(bb_2.atom_labels)

        # Changes the X atoms by the desired bond_atom
        final_label, final_pos = Tools.change_X_atoms(final_label, final_pos, bond_atom)
        
        # Cria a estrutura como entidade do pymatgen
        struct = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)

        # Remove os átomos duplicados
        struct.sort(reverse=True)
        struct.merge_sites(tol=.7, mode='delete')
        struct.translate_sites(range(len(struct.as_dict()['sites'])), [0, 0, 0.5], frac_coords=True, to_unit_cell=True)
        
        atom_labels = np.array([[i['label']] for i in struct.as_dict()['sites']]).flatten()
        atom_pos = np.array([i['xyz'] for i in struct.as_dict()['sites']])
        cell = np.array(struct.as_dict()['lattice']['matrix'])
        
        # Change the X atoms by the desired bond_atom
        final_label, final_pos = Tools.change_X_atoms(final_label, final_pos, bond_atom)

        # Cria a estrutura como entidade do pymatgen
        struct = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)
              
        # Remove os átomos duplicados
        struct.sort(reverse=True)
        struct.merge_sites(tol=.7, mode='delete')
        struct.translate_sites(range(len(struct.as_dict()['sites'])), [0, 0, 0.5], frac_coords=True, to_unit_cell=True)
        
        atom_labels = np.array([[i['label']] for i in struct.as_dict()['sites']]).flatten()
        atom_pos = np.array([i['xyz'] for i in struct.as_dict()['sites']])
        cell = np.array(struct.as_dict()['lattice']['matrix'])
    
        struct_symm_prim = struct
        
        if stacking not in self.available_stacking[self.topology]:
                raise Exception(f"""{stacking} is not in the available stack list for HCB net.
    Available options are: {self.available_stacking[self.topology]}""")
        
        # Create A stacking, a 2D isolated sheet with slab
        if stacking == 'A':
            self.stacking = 'A'
            # Simetriza a estrutura
            symm = SpacegroupAnalyzer(struct, symprec=self.symm_tol, angle_tolerance=self.angle_tol)
            self.symm_structure = symm.get_refined_structure()
        
        # Create AA staking. By default one sheet per unitary cell is used
        if stacking == 'AA':
            self.stacking = 'AA'
            # Simetriza a estrutura
            symm = SpacegroupAnalyzer(struct, symprec=self.symm_tol, angle_tolerance=self.angle_tol)
            self.symm_structure = symm.get_refined_structure()

        # Create AB1 staking. 
        if stacking == 'AB1':
            self.stacking = 'AB1'
            labels_conv_crystal = np.array([[i['label']] for i in struct.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz'] for i in struct.as_dict()['sites']])
            cell = np.array(struct.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal
            B = ion_conv_crystal

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_1 = Structure(lattice, AB_label + AB_label, AB, coords_are_cartesian=True)

            B_list = np.arange(len(AB_label)) + len(AB_label)

            AB_1.translate_sites(B_list, [2/3, 1/3, 0.5], frac_coords=True, to_unit_cell=True)
            symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        # Create AB2 stacking
        if stacking == 'AB2':
            self.stacking = 'AB2'
            labels_conv_crystal = np.array([[i['label']] for i in struct.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz'] for i in struct.as_dict()['sites']])
            cell = np.array(struct.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal
            B = ion_conv_crystal

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_2 = Structure(lattice, AB_label + AB_label, AB, coords_are_cartesian=True)

            B_list = np.arange(len(AB_label)) + len(AB_label)
            AB_2.translate_sites(B_list, [1/2, 0, 0.5], frac_coords=True, to_unit_cell=True)
            symm = SpacegroupAnalyzer(AB_2, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()
            
        # Create AAl stacking. Hexagonal cell with two sheets per cell shifited by the shift_vector in angstroms. 
        if stacking == 'AAl':
            self.stacking = 'AAl'
            labels_conv_crystal = np.array([[i['label']] for i in struct.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct.as_dict()['sites']])
            cell = np.array(struct.as_dict()['lattice']['matrix'])*(1, 1, 2)

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal * np.array([1, 1, 0.5])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(*Tools.cell_to_cellpar(cell))
            shift_vector = np.dot(r, np.array(shift_vector))

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal * np.array([1, 1, 1.5]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            symm = Structure(lattice, AB_label+AB_label, AB, coords_are_cartesian=False)
            
            self.symm_structure = symm.get_refined_structure()

        # Create AA tilted stacking. 
        #Tilted Hexagonal cell by tilt_angle with two sheets per cell shifited by the shift_vector in angstroms.
        if stacking == 'AAt':
            self.stacking = 'AAt'

            labels_conv_crystal = np.array([[i['label']] for i in struct.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz'] for i in struct.as_dict()['sites']])
            cell = struct.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta= cell['beta'] - tilt_angle
            gamma = cell['gamma']

            new_cell = Tools.cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + np.array([0, 0, 0.5*c_cell]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(new_cell)
            symm = Structure(lattice, AB_label+AB_label, AB, coords_are_cartesian=True)
            
            self.symm_structure = symm.get_refined_structure()

        if stacking == 'ABC1':
            self.stacking = 'ABC1'
            labels_conv_crystal = np.array([[i['label']] for i in struct.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct.as_dict()['sites']])
            cell = np.array(struct.as_dict()['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(ion_conv_crystal*(1, 1, 1) + (2/3, 1/3, 0))
            C = Tools.translate_inside(ion_conv_crystal*(1, 1, 1/3) + (4/3, 2/3, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label+ABC_label, ABC, coords_are_cartesian=False)
            symm = SpacegroupAnalyzer(ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_primitive_standard_structure()

        if stacking == 'ABC2':
            self.stacking = 'ABC2'
            labels_conv_crystal = np.array([[i['label']] for i in struct['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in struct['sites']])
            cell = np.array(struct['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(ion_conv_crystal*(1, 1, 1) + (1/3, 0, 0))
            C = Tools.translate_inside(ion_conv_crystal*(1, 1, 1/3) + (2/3, 0, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label+ABC_label, ABC, coords_are_cartesian=False)
            symm = SpacegroupAnalyzer(ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_primitive_standard_structure()
            
        dict_structure = self.symm_structure.as_dict()

        self.lattice = dict_structure['lattice']['matrix']

        self.atom_labels = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.n_atoms = len(self.symm_structure)
        self.composition = self.symm_structure.formula

        if self.verbosity is True:
            print(self.symm_structure)

        # Get the simmetry information of the generated structure
        self.lattice_type = symm.get_lattice_type()
        self.space_group = symm.get_space_group_symbol()
        self.space_group_n = symm.get_space_group_number()

        symm_op = symm.get_point_group_operations()
        self.hall = symm.get_hall()

        if print_result == True:
            Tools.print_result(self.name, str(self.lattice_type), str(self.hall[0:2]), str(self.space_group), str(self.space_group_n), len(symm_op))

        return [self.name, str(self.lattice_type), str(self.hall[0:2]), str(self.space_group), str(self.space_group_n), len(symm_op)]

################   Save file in different formats  ############################

    def save_cif(self, supercell=False, path=None):

        if path is not None:
            self.out_path = path

        os.makedirs(self.out_path, exist_ok=True)

        if supercell is not False:
            self.symm_structure.make_supercell([[supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, supercell[2]]])

        self.symm_structure.to(filename=os.path.join(self.out_path, self.name + '.cif'))
    
    def save_json(self, supercell=False, path=None):

        if path is not None:
            self.out_path = path
        
        os.makedirs(self.out_path, exist_ok=True)

        if supercell is not False:
            self.symm_structure.make_supercell([[supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, supercell[2]]])

        dict_sctructure = self.symm_structure.as_dict()

        a, b, c = dict_sctructure['lattice']['a'], dict_sctructure['lattice']['b'], dict_sctructure['lattice']['c']

        alpha = round(dict_sctructure['lattice']['alpha'], 3)
        beta = round(dict_sctructure['lattice']['beta'], 3)
        gamma = round(dict_sctructure['lattice']['gamma'], 3)

        atom_labels = [i['label'] for i in dict_sctructure['sites']]

        atom_pos = [i['xyz'] for i in dict_sctructure['sites']]

        Tools.save_json(self.out_path, self.name, [a, b, c, alpha, beta, gamma], atom_labels, atom_pos)

    def save_xsf(self, supercell=False, path=None):

        if path is not None:
            self.out_path = path
       
        os.makedirs(self.out_path, exist_ok=True)

        if supercell is not False:
            self.symm_structure.make_supercell([[supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, supercell[2]]])

        dict_sctructure = self.symm_structure.as_dict()

        a, b, c = dict_sctructure['lattice']['a'], dict_sctructure['lattice']['b'], dict_sctructure['lattice']['c']

        alpha = round(dict_sctructure['lattice']['alpha'], 3)
        beta = round(dict_sctructure['lattice']['beta'], 3)
        gamma = round(dict_sctructure['lattice']['gamma'], 3)

        atom_labels = [i['label'] for i in dict_sctructure['sites']]

        atom_pos = [i['xyz'] for i in dict_sctructure['sites']]

        Tools.save_xsf(self.out_path, self.name, [a, b, c, alpha, beta, gamma], atom_labels, atom_pos)
        
    def save_pdb(self, supercell=False, path=None):

        if path is not None:
            self.out_path = path

        os.makedirs(self.out_path, exist_ok=True)

        if supercell is not False:
            self.symm_structure.make_supercell([[supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, supercell[2]]])

        dict_sctructure = self.symm_structure.as_dict()

        a, b, c = dict_sctructure['lattice']['a'], dict_sctructure['lattice']['b'], dict_sctructure['lattice']['c']

        alpha = round(dict_sctructure['lattice']['alpha'], 3)
        beta = round(dict_sctructure['lattice']['beta'], 3)
        gamma = round(dict_sctructure['lattice']['gamma'], 3)

        atom_labels = [i['label'] for i in dict_sctructure['sites']]

        atom_pos = [i['xyz'] for i in dict_sctructure['sites']]

        Tools.save_pdb(self.out_path, self.name, [a, b, c, alpha, beta, gamma], atom_labels, atom_pos)
        
    def save_vasp(self, supercell=False, path=None):

        if path is not None:
            self.out_path = path
        
        os.makedirs(self.out_path, exist_ok=True)

        if supercell is not False:
            self.symm_structure.make_supercell([[supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, supercell[2]]])

        self.symm_structure.to(fmt='poscar', filename=os.path.join(self.out_path, self.name + '.vasp'))

    def save_turbomole(self, supercell=[1,1,2], path=None):

        if path is not None:
            self.out_path = path
        
        os.makedirs(self.out_path, exist_ok=True)

        if supercell is not False:
            self.symm_structure.make_supercell([[supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, supercell[2]]])

        dict_sctructure = self.symm_structure.as_dict()

        a, b, c = dict_sctructure['lattice']['a'], dict_sctructure['lattice']['b'], dict_sctructure['lattice']['c']

        alpha = round(dict_sctructure['lattice']['alpha'], 3)
        beta = round(dict_sctructure['lattice']['beta'], 3)
        gamma = round(dict_sctructure['lattice']['gamma'], 3)

        atom_labels = [i['label'] for i in dict_sctructure['sites']]

        atom_pos = [i['xyz'] for i in dict_sctructure['sites']]

        temp_file = open(os.path.join(self.out_path ,self.name + '.coord'), 'w')
        temp_file.write('$coord angs\n')

        for i in range(len(atom_labels)):
            temp_file.write('{:>15.7f}{:>15.7f}{:>15.7f}   {:<5s}\n'.format(atom_pos[i][0], atom_pos[i][1], atom_pos[i][2], atom_labels[i]))

        temp_file.write('$periodic 3\n')
        temp_file.write('$cell\n')
        temp_file.write(f'{a}  {b}  {c}  {alpha}  {beta}  {gamma}\n')
        temp_file.write('$opt\n')
        temp_file.write('   engine=inertial\n')
        temp_file.write('$end\n')

        temp_file.close()

    def save_xyz(self, supercell=False, path=None):

        if path is not None:
            self.out_path = path
        
        os.makedirs(self.out_path, exist_ok=True)

        if supercell is not False:
            self.symm_structure.make_supercell([[supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, supercell[2]]])

        dict_sctructure = self.symm_structure.as_dict()

        a, b, c = dict_sctructure['lattice']['a'], dict_sctructure['lattice']['b'], dict_sctructure['lattice']['c']

        alpha = round(dict_sctructure['lattice']['alpha'], 3)
        beta = round(dict_sctructure['lattice']['beta'], 3)
        gamma = round(dict_sctructure['lattice']['gamma'], 3)

        atom_labels = [i['label'] for i in dict_sctructure['sites']]

        atom_pos = [i['xyz'] for i in dict_sctructure['sites']]

        temp_file = open(os.path.join(self.out_path, self.name + '.xyz'), 'w')
        temp_file.write(f'{len(atom_labels)} \n')

        temp_file.write(f'{a}  {b}  {c}  {alpha}  {beta}  {gamma}\n')

        for i in range(len(atom_labels)):
            temp_file.write('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atom_labels[i], atom_pos[i][0], atom_pos[i][1], atom_pos[i][2]))

        temp_file.close()

    def save_qe(self, supercell=False, angs=False, path=None, ecut=40, erho=360, k_dist=0.3):

        if path is not None:
            self.out_path = path
        
        os.makedirs(self.out_path, exist_ok=True)

        if supercell is not False:
            self.symm_structure.make_supercell([[supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, supercell[2]]])

        dict_sctructure = self.symm_structure.as_dict()

        cell = dict_sctructure['lattice']['matrix']

        a = dict_sctructure['lattice']['a']
        b = dict_sctructure['lattice']['b']
        c = dict_sctructure['lattice']['c']

        celldm1 = a*1.8897259886  # 1 angstrom = 1.8897259886 bohr
        celldm2 = b/a
        celldm3 = c/a
        celldm4 = Tools.angle(cell[0], cell[1], unit='cos')
        celldm5 = Tools.angle(cell[0], cell[2], unit='cos')
        cellcm6 = Tools.angle(cell[1], cell[2], unit='cos')

        kx, ky, kz = Tools.get_kgrid(Tools.cellpar_to_cell(cell), k_dist)

        ion_conv_crystal = [[i['label']] + i['abc'] for i in dict_sctructure['sites']]
        ion_conv_angstrom = [[i['label']] + i['xyz'] for i in dict_sctructure['sites']]

        if angs is False:
            ion_pos = ion_conv_crystal
        if angs is True:
            ion_pos = ion_conv_angstrom

        out_file = open(os.path.join(self.out_path, self.name + '.in'), 'w', newline='\n')

        out_file.write('#!/bin/sh\n')
        out_file.write('\n')
        out_file.write('#PBS -l walltime=100:00:00\n')
        out_file.write('#PBS -l select=2:ncpus=48:mpiprocs=48\n')
        out_file.write(f'#PBS -N {self.name}\n')
        out_file.write('#PBS -m bea \n')
        out_file.write('#PBS -j oe\n')
        out_file.write('#PBS -V\n')
        out_file.write('\n')
        out_file.write('cd ${PBS_O_WORKDIR}\n')
        out_file.write('\n')
        out_file.write('module load intel/2018.4\n')
        out_file.write('\n')
        out_file.write('PSEUDO_DIR=\"/home/users/lipelopes/pseudo\"\n')
        out_file.write('CMD_PW="mpirun -np 96 /home/users/lipelopes/qe-6.4/bin/pw.x -nk 2"\n')
        out_file.write(f'PREFIX=\'{self.name}\'\n')
        out_file.write('CALC=\'vc-relax\'\n')
        out_file.write(f'SCRATCH_DIR=\'/scratch/31081a/lipelopes/{self.name}/\'\n')
        out_file.write('\n')
        out_file.write('cat>$PREFIX.$CALC.in<<EOF\n')

        out_file.write(' &control\n')
        out_file.write('    calculation = \'$CALC\' ,\n')
        out_file.write('    restart_mode = \'from_scratch\' ,\n')
        out_file.write('    wf_collect = .false. ,\n')
        out_file.write('    outdir = \'$SCRATCH_DIR\' ,\n')
        out_file.write('    pseudo_dir = \'$PSEUDO_DIR\' ,\n')
        out_file.write('    prefix = \'$PREFIX\' ,\n')
        out_file.write('    verbosity =\'high\' ,\n')
        out_file.write('    tstress = .true. ,\n')
        out_file.write('    tprnfor = .true. ,\n')
        out_file.write('    etot_conv_thr = 1.0d-4 ,\n')
        out_file.write('    forc_conv_thr = 1.0d-5 ,\n')
        out_file.write('    nstep=1000 ,\n')
        out_file.write(' / \n')

        out_file.write(' &system\n')

        if self.lattice_type == 'cubic' and 'P' in self.hall:
            out_file.write('    ibrav =   1\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')

        elif self.lattice_type == 'cubic' and 'F' in self.hall:
            out_file.write('    ibrav =   2\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')

        elif self.lattice_type == 'cubic' and 'I' in self.hall:
            out_file.write('    ibrav =   3\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')

        elif self.lattice_type == 'hexagonal':
            out_file.write('    ibrav =   4\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')
            out_file.write(f'    celldm(3) =   {celldm3:.8f}\n')

        elif self.lattice_type == 'trigonal' and 'P' in self.hall:
            out_file.write('    ibrav =   4\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')
            out_file.write(f'    celldm(3) =   {celldm3:.8f}\n')

        elif self.lattice_type == 'trigonal' and 'R' in self.hall:
            out_file.write('    ibrav =   5\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')
            out_file.write(f'    celldm(3) =   {celldm3:.8f}\n')

        elif self.lattice_type == 'tetragonal' and 'P' in self.hall:
            out_file.write('    ibrav =   6\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')
            out_file.write(f'    celldm(3) =   {celldm3:.8f}\n')

        elif self.lattice_type == 'tetragonal' and 'I' in self.hall:
            out_file.write('    ibrav =   7\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')
            out_file.write(f'    celldm(3) =   {celldm3:.8f}\n')

        elif self.lattice_type == 'orthorhombic' and 'P' in self.hall:
            out_file.write('    ibrav =   8\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')
            out_file.write(f'    celldm(2) =   {celldm2:.8f}\n')
            out_file.write(f'    celldm(3) =   {celldm3:.8f}\n')

        elif self.lattice_type == 'orthorhombic' and 'C' or 'A' in self.hall:
            out_file.write('    ibrav =   9\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')
            out_file.write(f'    celldm(2) =   {celldm2:.8f}\n')
            out_file.write(f'    celldm(3) =   {celldm3:.8f}\n')

        elif self.lattice_type == 'orthorhombic' and 'F' in self.hall:
            out_file.write('    ibrav =   10\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')
            out_file.write(f'    celldm(2) =   {celldm2:.8f}\n')
            out_file.write(f'    celldm(3) =   {celldm3:.8f}\n')

        elif self.lattice_type == 'orthorhombic' and 'I' in self.hall:
            out_file.write('    ibrav =   11\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')
            out_file.write(f'    celldm(2) =   {celldm2:.8f}\n')
            out_file.write(f'    celldm(3) =   {celldm3:.8f}\n')

        elif self.lattice_type == 'monoclinic' and round(c, 2) > round(b, 2) and 'P' in self.hall:
            out_file.write('    ibrav =   12\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')
            out_file.write(f'    celldm(2) =   {celldm2:.8f}\n')
            out_file.write(f'    celldm(3) =   {celldm3:.8f}\n')
            out_file.write(f'    celldm(4) =   {celldm4:.8f}\n')

        elif self.lattice_type == 'monoclinic' and round(b, 2) > round(c, 2) and 'P' in self.hall:
            out_file.write('    ibrav =  -12\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')
            out_file.write(f'    celldm(2) =   {celldm2:.8f}\n')
            out_file.write(f'    celldm(3) =   {celldm3:.8f}\n')
            out_file.write(f'    celldm(5) =   {celldm5:.8f}\n')

        elif self.lattice_type == 'monoclinic' and round(c, 2) > round(b, 2) and 'C' in self.hall:
            out_file.write('    ibrav =  13\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')
            out_file.write(f'    celldm(2) =   {celldm2:.8f}\n')
            out_file.write(f'    celldm(3) =   {celldm3:.8f}\n')
            out_file.write(f'    celldm(5) =   {celldm4:.8f}\n')

        elif self.lattice_type == 'monoclinic' and round(b, 2) > round(c, 2) and 'C' in self.hall:
            out_file.write('    ibrav =  -13\n')
            out_file.write(f'    celldm(1) =   {celldm1:.8f}\n')
            out_file.write(f'    celldm(2) =   {celldm2:.8f}\n')
            out_file.write(f'    celldm(3) =   {celldm3:.8f}\n')
            out_file.write(f'    celldm(5) =   {celldm5:.8f}\n')

        elif self.lattice_type == 'triclinic':
            out_file.write('    ibrav =  0\n')

        out_file.write(f'    nat =     {self.n_atoms}\n')
        out_file.write(f'    ntyp =     {len(set(self.atom_labels))}\n')
        out_file.write(f'    ecutwfc = {ecut} \n')
        out_file.write(f'    ecutrho = {erho} \n')
        out_file.write('    !occupations=\'smearing\' , \n')
        out_file.write('    !degauss=0.001 , \n')
        out_file.write('    !smearing=\'gaussian\' , \n')
        out_file.write('    vdw_corr=\'grimme-d3\' , \n')
        out_file.write(' / \n')
        out_file.write(' &electrons\n')
        out_file.write('    conv_thr =  1.0D-8 ,\n')
        out_file.write('    electron_maxstep = 100 ,\n')
        out_file.write('    mixing_beta = 0.3 ,\n')
        out_file.write(' / \n')

        out_file.write(' &IONS\n')
        out_file.write('    ion_dynamics = \'bfgs\'\n')
        out_file.write(' / \n')

        out_file.write(' &CELL\n')
        out_file.write('    cell_dynamics = \'bfgs\'\n')
        if self.lattice_type == 'triclinic':
            out_file.write('    cell_dofree = \'all\'\n')
        else:
            out_file.write('    cell_dofree = \'ibrav\'\n')
        out_file.write(' / \n')

        out_file.write('ATOMIC_SPECIES\n')
        for atom in set(self.atom_labels):
            out_file.write(f' {atom}   {Tools.elements_dict()[atom]:.4f}  {atom}.PSEUDO.UPF\n')
        out_file.write('\n')

        if self.lattice_type == 'triclinic':
            out_file.write('CELL_PARAMETERS (angstrom) \n')
            for v in cell:
                out_file.write(f'{v[0]:.9f}      {v[0]:.9f}      {v[0]:.9f}\n')
                out_file.write('\n')

        coords_type = 'angstrom'
        if angs is False:
            coords_type = 'crystal'
        out_file.write(f'ATOMIC_POSITIONS ({coords_type})\n')

        for atom in ion_pos:
            out_file.write('{:<5s}{:>15.9f}{:>15.9f}{:>15.9f}\n'.format(atom[0], atom[1], atom[2], atom[3]))

        out_file.write('K_POINTS automatic\n')
        out_file.write(f'   {kx} {ky} {kz}  1 1 1\n')

        out_file.write('EOF\n')
        out_file.write('$CMD_PW < $PREFIX.$CALC.in > $PREFIX.$CALC.out')

        out_file.close()
