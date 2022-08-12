# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020

@author: Felipe Lopes de Oliveira
"""

import os
import numpy as np
import pycofbuilder.tools as Tools


class Building_Block():

    def __init__(self, name=None, verbosity=False, save_dir=False):

        _ROOT = os.path.abspath(os.path.dirname(__file__))

        self.name = name
        self.verbosity = verbosity
        self.main_path = os.path.join(_ROOT, 'data')
        self.connectivity = None
        self.simetry = None
        self.size = 0
        self.mass = None
        self.composition = None
        self.charge = 0
        self.multiplicity = 1
        self.polarity = 0
        self.chirality = False
        self.atom_labels = None
        self.atom_pos = None

        # Check if save_dir exists and try to create it if not
        self.save_dir = save_dir
        if not os.path.exists(self.save_dir):
            os.makedirs(self.save_dir)

        if self.name is not None:
            self.get_BB()

    def get_BB(self):
        '''Automatically read or create a buiding block based on its name'''
        simm_check, nucleo_check, conector_check, radicals_check = self.check_existence()

        if simm_check and nucleo_check and conector_check and radicals_check:
            if self.name + '.xyz' in os.listdir(self.save_dir):
                self.read_structure()

            else:
                BB_name = self.name.split('_')
                simmetry = BB_name[0]
                nucleo = BB_name[1]
                conector = BB_name[2]
                radicals = BB_name[3:] + ['H']*(6 - len(BB_name[3:]))

                if simmetry == 'C2':
                    self.create_C2_BB(nucleo,
                                      conector,
                                      radicals[0],
                                      radicals[1],
                                      radicals[2],
                                      radicals[3],
                                      radicals[4],
                                      radicals[5])
                    self.save()
                if simmetry == 'C3':
                    self.create_C3_BB(nucleo,
                                      conector,
                                      radicals[0],
                                      radicals[1],
                                      radicals[2],
                                      radicals[3],
                                      radicals[4],
                                      radicals[5])
                    self.save()
                if simmetry == 'C4':
                    self.create_C4_BB(nucleo,
                                      conector,
                                      radicals[0],
                                      radicals[1],
                                      radicals[2],
                                      radicals[3],
                                      radicals[4],
                                      radicals[5])
                    self.save()
                if simmetry == 'C6':
                    self.create_C6_BB(nucleo,
                                      conector,
                                      radicals[0],
                                      radicals[1],
                                      radicals[2],
                                      radicals[3],
                                      radicals[4],
                                      radicals[5])
                    self.save()

    def n_atoms(self):
        ''' Returns the number of atoms in the unitary cell'''
        return len(self.atom_labels)

    def centralize_molecule(self, by_X=False):
        ''' Centralize the molecule on its geometrical center'''

        transposed = np.transpose(self.atom_pos)
        if by_X is True:
            x_transposed = np.transpose(self.get_X_points()[1])
        if by_X is False:
            x_transposed = transposed
        cm_x = transposed[0] - np.average(x_transposed[0])
        cm_y = transposed[1] - np.average(x_transposed[1])
        cm_z = transposed[2] - np.average(x_transposed[2])

        self.atom_pos = np.transpose([cm_x, cm_y, cm_z])
        return np.transpose([cm_x, cm_y, cm_z])

    def get_X_points(self):
        '''Get the X points in a molecule'''

        if 'X' in self.atom_labels:
            X_labels, X_pos = [], []
            for i in range(len(self.atom_labels)):
                if self.atom_labels[i] == 'X':
                    X_labels += [self.atom_labels[i]]
                    X_pos += [self.atom_pos[i]]

            return X_labels, np.array(X_pos)

        else:
            print('No X ponts could be found!')
            return self.atom_labels, self.atom_pos

    def get_Q_points(self, atom_labels, atom_pos):
        '''Get the Q points in a molecule'''

        Q_labels, Q_pos = [], []

        for i in range(len(atom_labels)):
            if atom_labels[i] == 'Q':
                Q_labels += [atom_labels[i]]
                Q_pos += [atom_pos[i]]

        return Q_labels, np.array(Q_pos)

    def get_R_points(self, atom_labels, atom_pos):
        """
        Get the R points in a molecule
        """

        # Create a dict with the R points
        R_dict = {'R': [],
                  'R1': [],
                  'R2': [],
                  'R3': [],
                  'R4': [],
                  'R5': [],
                  'R6': [],
                  'R7': [],
                  'R8': [],
                  'R9': []}

        # Iterate over the R keys
        for key in R_dict.keys():
            # Iterate over the atoms
            for i in range(len(atom_labels)):
                if atom_labels[i] == key:
                    R_dict[key] += [atom_pos[i]]

        return R_dict

    def add_X(self, label, pos, X='N'):

        label, pos = label, pos

        for i in range(len(label)):
            if label[i] == 'X':
                label[i] = X
        return label, pos

    def calculate_size(self):
        '''Calculate the size of the building block'''
        _, X_pos = self.get_X_points()
        self.size = [np.linalg.norm(i) for i in X_pos]

    def align_to(self, vec=[0, 1, 0]):
        '''Align the molecule to a given vector'''
        _, X_pos = self.get_X_points()
        R_matrix = Tools.rotation_matrix_from_vectors(X_pos[0], vec)

        self.atom_pos = np.dot(self.atom_pos, np.transpose(R_matrix))

    def rotate_to_xy_plane(self):
        '''Rotate the molecule to the xy plane'''
        _, X_pos = self.get_X_points()

        if len(X_pos) == 3:

            normal = np.cross(X_pos[0], X_pos[-1])
            if normal[0] != 0 and normal[1] != 0:
                R_matrix = Tools.rotation_matrix_from_vectors(normal, [0, 0, 1])
                self.atom_pos = np.dot(self.atom_pos, np.transpose(R_matrix))

        if len(X_pos) == 2:
            normal = np.cross(X_pos[0], self.atom_pos[1])
            if normal[0] != 0 and normal[1] != 0:
                R_matrix = Tools.rotation_matrix_from_vectors(normal, [0, 0, 1])
                self.atom_pos = np.dot(self.atom_pos, np.transpose(R_matrix))

    def print_structure(self):
        """
        Print the structure in the form:
        `atom_label     pos_x    pos_y    pos_z`
                """

        for i, _ in enumerate(self.atom_labels):
            print('{:<5s}{:>10.7f}{:>15.7f}{:>15.7f}'.format(self.atom_labels[i],
                                                             self.atom_pos[i][0],
                                                             self.atom_pos[i][1],
                                                             self.atom_pos[i][2]))

    def add_connection_group(self, conector_name):
        '''Adds the functional group by which the COF will be formed from the building blocks'''

        conector_label, conector_pos = Tools.read_gjf_file(os.path.join(self.main_path, 'conector'),
                                                           conector_name)

        # Get the position of the Q points in the structure
        location_Q_struct = self.get_Q_points(self.atom_labels, self.atom_pos)

        for i in range(len(location_Q_struct[0])):
            n_conector_label = conector_label.copy()
            n_conector_pos = conector_pos.copy()

            try:
                # Get the position of the closest atom to Q in the structure
                close_Q_struct = Tools.closest_atom('Q', location_Q_struct[1][i], self.atom_labels, self.atom_pos)[1]
            except Exception:
                # Set the closest position to origin, for building blocks as HDZ that only has two atoms
                close_Q_struct = [0, 0, 0]

            # Get the position of Q in the conection group
            location_Q_connector = self.get_Q_points(n_conector_label, n_conector_pos)  

            # Get the position of the closest atom to Q in the conection group
            close_Q_connector = Tools.closest_atom('Q', location_Q_connector[1][0], n_conector_label, n_conector_pos)[1]

            v1 = close_Q_struct - location_Q_struct[1][i]  # Create the vector Q in the structure
            v2 = np.array(close_Q_connector) - np.array(location_Q_connector[1][0])  # Create the vector Q in the conector

            Rot_m = Tools.rotation_matrix_from_vectors(v2, v1)  # Find the rotation matrix that align v2 with v1

            # Delete the "Q" atom position of the conector group and the structure
            n_conector_pos = np.delete(n_conector_pos, Tools.find_index(np.array([0., 0., 0.]), n_conector_pos), axis=0)

            self.atom_pos = np.delete(self.atom_pos, Tools.find_index(location_Q_struct[1][i], self.atom_pos), axis=0)

            # Rotate and translade the conector group to Q position in the strucutre
            rotated_translated_group = np.dot(n_conector_pos, -np.transpose(Rot_m)) + location_Q_struct[1][i]

            # Add the position of conector atoms to the main structure
            self.atom_pos = np.append(self.atom_pos, rotated_translated_group, axis=0)

            # Remove the Q atoms from structure
            self.atom_labels.remove('Q')
            n_conector_label.remove('Q')

            self.atom_labels = self.atom_labels + n_conector_label

    def add_R_group(self, R_name, R_type):
        '''Adds group R in building blocks'''
        
        # Read the R group
        group_label, group_pos = Tools.read_gjf_file(os.path.join(self.main_path, 'radical'), R_name)

        # Get the position of the R points in the structure
        location_R_struct = self.get_R_points(self.atom_labels, self.atom_pos)[R_type]  

        # Get the position of the R points in the R group
        for i, _ in enumerate(location_R_struct):
            n_group_label = group_label.copy()
            n_group_pos = group_pos.copy()

            # Get the position of the closest atom to R in the structure
            close_R_struct = Tools.closest_atom_struc(R_type,
                                                      location_R_struct[i],
                                                      self.atom_labels,
                                                      self.atom_pos)[1]

            # Get the position of R in the R group
            pos_R_group = self.get_R_points(n_group_label, n_group_pos)['R']

            # Get the position of the closest atom to R in the R group 
            close_R_group = Tools.closest_atom('R', pos_R_group[0], n_group_label, n_group_pos)[1]

            # Create the vector R in the structure
            v1 = close_R_struct - location_R_struct[i]

            # Create the vector R in the R group
            v2 = np.array(close_R_group) - np.array(pos_R_group[0])

            # Find the rotation matrix that align v2 with v1
            Rot_m = Tools.rotation_matrix_from_vectors(v2, v1)  

            # Delete the "R" atom position of the R group and the structure
            n_group_pos = np.delete(n_group_pos, Tools.find_index(np.array([0.0, 0.0, 0.0]), 
                                                                  n_group_pos),
                                                                  axis=0)

            # Rotate and translade the R group to R position in the strucutre
            rotated_translated_group = np.dot(n_group_pos, - np.transpose(Rot_m)) + location_R_struct[i]
            
            # Remove the R atoms from structure
            self.atom_pos = np.delete(self.atom_pos, Tools.find_index(location_R_struct[i], 
                                                                      self.atom_pos), 
                                                                      axis=0)

            # Add the position of rotated atoms to the main structure
            self.atom_pos = np.append(self.atom_pos, rotated_translated_group, axis=0)

            # Remove the R atoms from structure
            self.atom_labels.remove(R_type)

            # Remove the R atoms from R group
            n_group_label.remove('R')

            self.atom_labels = self.atom_labels + n_group_label

    def create_C2_BB(self,
                     nucleo_name='BENZ',
                     conector='CHO',
                     R1='H',
                     R2='H',
                     R3='H',
                     R4='H',
                     R5='H',
                     R6='H',
                     R7='H',
                     R8='H',
                     R9='H'):
        '''Create a building block with C2 simmetry'''

        self.name = f'C2_{nucleo_name}_{conector}'

        self.atom_labels, self.atom_pos = Tools.read_gjf_file(os.path.join(self.main_path, 'nucleo', 'C2'), nucleo_name)

        self.centralize_molecule()

        self.add_connection_group(conector)

        R_list_names = [R1, R2, R3, R4, R5, R6, R7, R8, R9]
        R_list_labels = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9']

        for i in range(len(R_list_names)):
            if R_list_labels[i] in self.atom_labels:
                self.add_R_group(R_list_names[i], R_list_labels[i])
                self.name += f'_{R_list_names[i]}'

        self.connectivity = len([i for i in self.atom_labels if 'X' in i])
        self.align_to()
        self.calculate_size()

    def create_C3_BB(self,
                     nucleo_name='BENZ',
                     conector='CHO',
                     R1='H',
                     R2='H',
                     R3='H',
                     R4='H',
                     R5='H',
                     R6='H',
                     R7='H',
                     R8='H',
                     R9='H'):
        '''Create a building block with C3 simmetry'''

        self.name = f'C3_{nucleo_name}_{conector}'

        self.atom_labels, self.atom_pos = Tools.read_gjf_file(os.path.join(self.main_path,
                                                                           'nucleo',
                                                                           'C3'),
                                                              nucleo_name)

        self.centralize_molecule()

        self.add_connection_group(conector)

        R_list_names = [R1, R2, R3, R4, R5, R6, R7, R8, R9]
        R_list_labels = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9']

        for i in range(len(R_list_names)):
            if R_list_labels[i] in self.atom_labels:
                self.add_R_group(R_list_names[i], R_list_labels[i])
                self.name += f'_{R_list_names[i]}'

        self.connectivity =  len([i for i in self.atom_labels if 'X' in i])
        self.align_to()
        self.calculate_size()

    def create_C4_BB(self,
                     nucleo_name='BENZ',
                     conector='CHO',
                     R1='H',
                     R2='H',
                     R3='H',
                     R4='H',
                     R5='H',
                     R6='H',
                     R7='H',
                     R8='H',
                     R9='H'):
        '''Create a building block with C4 simmetry'''

        self.name = f'C4_{nucleo_name}_{conector}'

        self.atom_labels, self.atom_pos = Tools.read_gjf_file(
            os.path.join(self.main_path, 
                         'nucleo',
                         'C4'),
            nucleo_name)

        self.centralize_molecule()

        self.add_connection_group(conector)

        R_list_names = [R1, R2, R3, R4, R5, R6, R7, R8, R9]
        R_list_labels = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9']

        for i in range(len(R_list_names)):
            if R_list_labels[i] in self.atom_labels:
                self.add_R_group(R_list_names[i], R_list_labels[i])
                self.name += f'_{R_list_names[i]}'

        self.connectivity =  len([i for i in self.atom_labels if 'X' in i])
        self.align_to()
        self.calculate_size()

    def create_C6_BB(self,
                     nucleo_name='BENZ',
                     conector='CHO',
                     R1='H',
                     R2='H',
                     R3='H',
                     R4='H',
                     R5='H',
                     R6='H',
                     R7='H',
                     R8='H',
                     R9='H'):
        '''Create a building block with C6 simmetry'''

        self.name = f'C6_{nucleo_name}_{conector}'

        self.atom_labels, self.atom_pos = Tools.read_gjf_file(os.path.join(self.main_path, 'nucleo', 'C6'), nucleo_name)

        self.centralize_molecule()

        self.add_connection_group(conector)

        R_list_names = [R1, R2, R3, R4, R5, R6, R7, R8, R9]
        R_list_labels = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9']

        for i in range(len(R_list_names)):
            if R_list_labels[i] in self.atom_labels:
                self.add_R_group(R_list_names[i], R_list_labels[i])
                self.name += f'_{R_list_names[i]}'

        self.connectivity =  len([i for i in self.atom_labels if 'X' in i])
        self.align_to()
        self.calculate_size()

    def save(self, extension='xyz'):

        if extension == 'xyz':
            Tools.save_xyz(self.save_dir, self.name + '.xyz', self.atom_labels, self.atom_pos)

    def read_structure(self):

        try:
            self.atom_labels, self.atom_pos = Tools.read_xyz_file(self.save_dir, self.name)
            self.connectivity = len([i for i in self.atom_labels if 'X' in i])

            self.align_to()
            self.calculate_size()
        except Exception:
            None
        
    def get_available_nucleo(self):

        C2_list = [i.rstrip('.gjf') for i in os.listdir(os.path.join(self.main_path, 'nucleo', 'C2')) if '.gjf' in i]
        C3_list = [i.rstrip('.gjf') for i in os.listdir(os.path.join(self.main_path, 'nucleo', 'C3')) if '.gjf' in i]
        C4_list = [i.rstrip('.gjf') for i in os.listdir(os.path.join(self.main_path, 'nucleo', 'C4')) if '.gjf' in i]
        C6_list = [i.rstrip('.gjf') for i in os.listdir(os.path.join(self.main_path, 'nucleo', 'C6')) if '.gjf' in i]

        return C2_list, C3_list, C4_list, C6_list
    
    def get_available_R(self):

        R_list = [i.rstrip('.gjf') for i in os.listdir(os.path.join(self.main_path, 'radical')) if '.gjf' in i]

        return R_list

    def get_available_conector(self):

        c_list = [i.rstrip('.gjf') for i in os.listdir(os.path.join(self.main_path, 'conector')) if '.gjf' in i]

        return c_list

    def check_existence(self):

        simm_check = True
        nucleo_check = True
        conector_check = True
        radicals_check = True

        if self.name != None:
            name = self.name.split('_')
            simm = name[0]
            nucleo = name[1]
            conector = name[2]
            radicals = name[3:]

            if simm not in ['C2', 'C3', 'C4', 'C6']:
                print('ERROR!: Building Block simmetry must be C2, C3, C4, or C6.')
                simm_check = False

            if simm == 'C2':
                list = self.get_available_nucleo()[0]
                if nucleo not in list:
                    print(f'ERROR!: {nucleo} not available! Available nucleos with C2 simmetry is {list}')
                    nucleo_check = False
            if simm == 'C3':
                list = self.get_available_nucleo()[1]
                if nucleo not in list:
                    print(f'ERROR!: {nucleo} not available! Available nucleos with C3 simmetry is {list}')
                    nucleo_check = False
            if simm == 'C4':
                list = self.get_available_nucleo()[2]
                if nucleo not in list:
                    print(f'ERROR!: {nucleo} not available! Available nucleos with C4 simmetry is {list}')
                    nucleo_check = False
            if simm == 'C6':
                list = self.get_available_nucleo()[3]
                if nucleo not in list:
                    print(f'ERROR!: {nucleo} not available! Available nucleos with C6 simmetry is {list}')
                    nucleo_check = False
            
            if conector not in self.get_available_conector():
                print(f'ERROR! {conector} is not a available conector. Available list: {self.get_available_conector()}')
                conector_check = False

            radicals_list = self.get_available_R()
            for rad in radicals:
                if rad not in radicals_list:
                    print(f'ERROR! Radical {rad} is not a available radical: Available list{radicals_list}')
                    radicals_check = False

        return simm_check, nucleo_check, conector_check, radicals_check

    def get_bipodal_NH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C2' == i.split('_')[0] and 'NH2' in i.split('_')[2]]

    def get_tripodal_NH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C3' == i.split('_')[0] and 'NH2' in i.split('_')[2]]

    def get_bipodal_CHO(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C2' == i.split('_')[0] and 'CHO' == i.split('_')[2]]

    def get_tripodal_CHO(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C3' == i.split('_')[0] and 'CHO' == i.split('_')[2]]

    def get_bipodal_BOH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C2' == i.split('_')[0] and 'BOH2' == i.split('_')[2]]

    def get_tripodal_BOH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C3' == i.split('_')[0] and 'BOH2' == i.split('_')[2]]

    def get_bipodal_OH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C2' == i.split('_')[0] and 'OH2' == i.split('_')[2]]

    def get_tripodal_OH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C3' == i.split('_')[0] and 'OH2' == i.split('_')[2]]
    
    def get_tetrapodal_squared_OH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C4' == i.split('_')[0] and 'OH2' == i.split('_')[2]]

    def get_tetrapodal_squared_CHO(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C4' == i.split('_')[0] and 'CHO' == i.split('_')[2]]

    def get_tetrapodal_squared_BOH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C4' == i.split('_')[0] and 'BOH2' == i.split('_')[2]]

    def get_tetrapodal_squared_NH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C4' == i.split('_')[0] and 'NH2' in i.split('_')[2]]
