# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020

@author: Felipe Lopes de Oliveira
"""

import os
import numpy as np
import pycofbuilder.tools as Tools


class Building_Block():

    def __init__(self, name=None, verbosity=False, save_dir=False, save_bb=True):

        _ROOTDIR = os.path.abspath(os.path.dirname(__file__))

        self.name = name
        self.verbosity = verbosity
        self.main_path = os.path.join(_ROOTDIR, 'data')
        self.save_dir = save_dir
        self.save_bb = save_bb
        self.connectivity = None
        self.size = 0
        self.mass = None
        self.composition = None
        self.atom_labels = None
        self.atom_pos = None
        self.smiles = None
        self.charge = 0
        self.multiplicity = 1
        self.chirality = None
        self.symmetry = None
        self.core = None
        self.conector = None
        self.radicals = None

        self.available_symmetry = ['C2', 'C3', 'C4', 'C6']

        # Check if save_dir exists and try to create it if not
        if self.save_dir is not False:
            os.makedirs(self.save_dir, exist_ok=True)

        # If a name is provided create building block from this name
        if self.name is not None:
            self.from_name(self.name)

    def __str__(self):
        return self._structure_as_string()

    def __repr__(self):
        return 'BuildingBlock({}, {}, {}, {})'.format(self.symmetry,
                                                      self.core,
                                                      self.conector,
                                                      self.radicals)

    def from_name(self, name):
        '''Automatically read or create a buiding block based on its name'''

        # Check the existence of the building block files
        symm_check, core_check, conector_check, radicals_check = self.check_existence(name)

        error_msg = "Building Block name is invalid!\n"
        error_msg += "Symm: {}, Core: {}, Connector: {}, Radical:{}".format(symm_check,
                                                                            core_check,
                                                                            conector_check,
                                                                            radicals_check)
        assert all([symm_check, core_check, conector_check, radicals_check]), error_msg

        # Make the BB name global
        self.name = name

        BB_name = self.name.split('_')
        self.symmetry = BB_name[0]
        self.core = BB_name[1]
        self.conector = BB_name[2]
        possible_radicals = BB_name[3:] + ['H'] * (9 - len(BB_name[3:]))

        self._create_BB_structure(self.symmetry,
                                  self.core,
                                  self.conector,
                                  *possible_radicals)
        if self.save_bb:
            self.save()

    def from_structure(self, path, file_name):
        '''
        Read the building block structure from a xyz file
        '''

        # Try to read the xyz file
        try:
            self.name = file_name.split('.')[0]

            self.atom_labels, self.atom_pos = Tools.read_xyz_file(path, file_name)
            self.connectivity = len([i for i in self.atom_labels if 'X' in i])

            self._align_to()
            self._calculate_size()

        except Exception:
            print('Error reading the structure!')

        # Try to decompose the building block nomenclature
        try:
            BB_name = self.name.split('_')
            self.symmetry = BB_name[0]
            self.core = BB_name[1]
            self.conector = BB_name[2]
            self.radicals = BB_name[3:]

        except Exception:
            print('Error on the definition of the BB nomenclature.')

    def n_atoms(self):
        ''' Returns the number of atoms in the unitary cell'''
        return len(self.atom_labels)

    def print_structure(self):
        """
        Print the structure in xyz format:
        `atom_label     pos_x    pos_y    pos_z`
        """

        print(self._structure_as_string())

    def _centralize_molecule(self, by_X=True):
        ''' Centralize the molecule on its geometrical center'''

        transposed = np.transpose(self.atom_pos)
        if by_X is True:
            x_transposed = np.transpose(self._get_X_points()[1])
        if by_X is False:
            x_transposed = transposed
        cm_x = transposed[0] - np.average(x_transposed[0])
        cm_y = transposed[1] - np.average(x_transposed[1])
        cm_z = transposed[2] - np.average(x_transposed[2])

        self.atom_pos = np.transpose([cm_x, cm_y, cm_z])
        return np.transpose([cm_x, cm_y, cm_z])

    def _get_X_points(self):
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

    def _get_Q_points(self, atom_labels, atom_pos):
        '''Get the Q points in a molecule'''

        Q_labels, Q_pos = [], []

        for i in range(len(atom_labels)):
            if atom_labels[i] == 'Q':
                Q_labels += [atom_labels[i]]
                Q_pos += [atom_pos[i]]

        return Q_labels, np.array(Q_pos)

    def _get_R_points(self, atom_labels, atom_pos):
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

    def _add_X(self, label, pos, X='N'):

        label, pos = label, pos

        for i in range(len(label)):
            if label[i] == 'X':
                label[i] = X
        return label, pos

    def _calculate_size(self):
        '''Calculate the size of the building block'''
        _, X_pos = self._get_X_points()
        self.size = [np.linalg.norm(i) for i in X_pos]

    def _align_to(self, vec=[0, 1, 0]):
        '''Align the molecule to a given vector'''
        _, X_pos = self._get_X_points()
        R_matrix = Tools.rotation_matrix_from_vectors(X_pos[0], vec)

        self.atom_pos = np.dot(self.atom_pos, np.transpose(R_matrix))

    def _rotate_to_xy_plane(self):
        '''Rotate the molecule to the xy plane'''
        _, X_pos = self._get_X_points()

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

    def _structure_as_string(self):
        struct_string = ''
        for i, _ in enumerate(self.atom_labels):
            struct_string += '{:<5s}{:>10.7f}{:>15.7f}{:>15.7f}\n'.format(self.atom_labels[i],
                                                                          self.atom_pos[i][0],
                                                                          self.atom_pos[i][1],
                                                                          self.atom_pos[i][2])

        return struct_string

    def _add_connection_group(self, conector_name):
        '''Adds the functional group by which the COF will be formed from the building blocks'''

        conector_chem_json = Tools.read_json(
            os.path.join(self.main_path, 'conector'),
            conector_name)

        conector_smiles = conector_chem_json['smiles'].replace('[Q]', '')
        self.smiles = self.smiles.replace('[Q]', f'({conector_smiles})')

        conector_label = conector_chem_json['atoms']['elements']['elementType']
        conector_pos = conector_chem_json['atoms']['coords']['3d']

        # Get the position of the Q points in the structure
        location_Q_struct = self._get_Q_points(self.atom_labels, self.atom_pos)

        for i in range(len(location_Q_struct[0])):

            n_conector_label = conector_label.copy()
            n_conector_pos = conector_pos.copy()

            try:
                # Get the position of the closest atom to Q in the structure
                close_Q_struct = Tools.closest_atom('Q',
                                                    location_Q_struct[1][i],
                                                    self.atom_labels,
                                                    self.atom_pos)[1]
            except Exception:
                # Set the closest position to origin for building blocks as HDZN
                close_Q_struct = [0, 0, 0]

            # Get the position of Q in the conection group
            location_Q_connector = self._get_Q_points(n_conector_label, n_conector_pos)

            # Get the position of the closest atom to Q in the conection group
            close_Q_connector = Tools.closest_atom('Q',
                                                   location_Q_connector[1][0],
                                                   n_conector_label,
                                                   n_conector_pos)[1]

            # Create the vector Q in the structure
            v1 = close_Q_struct - location_Q_struct[1][i]
            # Create the vector Q in the conector
            v2 = np.array(close_Q_connector) - np.array(location_Q_connector[1][0])

            # Find the rotation matrix that align v2 with v1
            Rot_m = Tools.rotation_matrix_from_vectors(v2, v1)

            # Delete the "Q" atom position of the conector group and the structure
            n_conector_pos = np.delete(
                n_conector_pos,
                Tools.find_index(np.array([0., 0., 0.]), n_conector_pos),
                axis=0)

            self.atom_pos = np.delete(
                self.atom_pos,
                Tools.find_index(location_Q_struct[1][i], self.atom_pos),
                axis=0)

            # Rotate and translade the conector group to Q position in the strucutre
            rotated_translated_group = np.dot(
                n_conector_pos,
                -np.transpose(Rot_m)) + location_Q_struct[1][i]

            # Add the position of conector atoms to the main structure
            self.atom_pos = np.append(self.atom_pos, rotated_translated_group, axis=0)

            # Remove the Q atoms from structure
            self.atom_labels.remove('Q')
            n_conector_label.remove('Q')

            self.atom_labels = self.atom_labels + n_conector_label

    def _add_R_group(self, R_name, R_type):
        '''Adds group R in building blocks'''

        # Read the R group
        R_chem_json = Tools.read_json(
            os.path.join(self.main_path, 'radical'),
            R_name)

        r_smiles = R_chem_json['smiles'].replace('[R]', '')
        self.smiles = self.smiles.replace(f'[{R_type}]', r_smiles)

        group_label = R_chem_json['atoms']['elements']['elementType']
        group_pos = np.array(R_chem_json['atoms']['coords']['3d']).astype(float)

        # Get the position of the R points in the structure
        location_R_struct = self._get_R_points(self.atom_labels, self.atom_pos)[R_type]

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
            pos_R_group = self._get_R_points(n_group_label, n_group_pos)['R']

            # Get the position of the closest atom to R in the R group
            close_R_group = Tools.closest_atom('R', pos_R_group[0], n_group_label, n_group_pos)[1]

            # Create the vector R in the structure
            v1 = close_R_struct - location_R_struct[i]

            # Create the vector R in the R group
            v2 = np.array(close_R_group) - np.array(pos_R_group[0])

            # Find the rotation matrix that align v2 with v1
            Rot_m = Tools.rotation_matrix_from_vectors(v2, v1)

            # Delete the "R" atom position of the R group and the structure
            n_group_pos = np.delete(
                n_group_pos,
                Tools.find_index(np.array([0.0, 0.0, 0.0]), n_group_pos),
                axis=0)

            # Rotate and translade the R group to R position in the strucutre
            rotated_translated_group = np.dot(
                n_group_pos,
                -np.transpose(Rot_m)
                ) + location_R_struct[i]

            # Remove the R atoms from structure
            self.atom_pos = np.delete(
                self.atom_pos,
                Tools.find_index(location_R_struct[i], self.atom_pos),
                axis=0)

            # Add the position of rotated atoms to the main structure
            self.atom_pos = np.append(self.atom_pos, rotated_translated_group, axis=0)

            # Remove the R atoms from structure
            self.atom_labels.remove(R_type)

            # Remove the R atoms from R group
            n_group_label.remove('R')

            self.atom_labels = self.atom_labels + n_group_label

    def _create_BB_structure(self,
                             symmetry='C2',
                             core_name='BENZ',
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
        '''Create a building block'''

        self.name = f'{symmetry}_{core_name}_{conector}'

        chem_json = Tools.read_json(
            os.path.join(self.main_path, 'core', symmetry),
            core_name
            )

        self.smiles = chem_json['smiles']

        self.atom_labels = chem_json['atoms']['elements']['elementType']
        self.atom_pos = chem_json['atoms']['coords']['3d']
        self.composition = chem_json['formula']

        self._add_connection_group(conector)

        R_list_names = [R1, R2, R3, R4, R5, R6, R7, R8, R9]
        R_list_labels = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9']

        radical_string = []
        for i in range(len(R_list_names)):
            if R_list_labels[i] in self.atom_labels:
                self._add_R_group(R_list_names[i], R_list_labels[i])
                self.name += f'_{R_list_names[i]}'
                radical_string.append(R_list_names[i])

        self.radicals = radical_string

        self.connectivity = len([i for i in self.atom_labels if 'X' in i])
        self._centralize_molecule()
        self._align_to()
        self._calculate_size()

    def save(self, extension='xyz'):

        if extension == 'xyz':
            Tools.save_xyz(path=self.save_dir,
                           file_name=self.name + '.xyz',
                           atom_labels=self.atom_labels,
                           atom_pos=self.atom_pos)

    def get_available_core(self):
        '''Get the list of available cores'''
        C2_PATH = os.path.join(self.main_path, 'core', 'C2')
        C2_list = [i.rstrip('.gjf') for i in os.listdir(C2_PATH) if '.gjf' in i]

        C3_PATH = os.path.join(self.main_path, 'core', 'C3')
        C3_list = [i.rstrip('.gjf') for i in os.listdir(C3_PATH) if '.gjf' in i]

        C4_PATH = os.path.join(self.main_path, 'core', 'C4')
        C4_list = [i.rstrip('.gjf') for i in os.listdir(C4_PATH) if '.gjf' in i]

        C6_PATH = os.path.join(self.main_path, 'core', 'C6')
        C6_list = [i.rstrip('.gjf') for i in os.listdir(C6_PATH) if '.gjf' in i]

        return C2_list, C3_list, C4_list, C6_list

    def get_available_R(self):
        '''Get the list of available radicals'''
        R_PATH = os.path.join(self.main_path, 'radical')
        R_list = [i.rstrip('.gjf') for i in os.listdir(R_PATH) if '.gjf' in i]

        return R_list

    def get_available_conector(self):
        '''Get the list of available conectores'''
        C_PATH = os.path.join(self.main_path, 'conector')
        C_list = [i.rstrip('.gjf') for i in os.listdir(C_PATH) if '.gjf' in i]

        return C_list

    def check_existence(self, name):

        symm_check = False
        core_check = False
        conector_check = False
        radicals_check = True

        name = name.split('_')
        symm = name[0]
        core = name[1]
        conector = name[2]
        radicals = name[3:]

        BB_dict = {s: self.get_available_core()[i] for i, s in enumerate(self.available_symmetry)}

        if symm in self.available_symmetry:
            symm_check = True
        else:
            print('ERROR!: Building Block symmetry must be C2, C3, C4, or C6.')
            symm_check = False

        if core in BB_dict[symm]:
            core_check = True
        else:
            print(f'ERROR!: {core} not available!')
            print(f'Available cores with {symm} symmetry are {BB_dict[symm]}')

        if conector in self.get_available_conector():
            conector_check = True
        else:
            print(f'ERROR! {conector} is not a available conector.')
            print(f'Available list: {self.get_available_conector()}')

        radicals_list = self.get_available_R()
        for rad in radicals:
            if rad not in radicals_list:
                print(f'ERROR! Radical {rad} is not a available.')
                print(f'Available list: {radicals_list}')
                radicals_check = False

        return symm_check, core_check, conector_check, radicals_check

    def get_bipodal_NH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C2' == i.split('_')[0]
                and 'NH2' in i.split('_')[2]]

    def get_tripodal_NH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C3' == i.split('_')[0]
                and 'NH2' in i.split('_')[2]]

    def get_bipodal_CHO(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C2' == i.split('_')[0]
                and 'CHO' == i.split('_')[2]]

    def get_tripodal_CHO(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C3' == i.split('_')[0]
                and 'CHO' == i.split('_')[2]]

    def get_bipodal_BOH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C2' == i.split('_')[0]
                and 'BOH2' == i.split('_')[2]]

    def get_tripodal_BOH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C3' == i.split('_')[0]
                and 'BOH2' == i.split('_')[2]]

    def get_bipodal_OH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C2' == i.split('_')[0]
                and 'OH2' == i.split('_')[2]]

    def get_tripodal_OH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C3' == i.split('_')[0]
                and 'OH2' == i.split('_')[2]]

    def get_tetrapodal_squared_OH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C4' == i.split('_')[0]
                and 'OH2' == i.split('_')[2]]

    def get_tetrapodal_squared_CHO(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C4' == i.split('_')[0]
                and 'CHO' == i.split('_')[2]]

    def get_tetrapodal_squared_BOH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C4' == i.split('_')[0]
                and 'BOH2' == i.split('_')[2]]

    def get_tetrapodal_squared_NH2(self):

        files_list = os.listdir(self.save_dir)

        return [i.rstrip('.xyz') for i in files_list if 'C4' == i.split('_')[0]
                and 'NH2' in i.split('_')[2]]
