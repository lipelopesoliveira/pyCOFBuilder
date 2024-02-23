# -*- coding: utf-8 -*-
# Created by Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

"""
The BuildingBlock class is used to create the building blocks for the Framework class.
"""

import os
import copy
import numpy as np
from scipy.spatial.transform import Rotation as R

from pycofbuilder.tools import (rotation_matrix_from_vectors,
                                closest_atom,
                                closest_atom_struc,
                                find_index,
                                unit_vector)
from pycofbuilder.io_tools import (read_xyz, save_xyz, read_gjf)
from pycofbuilder.cjson import ChemJSON

from pycofbuilder.logger import create_logger

from pycofbuilder.exceptions import MissingXError


class BuildingBlock():

    def __init__(self, name: str = '', **kwargs):

        self.name: str = name

        self.out_path: str = kwargs.get('out_dir', os.path.join(os.getcwd(), 'out'))
        self.save_bb: bool = kwargs.get('save_bb', True)
        self.bb_out_path: str = kwargs.get('bb_out_path', os.path.join(self.out_path, 'building_blocks'))

        self.logger = create_logger(level=kwargs.get('log_level', 'info'),
                                    format=kwargs.get('log_format', 'simple'),
                                    save_to_file=kwargs.get('save_to_file', False),
                                    log_filename=kwargs.get('log_filename', 'pycofbuilder.log'))

        _ROOTDIR = os.path.abspath(os.path.dirname(__file__))
        self.main_path = os.path.join(_ROOTDIR, 'data')

        self.connectivity = kwargs.get('connectivity', None)
        self.size = kwargs.get('size', None)
        self.mass = kwargs.get('mass', None)
        self.composition = kwargs.get('composition', None)

        self.atom_types = kwargs.get('atom_types', None)
        self.atom_pos = kwargs.get('atom_pos', None)
        self.atom_labels = kwargs.get('atom_labels', None)

        self.smiles = kwargs.get('smiles', None)
        self.charge = kwargs.get('charge', 0)
        self.multiplicity = kwargs.get('multiplicity', 1)
        self.chirality = kwargs.get('chirality', 0)
        self.symmetry = kwargs.get('symmetry', None)

        self.core = kwargs.get('core', None)
        self.conector = kwargs.get('conector', None)
        self.funcGroups = kwargs.get('funcGroups', None)

        self.available_symmetry = ['L2',
                                   'T3',
                                   'S4', 'D4',  # 'R4'
                                   'H6',  # 'O6', 'P6'
                                   # 'C8', 'A8', 'E8'
                                   # 'B12', 'I12', 'U12', 'X12'
                                   ]

        # Check if bb_out_path exists and try to create it if not
        if self.save_bb != '':
            os.makedirs(self.bb_out_path, exist_ok=True)

        # If a name is provided create building block from this name
        if self.name != '':
            if '.' not in self.name:
                self.from_name(self.name)

            # If a name is provided and it is a file, create building block from this file
            if '.' in self.name:
                self.from_file(self.out_path, self.name)

    def __str__(self):
        return self.structure_as_string()

    def __repr__(self):
        return 'BuildingBlock({}, {}, {}, {})'.format(self.symmetry,
                                                      self.core,
                                                      self.conector,
                                                      self.funcGroups)

    def copy(self):
        '''Return a deep copy of the BuildingBlock object'''
        return copy.deepcopy(self)

    def from_file(self, path, file_name):
        '''Read a building block from a file'''
        extension = file_name.split('.')[-1]

        read_func_dict = {'xyz': read_xyz,
                          'gjf': read_gjf}

        self.name = file_name.rstrip(f'.{extension}')
        self.atom_types, self.atom_pos = read_func_dict[extension](path, file_name)
        self.atom_labels = ['C']*len(self.atom_types)

        if any([i == 'X' for i in self.atom_types]):
            self.connectivity = len([i for i in self.atom_types if 'X' in i])
        else:
            raise MissingXError()

        self.centralize(by_X=True)
        self.calculate_size()
        pref_orientation = unit_vector(
            self.get_X_points()[1][0])

        self.align_to(pref_orientation)

    def from_name(self, name):
        '''Automatically read or create a buiding block based on its name'''

        # Check the existence of the building block files
        symm_check, core_check, conector_check, funcGroup_check = self.check_existence(name)

        error_msg = "Building Block name is invalid!\n"
        error_msg += "Symm: {}, Core: {}, Connector: {}, Functional Group:{}".format(symm_check,
                                                                                     core_check,
                                                                                     conector_check,
                                                                                     funcGroup_check)
        assert all([symm_check, core_check, conector_check, funcGroup_check]), error_msg

        self.name = name

        BB_name = self.name.split('_')
        self.symmetry = BB_name[0]
        self.core = BB_name[1]
        self.conector = BB_name[2]
        possible_funcGroups = BB_name[3:] + ['H'] * (9 - len(BB_name[3:]))

        self.create_BB_structure(self.symmetry,
                                 self.core,
                                 self.conector,
                                 *possible_funcGroups)
        if self.save_bb:
            self.save()

    def n_atoms(self):
        ''' Returns the number of atoms in the unitary cell'''
        return len(self.atom_types)

    def print_structure(self):
        """
        Print the structure in xyz format:
        `atom_label     pos_x    pos_y    pos_z`
        """

        print(self.structure_as_string())

    def centralize(self, by_X=True):
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

        if 'X' in self.atom_types:
            X_labels, X_pos = [], []
            for i in range(len(self.atom_types)):
                if self.atom_types[i] == 'X':
                    X_labels += [self.atom_types[i]]
                    X_pos += [self.atom_pos[i]]

            return X_labels, np.array(X_pos)

        else:
            print('No X ponts could be found!')
            return self.atom_types, self.atom_pos

    def get_Q_points(self, atom_types, atom_pos):
        '''Get the Q points in a molecule'''

        Q_labels, Q_pos = [], []

        for i in range(len(atom_types)):
            if atom_types[i] == 'Q':
                Q_labels += [atom_types[i]]
                Q_pos += [atom_pos[i]]

        return Q_labels, np.array(Q_pos)

    def get_R_points(self, atom_types, atom_pos):
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
            for i in range(len(atom_types)):
                if atom_types[i] == key:
                    R_dict[key] += [atom_pos[i]]

        return R_dict

    def calculate_size(self):
        '''Calculate the size of the building block'''
        _, X_pos = self.get_X_points()
        self.size = [np.linalg.norm(i) for i in X_pos]

    def align_to(self, vec: list = [0, 1, 0], n: int = 0):
        '''
        Align the first n-th X point to a given vector

        Parameters
        ----------
        vec : list
            Vector to align the molecule
        n : int
            Index of the X point to be aligned
        align_to_y : bool
            If True, the second point is aligned to the y axis
        '''
        _, X_pos = self.get_X_points()
        R_matrix = rotation_matrix_from_vectors(X_pos[n], vec)

        self.atom_pos = np.dot(self.atom_pos, np.transpose(R_matrix))

    def rotate_around(self, rotation_axis: list = [1, 0, 0], angle: float = 0.0, degree: bool = True):
        '''
        Rotate the molecule around a given axis

        Parameters
        ----------
        rotation_axis : list
            Rotation axis
        angle : float
            Rotation angle in degrees
        degree : bool
            If True, the angle is given in degrees.
            Default is True.
        '''

        rotation_axis = unit_vector(rotation_axis)
        rotation = R.from_rotvec(angle * rotation_axis, degrees=degree)

        self.atom_pos = rotation.apply(self.atom_pos)

    def rotate_to_xy_plane(self):
        '''Rotate the molecule to the xy plane'''
        _, X_pos = self.get_X_points()

        if len(X_pos) == 3:

            normal = np.cross(X_pos[0], X_pos[-1])
            if normal[0] != 0 and normal[1] != 0:
                R_matrix = rotation_matrix_from_vectors(normal, [0, 0, 1])
                self.atom_pos = np.dot(self.atom_pos, np.transpose(R_matrix))

        if len(X_pos) == 2:
            normal = np.cross(X_pos[0], self.atom_pos[1])
            if normal[0] != 0 and normal[1] != 0:
                R_matrix = rotation_matrix_from_vectors(normal, [0, 0, 1])
                self.atom_pos = np.dot(self.atom_pos, np.transpose(R_matrix))

    def shift(self, shift_vector: list):
        '''
        Shift the molecule by a given vector

        Parameters
        ----------
        shift_vector : list
            Shift vector
        '''

        self.atom_pos = np.array(self.atom_pos) + np.array(shift_vector)

    def structure_as_string(self):
        struct_string = ''
        for i, _ in enumerate(self.atom_types):
            struct_string += '{:<5s}{:>10.7f}{:>15.7f}{:>15.7f}\n'.format(self.atom_types[i],
                                                                          self.atom_pos[i][0],
                                                                          self.atom_pos[i][1],
                                                                          self.atom_pos[i][2])

        return struct_string

    def add_connection_group(self, conector_name):
        '''Adds the functional group by which the COF will be formed from the building blocks'''

        connector = ChemJSON()
        connector.from_cjson(os.path.join(self.main_path, 'conector'), conector_name)

        self.smiles = self.smiles.replace('[Q]',
                                          f"{connector.properties['smiles'].replace('[Q]', '')}")

        conector_label = connector.atomic_types
        conector_pos = connector.cartesian_positions

        # Get the position of the Q points in the structure
        location_Q_struct = self.get_Q_points(self.atom_types, self.atom_pos)

        for i in range(len(location_Q_struct[0])):

            n_conector_label = conector_label.copy()
            n_conector_pos = conector_pos.copy()

            # Get the position of the closest atom to Q in the structure
            close_Q_struct = closest_atom('Q',
                                          location_Q_struct[1][i],
                                          self.atom_types,
                                          self.atom_pos)[1]

            # Get the position of Q in the conection group
            location_Q_connector = self.get_Q_points(n_conector_label, n_conector_pos)

            # Get the position of the closest atom to Q in the conection group
            close_Q_connector = closest_atom('Q',
                                             location_Q_connector[1][0],
                                             n_conector_label,
                                             n_conector_pos)[1]

            # Create the vector Q in the structure
            v1 = close_Q_struct - location_Q_struct[1][i]
            # Create the vector Q in the conector
            v2 = np.array(close_Q_connector) - np.array(location_Q_connector[1][0])

            # Find the rotation matrix that align v2 with v1
            Rot_m = rotation_matrix_from_vectors(v2, v1)

            # Delete the "Q" atom position of the conector group and the structure
            n_conector_pos = np.delete(
                n_conector_pos,
                find_index(np.array([0., 0., 0.]), n_conector_pos),
                axis=0)

            self.atom_labels = np.delete(
                self.atom_labels,
                find_index(location_Q_struct[1][i], self.atom_pos),
                axis=0
            )

            self.atom_pos = np.delete(
                self.atom_pos,
                find_index(location_Q_struct[1][i], self.atom_pos),
                axis=0
            )

            # Rotate and translade the conector group to Q position in the strucutre
            rotated_translated_group = np.dot(n_conector_pos, -np.transpose(Rot_m)) + location_Q_struct[1][i]

            # Add the position of conector atoms to the main structure
            self.atom_pos = np.append(self.atom_pos, rotated_translated_group, axis=0)

            # Remove the Q atoms from structure
            self.atom_types.remove('Q')
            n_conector_label.remove('Q')

            self.atom_types = self.atom_types + n_conector_label

            self.atom_labels = np.append(self.atom_labels, [['Q'] * len(n_conector_label)])

    def add_connection_group_symm(self, conector_name):
        '''Adds the functional group by which the COF will be formed from the building blocks'''

        connector = ChemJSON()
        connector.from_cjson(os.path.join(self.main_path, 'conector'), conector_name)

        self.smiles = self.smiles.replace('[Q]',
                                          f"{connector.properties['smiles'].replace('[Q]', '')}")

        conector_types = connector.atomic_types
        conector_pos = connector.cartesian_positions

        # Get the position of the Q points in the structure
        _, Q_vec = self.get_Q_points(self.atom_types, self.atom_pos)

        # Remove the Q atoms from structure
        self.atom_types = self.atom_types[:-4]
        self.atom_pos = self.atom_pos[:-4]
        self.atom_labels = self.atom_labels[:-4]

        # Create the vector Q in the structure
        QS_vector = Q_vec[0]

        # Create the vector Q in the conector
        QC_vector = np.array(conector_pos[1])

        # Find the rotation matrix that align the connector with the structure
        Rot_m = rotation_matrix_from_vectors(QC_vector, QS_vector)

        # Rotate and translade the conector group to Q position in the strucutre
        conector_pos = np.dot(conector_pos, np.transpose(Rot_m)) + Q_vec[0]
        conector_pos = R.from_rotvec(
            -90 * unit_vector(conector_pos[1] - conector_pos[0]), degrees=True).apply(conector_pos)

        # Add the position of conector atoms to the main structure
        self.atom_types = list(self.atom_types) + conector_types[1:]
        self.atom_pos = list(self.atom_pos) + list(conector_pos[1:])
        self.atom_labels = np.append(self.atom_labels, [['Q'] * len(conector_types[1:])])

        # Apply the improper rotations to match S4 symmetry
        Q_vec = [unit_vector(i) for i in Q_vec]

        # First S4 axis is location_Q_struct[1]
        R1 = R.from_rotvec(120 * Q_vec[1], degrees=True).apply(conector_pos)
        R1 = R.from_rotvec(60 * unit_vector(R1[1] - R1[0]), degrees=True).apply(R1)

        # Add the position of conector atoms to the main structure
        self.atom_types = list(self.atom_types) + conector_types[1:]
        self.atom_pos = list(self.atom_pos) + list(R1[1:])
        self.atom_labels = np.append(self.atom_labels, [['Q'] * len(conector_types[1:])])

        # Second S4 axis is location_Q_struct[2]
        R2 = R.from_rotvec(120 * Q_vec[2], degrees=True).apply(conector_pos)
        R2 = R.from_rotvec(-120 * unit_vector(R2[1]-R2[0]), degrees=True).apply(R2)

        # Add the position of conector atoms to the main structure
        self.atom_types = list(self.atom_types) + conector_types[1:]
        self.atom_pos = list(self.atom_pos) + list(R2)[1:]
        self.atom_labels = np.append(self.atom_labels, [['Q'] * len(conector_types[1:])])

        # Third S4 axis is location_Q_struct[3]
        R3 = R.from_rotvec(120 * Q_vec[3], degrees=True).apply(conector_pos)
        R3 = R.from_rotvec(60 * unit_vector(R3[1] - R3[0]), degrees=True).apply(R3)

        # Add the position of conector atoms to the main structure
        self.atom_types = list(self.atom_types) + conector_types[1:]
        self.atom_pos = list(self.atom_pos) + list(R3[1:])
        self.atom_labels = np.append(self.atom_labels, [['Q'] * len(conector_types[1:])])

    def add_R_group(self, R_name, R_type):
        '''Adds group R in building blocks'''

        rgroup = ChemJSON()
        rgroup.from_cjson(os.path.join(self.main_path, 'func_groups'), R_name)

        self.smiles = self.smiles.replace(f'[{R_type}]',
                                          f"{rgroup.properties['smiles'].replace('[Q]', '')}")

        group_label = rgroup.atomic_types
        group_pos = rgroup.cartesian_positions

        # Get the position of the R points in the structure
        location_R_struct = self.get_R_points(self.atom_types, self.atom_pos)[R_type]

        # Get the position of the R points in the R group
        for i, _ in enumerate(location_R_struct):
            n_group_label = group_label.copy()
            n_group_pos = group_pos.copy()

            # Get the position of the closest atom to R in the structure
            close_R_struct = closest_atom_struc(R_type,
                                                location_R_struct[i],
                                                self.atom_types,
                                                self.atom_pos)[1]

            # Get the position of R in the R group
            pos_R_group = self.get_R_points(n_group_label, n_group_pos)['R']

            # Get the position of the closest atom to R in the R group
            close_R_group = closest_atom('R', pos_R_group[0], n_group_label, n_group_pos)[1]

            # Create the vector R in the structure
            v1 = close_R_struct - location_R_struct[i]

            # Create the vector R in the R group
            v2 = np.array(close_R_group) - np.array(pos_R_group[0])

            # Find the rotation matrix that align v2 with v1
            Rot_m = rotation_matrix_from_vectors(v2, v1)

            # Delete the "R" atom position of the R group and the structure
            n_group_pos = np.delete(
                n_group_pos,
                find_index(np.array([0.0, 0.0, 0.0]), n_group_pos),
                axis=0
            )

            # Rotate and translade the R group to R position in the strucutre
            rotated_translated_group = np.dot(
                n_group_pos,
                -np.transpose(Rot_m)
                ) + location_R_struct[i]

            # Remove the R atoms from the labels list
            # Remove the R atoms from structure
            self.atom_labels = np.delete(
                self.atom_labels,
                find_index(location_R_struct[i], self.atom_pos),
                axis=0
            )

            # Remove the R atoms from structure
            self.atom_pos = np.delete(
                self.atom_pos,
                find_index(location_R_struct[i], self.atom_pos),
                axis=0
            )

            # Add the position of rotated atoms to the main structure
            self.atom_pos = np.append(self.atom_pos, rotated_translated_group, axis=0)

            # Remove the R atoms from structure
            self.atom_types.remove(R_type)

            # Remove the R atoms from R group
            n_group_label.remove('R')

            self.atom_types = self.atom_types + n_group_label
            self.atom_labels = np.append(self.atom_labels, ['R'] * len(n_group_label))

    def create_BB_structure(self,
                            symmetry='L2',
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

        core = ChemJSON()
        core.from_cjson(os.path.join(self.main_path, 'core', symmetry), core_name)

        self.smiles = core.properties['smiles']

        self.atom_types = core.atomic_types
        self.atom_pos = core.cartesian_positions
        self.composition = core.formula
        self.atom_labels = ['C']*len(self.atom_types)

        pref_orientation = unit_vector(
            self.get_Q_points(core.atomic_types, core.cartesian_positions)[1][0])

        if symmetry == 'D4':
            self.add_connection_group_symm(conector)
        else:
            self.add_connection_group(conector)

        R_list_names = [R1, R2, R3, R4, R5, R6, R7, R8, R9]
        R_list_labels = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9']

        funcGroup_string = []
        for i in range(len(R_list_names)):
            if R_list_labels[i] in self.atom_types:
                self.add_R_group(R_list_names[i], R_list_labels[i])
                self.name += f'_{R_list_names[i]}'
                funcGroup_string.append(R_list_names[i])

        self.funcGroups = funcGroup_string

        self.connectivity = len([i for i in self.atom_types if 'X' in i])
        self.centralize()
        self.align_to(pref_orientation)
        self.calculate_size()

    def replace_X(self, target_type):
        for i in range(len(self.atom_types)):
            if self.atom_types[i] == "X":
                self.atom_types[i] = target_type

    def remove_X(self):

        atom_types, atom_pos, atom_labels = [], [], []
        for i in range(len(self.atom_types)):
            if self.atom_types[i] != "X":
                atom_types.append(self.atom_types[i])
                atom_pos.append(self.atom_pos[i])
                atom_labels.append(self.atom_labels[i])

        self.atom_types = atom_types
        self.atom_pos = np.array(atom_pos)
        self.atom_labels = atom_labels

    def save(self, extension='xyz'):

        if extension == 'xyz':
            save_xyz(path=self.bb_out_path,
                     file_name=self.name + '.xyz',
                     atom_types=self.atom_types,
                     atom_pos=self.atom_pos)

    def get_available_core(self):
        '''Get the list of available cores'''

        available_cores = {}

        for symm in self.available_symmetry:
            symm_path = os.path.join(self.main_path, 'core', symm)
            available_cores[symm] = [i.rstrip('.cjson') for i in os.listdir(symm_path) if '.cjson' in i]

        return available_cores

    def get_available_R(self):
        '''Get the list of available functional groups'''
        R_PATH = os.path.join(self.main_path, 'func_groups')
        R_list = [i.rstrip('.cjson') for i in os.listdir(R_PATH) if '.cjson' in i]

        return R_list

    def get_available_conector(self):
        '''Get the list of available conectores'''
        C_PATH = os.path.join(self.main_path, 'conector')
        C_list = [i.rstrip('.cjson') for i in os.listdir(C_PATH) if '.cjson' in i]

        return C_list

    def check_existence(self, name):

        symm_check = False
        core_check = False
        conector_check = False
        funcGroup_check = True

        name = name.split('_')
        symm = name[0]
        core = name[1]
        conector = name[2]
        funcGroups = name[3:]

        BB_dict = self.get_available_core()

        if symm in self.available_symmetry:
            symm_check = True
        else:
            print('ERROR!: Building Block symmetry must be L2, T3, S4, or H6.')
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

        possible_funcGroups_list = self.get_available_R()
        for func in funcGroups:
            if func not in possible_funcGroups_list:
                print(f'ERROR! Functional group {func} is not a available.')
                print(f'Available list: {possible_funcGroups_list}')
                funcGroup_check = False

        return symm_check, core_check, conector_check, funcGroup_check

    def get_buildingblock_list(self, shape, connector_group):

        files_list = os.listdir(self.bb_out_path)

        return [i.rstrip('.xyz') for i in files_list if shape == i.split('_')[0] and connector_group in i.split('_')[2]]
