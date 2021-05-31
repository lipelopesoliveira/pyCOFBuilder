# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020

@author: lipel
"""

import os
import numpy as np
import pycofbuilder.tools as Tools

from scipy.spatial import distance

class Building_Block():

    def __init__(self, lib='bb_lib', name=None, verbosity=False):

        self.name = name
        self.verbosity = verbosity
        self.main_path = 'data'
        self.lib_path = os.path.join(self.main_path, lib)
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

    def n_atoms(self):
        return len(self.atom_labels)

    def centralize_molecule(self, by_X=False):
        ''' Centralize the molecule in its geometrical center'''

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

    
        Q_labels, Q_pos = [], []

        for i in range(len(atom_labels)):
            if atom_labels[i] == 'Q':
                Q_labels += [atom_labels[i]]
                Q_pos += [atom_pos[i]]

        return Q_labels, np.array(Q_pos)

    def get_R_points(self, atom_labels, atom_pos):

        R_pos = []
        R1_pos = []
        R2_pos = []
        R3_pos = []
        R4_pos = []
        R5_pos = []
        R6_pos = []

        if 'R' in atom_labels:
            for i in range(len(atom_labels)):
                if atom_labels[i] == 'R':
                    R_pos += [atom_pos[i]]

        if 'R1' in atom_labels:
            for i in range(len(atom_labels)):
                if atom_labels[i] == 'R1':
                    R1_pos += [atom_pos[i]]

        if 'R2' in atom_labels:
            for i in range(len(atom_labels)):
                if atom_labels[i] == 'R2':
                    R2_pos += [atom_pos[i]]

        if 'R3' in atom_labels:
            for i in range(len(atom_labels)):
                if atom_labels[i] == 'R3':
                    R3_pos += [atom_pos[i]]

        if 'R4' in atom_labels:
            for i in range(len(atom_labels)):
                if atom_labels[i] == 'R4':
                    R4_pos += [atom_pos[i]]

        if 'R5' in atom_labels:
            for i in range(len(atom_labels)):
                if atom_labels[i] == 'R5':
                    R5_pos += [atom_pos[i]]

        if 'R6' in atom_labels:
            for i in range(len(atom_labels)):
                if atom_labels[i] == 'R6':
                    R6_pos += [atom_pos[i]]

        return {'R': R_pos, 'R1': R1_pos, 'R2': R2_pos, 'R3': R3_pos, 'R4': R4_pos, 'R5': R5_pos, 'R6': R6_pos}

    def add_R(self, label, pos, R1='H', R2='H', R3='H', R4='H', R5='H', R6='H'):

        for i in range(len(label)):
            if label[i] == 'R1':
                label[i] = R1

            if label[i] == 'R2':
                label[i] = R2

            if label[i] == 'R3':
                label[i] = R3

            if label[i] == 'R4':
                label[i] = R4

            if label[i] == 'R5':
                label[i] = R5

            if label[i] == 'R6':
                label[i] = R6

        return label, pos

    def add_X(self, label, pos, X='N'):

        label, pos = label, pos

        for i in range(len(label)):
            if label[i] == 'X':
                label[i] = X
        return label, pos

    def find_index(self, element, e_list):
        for i in range(len(e_list)):
            if np.array_equal(e_list[i], element):
                return i

    def closest_atom(self, label_1, pos_1, labels, pos):

        list_labels = []
        list_pos = []

        for i in range(len(labels)):
            if labels[i] != label_1:
                list_labels += [labels[i]]
                list_pos += [pos[i]]

        closest_index = distance.cdist([pos_1], list_pos).argmin()

        return list_labels[closest_index], list_pos[closest_index], np.linalg.norm(pos_1 - list_pos[closest_index])

    def closest_atom_struc(self, label_1, pos_1, labels, pos):
        list_labels = []
        list_pos = []
        for i in range(len(labels)):
            if labels[i] != label_1:
                if 'C' in labels[i]:
                    list_labels += [labels[i]]
                    list_pos += [pos[i]]
        closest_index = distance.cdist([pos_1], list_pos).argmin()

        return list_labels[closest_index], list_pos[closest_index], np.linalg.norm(pos_1-list_pos[closest_index])

    def calculate_size(self):
        X_labels, X_pos = self.get_X_points()
        self.size = [np.linalg.norm(i) for i in X_pos]

    def align_to(self, vec=[0, 1, 0]):

        X_labels, X_pos = self.get_X_points()
        R_matrix = Tools.rotation_matrix_from_vectors(X_pos[0], vec)

        self.atom_pos = np.dot(self.atom_pos, np.transpose(R_matrix))

    def rotate_to_xy_plane(self):

        X_labels, X_pos = self.get_X_points()

        if len(X_pos) == 3:

            normal = np.cross(X_pos[0], X_pos[-1])
            if normal[0] != 0 and normal[1] != 0:
                R_matrix = self.rotation_matrix_from_vectors(normal, [0, 0, 1])
                self.atom_pos = np.dot(self.atom_pos, np.transpose(R_matrix))

        if len(X_pos) == 2:
            normal = np.cross(X_pos[0], self.atom_pos[1])
            if normal[0] != 0 and normal[1] != 0:
                R_matrix = self.rotation_matrix_from_vectors(normal, [0, 0, 1])
                self.atom_pos = np.dot(self.atom_pos, np.transpose(R_matrix))

    def print_structure(self):

        for i in range(len(self.atom_labels)):
            print('{:<5s}{:>10.7f}{:>15.7f}{:>15.7f}'.format(self.atom_labels[i], self.atom_pos[i][0], self.atom_pos[i][1], self.atom_pos[i][2]))

    def add_connection_group(self, conector_name):
        '''Adiciona o grupo funcional pelo qual o COF será formado em blocos de construção'''

        conector_label, conector_pos = Tools.read_gjf_file(os.path.join(self.main_path, 'conector'), conector_name)

        location_Q_struct = self.get_Q_points(self.atom_labels, self.atom_pos)  # Get the position of the Q points in the structure

        for i in range(len(location_Q_struct[0])):
            n_conector_label = conector_label.copy()
            n_conector_pos = conector_pos.copy()

            # Pega a posição do átomo mais próximo ao primeiro Q na estrutrua
            close_Q_struct = self.closest_atom('Q', location_Q_struct[1][i], self.atom_labels, self.atom_pos)[1]

            location_Q_connector = self.get_Q_points(n_conector_label, n_conector_pos)  # Pega a posição de Q no grupo conector

            # Pega a posição do átomo mais próximo a Q no grupo conector
            close_Q_connector = self.closest_atom('Q', location_Q_connector[1][0], n_conector_label, n_conector_pos)[1]

            v1 = close_Q_struct - location_Q_struct[1][i]  # Cria o vetor Q na estrutura
            v2 = np.array(close_Q_connector) - np.array(location_Q_connector[1][0])  # Cria o vetor Q no conector

            Rot_m = Tools.rotation_matrix_from_vectors(v2, v1)  # Determina a matriz de rotação que alinha V2 com V1

            # Deleta o átomo "Q" da lista de átomos e posições do conector
            n_conector_pos = np.delete(n_conector_pos, self.find_index(np.array([0., 0., 0.]), n_conector_pos), axis=0)

            # Rotaciona e translada o grupo radical para a posição Q
            rotated_translated_group = np.dot(n_conector_pos, -np.transpose(Rot_m)) + location_Q_struct[1][i]

            self.atom_pos = np.delete(self.atom_pos, self.find_index(location_Q_struct[1][i], self.atom_pos), axis=0)

            self.atom_pos = np.append(self.atom_pos, rotated_translated_group, axis=0)

            self.atom_labels.remove('Q')

            n_conector_label.remove('Q')

            self.atom_labels = self.atom_labels + n_conector_label

    def add_R_group(self, R_name, R_type):
        '''Adiciona o grupo em blocos de construção com simetria C4. '''

        group_label, group_pos = Tools.read_gjf_file(os.path.join(self.main_path, 'radical'), R_name)

        location_R_struct = self.get_R_points(self.atom_labels, self.atom_pos)[R_type]  # Pega a posição dos R na estrutura

        for i in range(len(location_R_struct)):
            n_group_label = group_label.copy()
            n_group_pos = group_pos.copy()

            # Pega a posição do átomo mais próximo ao R na estrutrua
            close_R_struct = self.closest_atom_struc(R_type, location_R_struct[i], self.atom_labels, self.atom_pos)[1]

            pos_R_group = self.get_R_points(n_group_label, n_group_pos)['R']  # Pega a posição de Q no grupo radical
            close_R_group = self.closest_atom('R', pos_R_group[0], n_group_label, n_group_pos)[1]  # Pega a posição do átomo mais próximo a Q no grupo radical

            v1 = close_R_struct - location_R_struct[i]  # Cria o vetor R na estrutura
            v2 = np.array(close_R_group) - np.array(pos_R_group[0])  # Cria o vetor R do grupo

            Rot_m = Tools.rotation_matrix_from_vectors(v2, v1)  # Determina a matriz de rotação que alinha V2 com V1

            n_group_pos = np.delete(n_group_pos, self.find_index(np.array([0.0, 0.0, 0.0]), n_group_pos), axis=0)

            # Rotaciona e translada o grupo radical para a posição R
            rotated_translated_group = np.dot(n_group_pos, - np.transpose(Rot_m)) + location_R_struct[i]

            self.atom_pos = np.delete(self.atom_pos, self.find_index(location_R_struct[i], self.atom_pos), axis=0)

            self.atom_pos = np.append(self.atom_pos, rotated_translated_group, axis=0)

            self.atom_labels.remove(R_type)

            n_group_label.remove('R')

            self.atom_labels = self.atom_labels + n_group_label

    def create_C2_BB(self, nucleo_name='BENZ', conector='CHO', R1='H', R2='H', R3='H', R4='H', R5='H', R6='H'):
        '''Create a building block with C2 simmetry'''

        self.name = f'C2_{nucleo_name}_{conector}'

        self.atom_labels, self.atom_pos = Tools.read_gjf_file(os.path.join(self.main_path, 'nucleo', 'C2'), nucleo_name)

        self.centralize_molecule()

        self.add_connection_group(conector)

        R_list_names = [R1, R2, R3, R4, R5, R6]
        R_list_labels = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']

        for i in range(len(R_list_names)):
            if R_list_labels[i] in self.atom_labels:
                self.add_R_group(R_list_names[i], R_list_labels[i])
                self.name += f'_{R_list_names[i]}'

        self.align_to()

    def create_C3_BB(self, nucleo_name='BENZ', conector='CHO', R1='H', R2='H', R3='H', R4='H', R5='H', R6='H'):
        '''Cria o bloco de construção com simetria C3'''

        self.name = f'C3_{nucleo_name}_{conector}'

        self.atom_labels, self.atom_pos = Tools.read_gjf_file(os.path.join(self.main_path, 'Nucleo', 'C3'), nucleo_name)

        self.centralize_molecule()

        self.add_connection_group(conector)

        R_list_names = [R1, R2, R3, R4, R5, R6]
        R_list_labels = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']

        for i in range(len(R_list_names)):
            if R_list_labels[i] in self.atom_labels:
                self.add_R_group(R_list_names[i], R_list_labels[i])
                self.name += f'_{R_list_names[i]}'

        self.align_to()

    def create_C4_BB(self, nucleo_name='BENZ', conector='CHO', R1='H', R2='H', R3='H', R4='H', R5='H', R6='H'):
        '''Cria o bloco de construção com simetria C3'''

        self.name = f'C4_{nucleo_name}_{conector}'

        self.atom_labels, self.atom_pos = Tools.read_gjf_file(os.path.join(self.main_path, 'Nucleo', 'C4'), nucleo_name)

        self.centralize_molecule()

        self.add_connection_group(conector)

        R_list_names = [R1, R2, R3, R4, R5, R6]
        R_list_labels = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']

        for i in range(len(R_list_names)):
            if R_list_labels[i] in self.atom_labels:
                self.add_R_group(R_list_names[i], R_list_labels[i])
                self.name += f'_{R_list_names[i]}'

        self.align_to()

    def save(self, extension='xyz'):

        if extension == 'xyz':
            Tools.save_xyz(self.lib_path, self.name + '.xyz', self.atom_labels, self.atom_pos)

    def read_structure(self):

        try:
            self.atom_labels, self.atom_pos, self.n_atoms, self.connectivity = Tools.read_xyz_file(self.lib_path, self.name)
        except Exception:
            print(f'Unable to load file {self.lib_path}\\{self.name}.xyz!')

        # self.centralize_molecule(True)
        # self.rotate_to_xy_plane()
        self.align_to()
        self.calculate_size()

    def get_bipodal_NH2(self):

        files_list = os.listdir(self.lib_path)

        return [i.rstrip('.xyz') for i in files_list if 'C2' == i.split('_')[0] and 'NH2' in i.split('_')[2]]

    def get_tripodal_NH2(self):

        files_list = os.listdir(self.lib_path)

        return [i.rstrip('.xyz') for i in files_list if 'C3' == i.split('_')[0] and 'NH2' in i.split('_')[2]]

    def get_bipodal_CHO(self):

        files_list = os.listdir(self.lib_path)

        return [i.rstrip('.xyz') for i in files_list if 'C2' == i.split('_')[0] and 'CHO' == i.split('_')[2]]

    def get_tripodal_CHO(self):

        files_list = os.listdir(self.lib_path)

        return [i.rstrip('.xyz') for i in files_list if 'C3' == i.split('_')[0] and 'CHO' == i.split('_')[2]]

    def get_bipodal_BOH2(self):

        files_list = os.listdir(self.lib_path)

        return [i.rstrip('.xyz') for i in files_list if 'C2' == i.split('_')[0] and 'BOH2' == i.split('_')[2]]

    def get_tripodal_BOH2(self):

        files_list = os.listdir(self.lib_path)

        return [i.rstrip('.xyz') for i in files_list if 'C3' == i.split('_')[0] and 'BOH2' == i.split('_')[2]]

    def get_bipodal_OH2(self):

        files_list = os.listdir(self.lib_path)

        return [i.rstrip('.xyz') for i in files_list if 'C2' == i.split('_')[0] and 'OH2' == i.split('_')[2]]

    def get_tripodal_OH2(self):

        files_list = os.listdir(self.lib_path)

        return [i.rstrip('.xyz') for i in files_list if 'C3' == i.split('_')[0] and 'OH2' == i.split('_')[2]]

    def get_tetrapodal_squared_CHO(self):

        files_list = os.listdir(self.lib_path)

        return [i.rstrip('.xyz') for i in files_list if 'C2p' == i.split('_')[0] and 'CHO' == i.split('_')[2]]

    def get_tetrapodal_squared_BOH2(self):

        files_list = os.listdir(self.lib_path)

        return [i.rstrip('.xyz') for i in files_list if 'C2p' == i.split('_')[0] and 'BOH2' == i.split('_')[2]]

    def get_tetrapodal_squared_NH2(self):

        files_list = os.listdir(self.lib_path)

        return [i.rstrip('.xyz') for i in files_list if 'C2p' == i.split('_')[0] and 'NH2' in i.split('_')[2]]