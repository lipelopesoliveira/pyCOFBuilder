# -*- coding: utf-8 -*-
# Created by Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

"""
The CJSON package implements functions to read, create and manipulate Chemical JSON objects.
"""

import os
import simplejson
import numpy as np

from pycofbuilder.tools import elements_dict

import gemmi
from ase.cell import Cell


class ChemJSON:
    '''
    Class to read, create and manupulate ChemJSON files.

    Attributes
    ----------

    file_name : str
        The name of the file.
    name : str
        The name of the structure.
    cell_parameters : list
        The cell parameters of the structure as a (1,6) list.
    cell_matrix : list
        The cell matrix of the structure as a (3,3) list.
    cartesian_positions : list
        The cartesian positions of the structure as a (n,3) list.
    fractional_positions : list
        The fractional positions of the structure as a (n,3) list.
    atomic_numbers : list
        The atomic numbers of the structure as a (n,1) list.
    atomic_types : list
        The atomic types of the structure as a (n,1) list.
    atomic_labels : list
        The atomic labels of the structure as a (n,1) list.
    formula : str
        The formula of the structure.
    properties : dict
        The properties of the structure.
    partial_charges : dict
        A dictionary contaning the partial charges of the atoms on the structure.
        Example: {'DDEC': [0.1, 0.2, 0.15], 'EQeq': [0.05, 0.15, 0.19]}
    '''
    def __init__(self):
        self.file_name = ''
        self.name = ''

        # Structure properties
        self.cell_parameters = None
        self.cell_matrix = None
        self.cartesian_positions = None
        self.fractional_positions = None
        self.atomic_numbers = None
        self.atomic_types = None
        self.atomic_labels = None
        self.formula = ''
        self.partial_charges = None

        self.properties = None
        self.results = []

    # Create a custom representation of the class
    def __repr__(self):
        '''
        Returns a custom representation of the class.
        '''

        repr_string = "ChemJSON(name='{}', formula='{}', number of atoms={}".format(self.name,
                                                                                    self.formula,
                                                                                    len(self.atomic_types))

        return repr_string

    # Create a custom print of the class
    def __str__(self):
        '''
        Returns a custom print of the class.
        '''

        string_string = "ChemJSON(name='{}', formula='{}', number of atoms={})\n".format(self.name,
                                                                                         self.formula,
                                                                                         len(self.atomic_types))

        if self.cell_parameters is not None:
            string_string += f"""Cell parameters:
    a = {self.cell_parameters[0]:>12.7f} Å
    b = {self.cell_parameters[1]:>12.7f} Å
    c = {self.cell_parameters[2]:>12.7f} Å
    α = {self.cell_parameters[3]:>12.7f} °
    β = {self.cell_parameters[4]:>12.7f} °
    γ = {self.cell_parameters[5]:>12.7f} °

Cell matrix:
A   {self.cell_matrix[0][0]:>12.7f}  {self.cell_matrix[0][1]:>12.7f} {self.cell_matrix[0][2]:>12.7f}
B   {self.cell_matrix[1][0]:>12.7f}  {self.cell_matrix[1][1]:>12.7f} {self.cell_matrix[1][2]:>12.7f}
C   {self.cell_matrix[2][0]:>12.7f}  {self.cell_matrix[2][1]:>12.7f} {self.cell_matrix[2][2]:>12.7f}
"""
        if self.cartesian_positions is not None:
            string_string += "Cartesian positions:\n"
            for i, position in enumerate(self.cartesian_positions):
                string_string += "    {:3} {:>9.5f}   {:>9.5f}   {:>9.5f}\n".format(self.atomic_types[i],
                                                                                    position[0],
                                                                                    position[1],
                                                                                    position[2]
                                                                                    )

        if self.fractional_positions is not None:
            string_string += "Fractional positions:\n"
            for i, position in enumerate(self.fractional_positions):
                string_string += "    {:3} {:>9.5f}   {:>9.5f}   {:>9.5f}\n".format(self.atomic_types[i],
                                                                                    position[0],
                                                                                    position[1],
                                                                                    position[2])

        return string_string

    def set_properties(self, properties):
        '''
        Sets the properties of the structure.
        '''
        self.properties = properties

    def set_results(self, results):
        '''
        Sets the results of the structure.
        '''
        self.results = results

    def set_cell_parameters(self, cell_parameters):
        '''
        Sets the cell parameters of the structure.
        '''
        self.cell_parameters = cell_parameters

        aseCell = Cell.fromcellpar(cell_parameters)
        self.cell_matrix = np.array(aseCell)

    def set_cell_matrix(self, cell_matrix):
        '''
        Sets the cell matrix of the structure. The cell
        parameters will be calculated and also updated.
        '''
        self.cell_matrix = cell_matrix

        aseCell = Cell(cell_matrix)
        self.cell_parameters = aseCell.cellpar()

    def set_cartesian_positions(self, cartesian_positions):
        '''
        Sets the cartesian positions of the structure. The fractional
        positions will be calculated and also updated.
        '''
        self.cartesian_positions = np.array(cartesian_positions).astype(float)

        if self.cell_parameters is not None:
            aseCell = Cell.fromcellpar(self.cell_parameters)
            self.fractional_positions = aseCell.scaled_positions(cartesian_positions)

    def set_fractional_positions(self, fractional_positions):
        '''
        Sets the fractional positions of the structure. The cartesian
        positions will be calculated and also updated.
        '''
        self.fractional_positions = np.array(fractional_positions).astype(float)

        if self.cell_parameters is not None:
            aseCell = Cell.fromcellpar(self.cell_parameters)
            self.cartesian_positions = aseCell.cartesian_positions(fractional_positions)

    def set_atomic_types(self, atomic_types):
        '''
        Sets the atomic labels of the structure.
        '''
        self.atomic_types = atomic_types

        self.atomic_labels = [f"{atom}{i+1}" for i, atom in enumerate(self.atomic_types)]

        symbol_dict = elements_dict('atomic_number')
        self.atomic_numbers = [symbol_dict[i] for i in atomic_types]

        self.formula = ''.join([f'{atom}{self.atomic_types.count(atom)}' for atom in set(self.atomic_types)])

    def set_atomic_numbers(self, atomic_numbers):
        '''
        Sets the atomic numbers of the structure. The atomic types and formula
        will be calculated and also updated.
        '''
        self.atomic_numbers = atomic_numbers

        symbol_dict = elements_dict('atomic_number')
        number_dict = {j: i for i, j in zip(symbol_dict.keys(), symbol_dict.values())}

        self.atomic_types = [number_dict[i] for i in atomic_numbers]

        self.atomic_labels = [f"{atom}{i+1}" for i, atom in enumerate(self.atomic_types)]

        self.formula = ''.join([f'{atom}{self.atomic_types.count(atom)}' for atom in set(self.atomic_types)])

    def from_cjson(self, path, file_name):
        '''
        Reads a ChemJSON file from a given path and file_name.
        '''
        self.file_name = os.path.join(path, file_name.split('.')[0] + '.cjson')

        with open(self.file_name, 'r') as file:
            cjson_data = simplejson.load(file)

        if "name" in cjson_data:
            self.name = cjson_data['name']

        if "unitCell" in cjson_data:
            if 'a' in cjson_data['unitCell']:
                self.set_cell_parameters(
                    [cjson_data['unitCell'][i] for i in ['a', 'b', 'c', 'alpha', 'beta', 'gamma']]
                )
            elif 'cellVectors' in cjson_data['unitCell']:
                self.set_cell_matrix(
                    np.array(cjson_data['unitCell']['cellVectors']).reshape(3, 3)
                )

        if "atoms" in cjson_data:
            if 'coords' in cjson_data['atoms']:
                if '3d' in cjson_data['atoms']['coords']:
                    self.set_cartesian_positions(
                        np.array(cjson_data['atoms']['coords']['3d']).reshape(-1, 3)
                    )
                elif '3dFractional' in cjson_data['atoms']['coords']:
                    self.set_fractional_positions(
                        np.array(cjson_data['atoms']['coords']['3dFractional']).reshape(-1, 3)
                    )

            if "elements" in cjson_data['atoms']:
                if 'type' in cjson_data['atoms']['elements']:
                    self.set_atomic_types(cjson_data['atoms']['elements']['type'])

                elif 'number' in cjson_data['atoms']['elements']:
                    self.set_atomic_numbers(cjson_data['atoms']['elements']['number'])

                if 'label' in cjson_data['atoms']['elements']:
                    self.atomic_labels = cjson_data['atoms']['elements']['label']
                else:
                    self.atomic_labels = [f"{atom}{i+1}" for i, atom in enumerate(self.atomic_types)]

        if 'properties' in cjson_data:
            self.set_properties(cjson_data['properties'])
        if 'results' in cjson_data:
            self.set_results(cjson_data['results'])
        if 'partialCharges' in cjson_data:
            self.partial_charges = cjson_data['partialCharges']

    def from_xyz(self, path, file_name):
        '''
        Reads a XYZ file from a given path and file_name.
        '''

        self.file_name = os.path.join(path, file_name.split('.')[0] + '.xyz')
        self.name = file_name.split('.')[0]

        with open(self.file_name, 'r') as file:
            xyz_data = file.read().splitlines()
        n_atoms = int(xyz_data[0])

        atomic_types = []
        cartesian_positions = []

        for line in xyz_data[2: n_atoms + 3]:
            atomic_types.append(line.split()[0])
            cartesian_positions.append([float(i) for i in line.split()[1:]])

        self.set_atomic_types(atomic_types)

        self.set_cartesian_positions(np.array(cartesian_positions))

    def from_gjf(self, path, file_name):
        '''
        Reads a Gaussian input file from a given path and file_name.
        '''

        self.file_name = os.path.join(path, file_name.split('.')[0] + '.gjf')
        self.name = file_name.split('.')[0]

        with open(self.file_name, 'r') as file:
            gjf_data = file.read().splitlines()

        # Remove empty lines
        gjf_data = [line for line in gjf_data if line != '']

        atomic_types = []
        cartesian_positions = []

        for line in gjf_data:
            if line.split()[0] in elements_dict('atomic_number').keys():
                atomic_types.append(line.split()[0])
                cartesian_positions.append([float(i) for i in line.split()[1:4]])

        cell_matrix = []
        for line in gjf_data:
            if line.split()[0] == 'Tv':
                cell_matrix.append([float(i) for i in line.split()[1:4]])

        if cell_matrix != []:
            self.set_cell_matrix(np.array(cell_matrix))

        self.set_atomic_types(atomic_types)

        self.set_cartesian_positions(np.array(cartesian_positions))

    def from_cif(self, path, file_name):
        '''
        Reads a CIF file from a given path and file_name.
        '''

        # Read the cif file and get the lattice parameters and atomic positions
        cif_filename = os.path.join(path, file_name.split('.')[0] + '.cif')

        cif = gemmi.cif.read_file(cif_filename).sole_block()

        a = float(cif.find_value('_cell_length_a').split('(')[0])
        b = float(cif.find_value('_cell_length_b').split('(')[0])
        c = float(cif.find_value('_cell_length_c').split('(')[0])
        beta = float(cif.find_value('_cell_angle_beta').split('(')[0])
        gamma = float(cif.find_value('_cell_angle_gamma').split('(')[0])
        alpha = float(cif.find_value('_cell_angle_alpha').split('(')[0])

        CellParameters = [a, b, c, alpha, beta, gamma]

        AtomicTypes = list(cif.find_values('_atom_site_type_symbol'))
        PosX = np.array(cif.find_values('_atom_site_fract_x')).astype(float)
        PosY = np.array(cif.find_values('_atom_site_fract_y')).astype(float)
        PosZ = np.array(cif.find_values('_atom_site_fract_z')).astype(float)
        try:
            charges = np.array(cif.find_values('_atom_site_charge')).astype(float)
            charge_type = 'DDEC'
        except Exception:
            charges = None
            charge_type = None

        self.set_cell_parameters(CellParameters)

        self.set_atomic_types(AtomicTypes)

        self.set_fractional_positions(np.array([PosX, PosY, PosZ]).T)

        if charges is not None:
            self.partial_charges = {charge_type: charges}

    def as_dict(self):
        '''
        Returns the structure as a dictionary.
        '''
        structure_dict = {
            'chemical json': 1,
            'name': self.name,
            'formula': self.formula,
        }
        if self.cell_parameters is not None:
            structure_dict['unit cell'] = {
                'a': self.cell_parameters[0],
                'b': self.cell_parameters[1],
                'c': self.cell_parameters[2],
                'alpha': self.cell_parameters[3],
                'beta': self.cell_parameters[4],
                'gamma': self.cell_parameters[5],
                'cellVectors': self.cell_matrix.flatten().tolist()
            }

        structure_dict['atoms'] = {
                'elements': {
                    'type': self.atomic_types,
                    'number': self.atomic_numbers,
                },
                'coords': {
                    '3d': self.cartesian_positions.flatten().tolist(),
                }
            }

        if self.cell_parameters is not None:
            structure_dict['atoms']['coords']['3dFractional'] = self.fractional_positions.flatten().tolist()

        if self.partial_charges is not None:
            structure_dict['partialCharges'] = self.partial_charges

        structure_dict['properties'] = self.properties

        structure_dict['results'] = self.results

        return structure_dict

    def write_cjson(self, path, file_name):
        '''
        Writes a ChemJSON file to a given path and file_name.
        '''
        self.file_name = os.path.join(path, file_name.split('.')[0] + '.cjson')
        with open(self.file_name, 'w') as file:
            simplejson.dump(self.as_dict(), file, indent=4)
