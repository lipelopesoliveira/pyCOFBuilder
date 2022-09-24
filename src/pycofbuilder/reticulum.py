# -*- coding: utf-8 -*-
# Copyright (c) Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

"""
This class implements definitions for Reticulum buiding
"""

import os
import numpy as np

from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from scipy.spatial.transform import Rotation as R

import pycofbuilder.tools as Tools
from pycofbuilder.building_block import Building_Block

class Reticulum():
    """
    A class used to represent a Reticulum

    ...

    Attributes
    ----------
    verbosity : bool
        control the printing options
    available_2D_topologies : list
        List of available 2D topologies
    available_3D_topologies : list
        List of available 3D topologies
    available_topologies : list
        List of all available topologies
    available_stacking : list
        List of available stakings for all 2D topologies
    lib_bb : str
        String with the name of the folder containing the building block files
        Default: bb_lib
    main_path : str
        String containing the data folder.
        Defailt: os.path.join(_ROOT, 'data')
    lib_path : str
        Path for the building block files.
        Default: os.path.join(self.main_path, bb_lib)
    out_path : str
        Path to save the results.
        Default: os.path.join(os.getcwd(), 'out')
    name : str
        Name of the material
    topology : str = None
    dimention : str = None
    lattice : str = None
    lattice_sgs : str = None
    space_group : str = None
    space_group_n : str = None
    stacking : str = None
    mass : str = None
    composition : str = None
    charge : int  = 0
    multiplicity : int = 1
    chirality : bool = False
    atom_labels : list = []
    atom_pos : list = []
    lattice : list = [[], [], []]
    symm_tol : float = 0.2
    angle_tol : float = 0.2

    Methods
    -------
    n_atoms()
        Returns the number of atoms in the unitary cell
    print_available_topologies()
        Print all available topologies
    create_hcb_structure()
        Creates a COF with HCB network.
    create_hcb_a_structure()
        Creates a COF with HCB-A network.
    create_sql_structure()
        Creates a COF with SQL network.
    create_sql_a_structure()
        Creates a COF with SQL-A network.
    create_kgd_structure()
        Creates a COF with KGD network.
    create_hxl_a_structure()
        Creates a COF with HXL-A network.
    create_kgm_structure()
        Creates a COF with KGM network.
    create_kgm_a_structure()
        Creates a COF with KGM-A network.
    save_cif()
        Save the structure in .cif format
    save_json()
        Save the structure in .json format
    save_xsf()
        Save the structure in .xsf format
    save_pdb()
        Save the structure in .pdb format
    save_vasp()
        Save the structure in .vasp format
    save_turbomole()
        Save the structure in .coord format
    save_xyz()
        Save the structure in .xyz format
    save_qe()
        Save the structure in .in format
    """

    def __init__(self, verbosity=False, out_dir=None):

        _ROOTDIR = os.path.abspath(os.path.dirname(__file__))

        self.verbosity = verbosity
        # Falta adicionar: 'HXL', 'KGD'
        self.available_2D_topologies = ['HCB', 'HCB_A',
                                        'SQL', 'SQL_A',
                                        'KGM', 'KGM_A',
                                        'KGD',
                                        'HXL_A']

        # Falta adicionar: ['dia', 'bor', 'srs', 'pts', 'ctn', 'rra', 'fcc',
        # 'lon', 'stp', 'acs', 'tbo', 'bcu', 'fjh', 'ceq']
        self.available_3D_topologies = ['DIA', 'BOR'] # Temporary
        self.available_topologies = self.available_2D_topologies + self.available_3D_topologies

        # Define available stackings for all 2D topologies
        self.available_stacking = {'HCB': ['AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'HCB_A': ['AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'SQL':['AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'SQL_A':['AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'KGD': ['A', 'AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'HXL_A': ['A', 'AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'KGM': ['AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'KGM_A': ['AA', 'AB1x', 'AB1y', 'AB1xy', 'AB2', 'AAl', 'AAt'],
                                   'DIA': [1, 2, 3, 4],  # Temporary
                                   'BOR': [5, 8, 6, 7]  # Temporary
                                   }

        self.main_path = os.path.join(_ROOTDIR, 'data')

        if out_dir is None:
            self.out_path = os.path.join(os.getcwd(), 'out')
        else:
            self.out_path = out_dir

        self.lib_path = os.path.join(self.out_path, 'building_blocks')

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
        self.symm_tol = 0.1
        self.angle_tol = 0.1
        self.n_atoms = self.get_n_atoms()

    def get_n_atoms(self):
        ''' Returns the number of atoms in the unitary cell'''
        return len(self.atom_labels)

    def get_available_topologies(self, dimensionality='all'):

        if dimensionality == 'all' or dimensionality == '2D':
            print('Available 2D Topologies:')
            for i in self.available_2D_topologies:
                print(i.upper())

        if dimensionality == 'all' or dimensionality == '3D':
            print('Available 3D Topologies:')
            for i in self.available_3D_topologies:
                print(i.upper())

################   Net creation methods  ############################

    def create_hcb_structure(self,
                             name_bb_a: str,
                             name_bb_b: str,
                             stacking: str = 'AA',
                             bond_atom: str = 'N',
                             c_parameter_base: float = 3.6,
                             print_result: bool = True,
                             slab: float = 10.0,
                             shift_vector: list = [1.0 ,1.0 ,0],
                             tilt_angle: float = 5.0):
        """Creates a COF with HCB network.

        The HCB net is composed of two tripodal building blocks.

        Parameters
        ----------
        name_bb_a : str, required
            The 4 letter code for the tripodal Buiding Block A
        name_bb_b : str, required
            The 4 letter code for the tripodal Buiding Block B
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
        bond_atom : str, optional
            The atom that connect the buiding blocks (default is 'N')
        c_parameter_base : float, optional
            The base value for interlayer distance in angstroms (default is 3.6)
        print_result : bool, optional
            Parameter for the control for printing the result (default is True)
        slab : float, optional
            Default parameter for the interlayer slab (default is 10.0)
        shift_vector: list, optional
            Shift vector for the AAl and AAt stakings (defatult is [1.0,1.0,0])
        tilt_angle: float, optional
            Tilt angle for the AAt staking in degrees (default is 5.0)

        Returns
        -------
        list
            A list of strings containing:
                1. the structure name,
                2. lattice type,
                3. hall symbol of the cristaline structure,
                4. space group,
                5. number of the space group,
                6. number of operation symmetry
        """

        self.topology = 'HCB'
        self.dimension = 2

        bb_1 = Building_Block(name_bb_a,
                              save_dir=self.lib_path,
                              verbosity=self.verbosity)

        bb_2 = Building_Block(name_bb_b,
                              save_dir=self.lib_path,
                              verbosity=self.verbosity)

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
        r_pos_a_2 = np.dot(
            bb_2.atom_pos, 
            R.from_euler('z', 180, degrees=True).as_matrix()
            ) + np.array([0, np.sqrt(3)/3, 0])*a

        final_pos = np.vstack((final_pos, r_pos_a_2))
        final_label += bb_2.atom_labels

        # Changes the X atoms by the desirede bond_atom
        final_label, final_pos = Tools.change_X_atoms(final_label, final_pos, bond_atom)

        # Creates the COF as a pymatgen structure
        struct = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)

        # Remove duplicate atoms and translate the structure to the center of the cell
        struct.sort(reverse=True)
        struct.merge_sites(tol=.5, mode='delete')
        struct.translate_sites(range(len(struct.as_dict()['sites'])),
                               [0, 0, 0.5],
                               frac_coords=True,
                               to_unit_cell=True)

        # Simetrizes the structure using pymatgen
        symm = SpacegroupAnalyzer(struct, symprec=self.symm_tol, angle_tolerance=self.angle_tol)
        struct_symm_prim = symm.get_primitive_standard_structure()

        if stacking not in self.available_stacking[self.topology]:
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
            AB_1_symm = SpacegroupAnalyzer(AB_1, 
                                           symprec=self.symm_tol, 
                                           angle_tolerance=self.angle_tol)

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
            
        # Create AAl stacking. 
        # Hexagonal cell with two sheets per cell shifited by the shift_vector in angstroms. 
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

        # Create AA tilted stacking. 
        # Tilted Hexagonal cell by tilt_angle with two sheets per cell 
        # shifited by the shift_vector in angstroms.
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

        if print_result is True:
            Tools.print_result(self.name,
                               str(self.lattice_type),
                               str(self.hall[0:2]),
                               str(self.space_group),
                               str(self.space_group_n),
                               len(symm_op))

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]


    def create_hcb_a_structure(self,
                               name_bb_a: str,
                               name_bb_b: str,
                               stacking: str = 'AA',
                               bond_atom: str = 'N',
                               c_parameter_base: float = 3.6,
                               print_result: bool = True,
                               slab: float = 10.0,
                               shift_vector: list = [1.0, 1.0, 0],
                               tilt_angle: float = 5.0):
        """Creates a COF with HCB-A network.

        The HCB-A net is composed of one tripodal and one linear building blocks.

        Parameters
        ----------
        name_bb_a : str, required
            The 4 letter code for the tripodal Buiding Block A
        name_bb_b : str, required
            The 4 letter code for the linear Buiding Block B
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
        bond_atom : str, optional
            The atom that connect the buiding blocks (default is 'N')
        c_parameter_base : float, optional
            The base value for interlayer distance in angstroms (default is 3.6)
        print_result : bool, optional
            Parameter for the control for printing the result (default is True)
        slab : float, optional
            Default parameter for the interlayer slab (default is 10.0)
        shift_vector: list, optional
            Shift vector for the AAl and AAt stakings (defatult is [1.0,1.0,0])
        tilt_angle: float, optional
            Tilt angle for the AAt staking in degrees (default is 5.0)

        Returns
        -------
        list
            A list of strings containing:
                1. the structure name
                2. lattice type
                3. hall symbol of the cristaline structure
                4. space group
                5. number of the space group,
                6. number of operation symmetry
        """

        self.topology = 'HCB_A'
        self.dimension = 2

        bb_triangular = Building_Block(name_bb_a,
                                       save_dir=self.lib_path,
                                       verbosity=self.verbosity)

        bb_linear = Building_Block(name_bb_b,
                                   save_dir=self.lib_path,
                                   verbosity=self.verbosity)

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
        struct.translate_sites(range(len(struct.as_dict()['sites'])), 
                              [0, 0, 0.5], 
                              frac_coords=True, 
                              to_unit_cell=True)

        # Simetriza a estrutura
        try:
            symm = SpacegroupAnalyzer(struct, symprec=0.3, angle_tolerance=3.0)
            struct_symm_prim = symm.get_refined_structure()
        except Exception:
            return None

        if stacking not in self.available_stacking[self.topology]:
            raise Exception(f"""{stacking} is not in the available stack list for HCB-A net.
    Available options are: {self.available_stacking[self.topology]}""")

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

        # Create AAl stacking. 
        # Hexagonal cell with two sheets per cell shifited by the shift_vector in angstroms.
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

        if print_result is True:
            Tools.print_result(self.name,
                               str(self.lattice_type),
                               str(self.hall[0:2]),
                               str(self.space_group),
                               str(self.space_group_n),
                               len(symm_op))

        return [self.name,
                str(self.lattice_type), 
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]


    def create_sql_structure(self,
                             name_bb_a: str,
                             name_bb_b: str,
                             stacking: str ='AA',
                             bond_atom: str ='N',
                             c_parameter_base: float = 3.6,
                             print_result: bool = True,
                             slab: float = 10.0,
                             shift_vector: list = [1.0, 1.0, 0],
                             tilt_angle: float = 5.0):
        """Creates a COF with SQL network.

        The SQL net is composed of two tetrapodal building blocks.

        Parameters
        ----------
        name_bb_a : str, required
            The 4 letter code for the tetrapodal Buiding Block A
        name_bb_b : str, required
            The 4 letter code for the tetrapodal Buiding Block B
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
        bond_atom : str, optional
            The atom that connect the buiding blocks (default is 'N')
        c_parameter_base : float, optional
            The base value for interlayer distance in angstroms (default is 3.6)
        print_result : bool, optional
            Parameter for the control for printing the result (default is True)
        slab : float, optional
            Default parameter for the interlayer slab (default is 10.0)
        shift_vector: list, optional
            Shift vector for the AAl and AAt stakings (defatult is [1.0,1.0,0])
        tilt_angle: float, optional
            Tilt angle for the AAt staking in degrees (default is 5.0)

        Returns
        -------
        list
            A list of strings containing:
                1. the structure name,
                2. lattice type,
                3. hall symbol of the cristaline structure,
                4. space group,
                5. number of the space group,
                6. number of operation symmetry
        """

        self.topology = 'SQL'
        self.dimension = 2

        bb_1 = Building_Block(name_bb_a,
                              save_dir=self.lib_path,
                              verbosity=self.verbosity)

        bb_2 = Building_Block(name_bb_b,
                              save_dir=self.lib_path,
                              verbosity=self.verbosity)

        self.charge = bb_1.charge + bb_2.charge
        self.chirality = bb_1.chirality or bb_2.chirality

        self.name = f'{bb_1.name}-{bb_2.name}-SQL-{stacking}'

        if self.verbosity is True:
            print('Starting the creation of a SQL net')

        if bb_1.connectivity != 4:
            print(f'Building block A ({name_bb_a}) must present connectivity 4 insted of',
                  bb_1.connectivity)
            return None
        if bb_2.connectivity != 4:
            print(f'Building block B ({name_bb_b}) must present connectivity 4 insted of',
            bb_2.connectivity)
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

        # Build the matrix of unitary hexagonal cell
        if stacking == 'A':
            c = slab
        else:
            c = c_parameter_base + max([delta_a, delta_b])

        # Constrói a matriz da célula unitária tetragonal
        lattice = [[a, 0, 0], [0, a, 0], [0, 0, c]]

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
        struct.translate_sites(range(len(struct.as_dict()['sites'])), 
                               [0, 0, 0.5], 
                               frac_coords=True, 
                               to_unit_cell=True)

        # Simetriza a estrutura
        symm = SpacegroupAnalyzer(struct, symprec=1, angle_tolerance=5.0)
        struct_symm_prim = symm.get_refined_structure()

        if stacking not in self.available_stacking[self.topology]:
                raise Exception(f"""{stacking} is not in the available stack list for HCB net.
    Available options are: {self.available_stacking[self.topology]}""")

        # Create A stacking, a 2D isolated sheet with slab
        if stacking == 'A':
            self.stacking = 'A'
            self.symm_structure = struct_symm_prim

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
            AB_2 = Structure(lattice, 
                             AB_label + AB_label,
                             AB,
                             coords_are_cartesian=False)

            AB_2_symm = SpacegroupAnalyzer(AB_2,
                                           symprec=self.symm_tol,
                                           angle_tolerance=self.angle_tol)

            self.symm_structure = AB_2_symm.get_refined_structure()

        # Create AAl stacking. Tetragonal cell with two sheets
        # per cell shifited by the shift_vector in angstroms.
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

        if print_result is True:
            Tools.print_result(self.name,
                               str(self.lattice_type),
                               str(self.hall[0:2]),
                               str(self.space_group),
                               str(self.space_group_n),
                               len(symm_op))

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]


    def create_sql_a_structure(self,
                               name_bb_a: str,
                               name_bb_b: str,
                               stacking: str = 'AA',
                               bond_atom: str = 'N',
                               c_parameter_base: float = 3.6,
                               print_result: bool = True,
                               slab: float = 10.0,
                               shift_vector: list = [1.0, 1.0, 0],
                               tilt_angle: float = 5.0):
        """Creates a COF with SQL-A network.

        The SQL-A net is composed of one tetrapodal and one linear building blocks.

        Parameters
        ----------
        name_bb_a : str, required
            The 4 letter code for the tripodal Buiding Block A
        name_bb_b : str, required
            The 4 letter code for the linear Buiding Block B
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
        bond_atom : str, optional
            The atom that connect the buiding blocks (default is 'N')
        c_parameter_base : float, optional
            The base value for interlayer distance in angstroms (default is 3.6)
        print_result : bool, optional
            Parameter for the control for printing the result (default is True)
        slab : float, optional
            Default parameter for the interlayer slab (default is 10.0)
        shift_vector: list, optional
            Shift vector for the AAl and AAt stakings (defatult is [1.0,1.0,0])
        tilt_angle: float, optional
            Tilt angle for the AAt staking in degrees (default is 5.0)

        Returns
        -------
        list
            A list of strings containing:
                1. the structure name,
                2. lattice type,
                3. hall symbol of the cristaline structure,
                4. space group,
                5. number of the space group,
                6. number of operation symmetry
        """

        self.topology = 'SQL_A'
        self.dimension = 2

        bb_1 = Building_Block(name_bb_a,
                              save_dir=self.lib_path,
                              verbosity=self.verbosity)

        bb_2 = Building_Block(name_bb_b,
                              save_dir=self.lib_path,
                              verbosity=self.verbosity)

        self.charge = bb_1.charge + bb_2.charge
        self.chirality = bb_1.chirality or bb_2.chirality

        self.name = f'{bb_1.name}-{bb_2.name}-SQL_A-{stacking}'

        if self.verbosity is True:
            print('Starting the creation of a sql net')

        if bb_1.connectivity != 4:
            print(f'Building block A ({name_bb_a}) must present connectivity 4 insted of', bb_1.connectivity)
            return None
        if bb_2.connectivity != 2:
            print(f'Building block B ({name_bb_b}) must present connectivity 2 insted of', bb_2.connectivity)
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

        # Build the matrix of unitary hexagonal cell
        if stacking == 'A':
            c = slab
        else:
            c = c_parameter_base + max([delta_a, delta_b])

        # Constrói a matriz da célula unitária tetragonal
        lattice = [[a, 0, 0], [0, a, 0], [0, 0, c]]

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

        if stacking not in self.available_stacking[self.topology]:
                raise Exception(f"""{stacking} is not in the available stack list for HCB net.
    Available options are: {self.available_stacking[self.topology]}""")

        # Create A stacking, a 2D isolated sheet with slab
        if stacking == 'A':
            self.stacking = 'A'
            self.symm_structure = struct_symm_prim
            
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
            self.symm_structure = Structure(lattice, 
                                            AB_label+AB_label, 
                                            AB, 
                                            coords_are_cartesian=False)

        # Create AA tilted stacking. 
        # Tilted tetragonal cell with two sheets per cell tilted by tilt_angle.
        if stacking == 'AAt':
            self.stacking = 'AAt'

            labels_conv_crystal = np.array([[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz'] for i in struct_symm_prim.as_dict()['sites']])
            cell = struct_symm_prim.as_dict()['lattice']

            # Shift the cell by the tilt angle
            print(cell)
            a_cell = cell['a']
            b_cell = cell['b']
            c_parameter_base = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta= cell['beta'] - tilt_angle
            gamma = cell['gamma']

            new_cell = Tools.cellpar_to_cell([a_cell, b_cell, c_parameter_base, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + np.array([0, 0, 0.5*c_parameter_base]) + shift_vector

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

        if print_result is True:
            Tools.print_result(self.name,
                               str(self.lattice_type),
                               str(self.hall[0:2]),
                               str(self.space_group),
                               str(self.space_group_n),
                               len(symm_op))

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]


    def create_kgd_structure(self, 
                             name_bb_a: str,
                             name_bb_b: str,
                             stacking: str ='AA',
                             bond_atom: str ='N',
                             c_parameter_base: float = 3.6,
                             print_result: bool = True,
                             slab: float = 10.0,
                             shift_vector: list = [1.0, 1.0, 0],
                             tilt_angle: float = 5.0):

        """Creates a COF with KGD network.

        The KGD net is composed of one hexapodal and one tripodal building blocks.

        Parameters
        ----------
        name_bb_a : str, required
            The 4 letter code for the hexapodal Buiding Block A
        name_bb_b : str, required
            The 4 letter code for the tripodal Buiding Block B
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
        bond_atom : str, optional
            The atom that connect the buiding blocks (default is 'N')
        c_parameter_base : float, optional
            The base value for interlayer distance in angstroms (default is 3.6)
        print_result : bool, optional
            Parameter for the control for printing the result (default is True)
        slab : float, optional
            Default parameter for the interlayer slab (default is 10.0)
        shift_vector: list, optional
            Shift vector for the AAl and AAt stakings (defatult is [1.0,1.0,0])
        tilt_angle: float, optional
            Tilt angle for the AAt staking in degrees (default is 5.0)

        Returns
        -------
        list
            A list of strings containing:
                1. the structure name,
                2. lattice type,
                3. hall symbol of the cristaline structure,
                4. space group,
                5. number of the space group,
                6. number of operation symmetry
        """

        self.topology = 'KGD'
        self.dimension = 2

        bb_1 = Building_Block(name_bb_a,
                              save_dir=self.lib_path,
                              verbosity=self.verbosity)

        bb_2 = Building_Block(name_bb_b,
                              save_dir=self.lib_path,
                              verbosity=self.verbosity)

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
        struct.translate_sites(range(len(struct.as_dict()['sites'])), 
                               [0, 0, 0.5],
                               frac_coords=True,
                               to_unit_cell=True)

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

        if print_result is True:
            Tools.print_result(self.name,
                               str(self.lattice_type),
                               str(self.hall[0:2]),
                               str(self.space_group),
                               str(self.space_group_n),
                               len(symm_op))

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]


    def create_hxl_a_structure(self,
                               name_bb_a: str,
                               name_bb_b: str,
                               stacking: str = 'AA',
                               bond_atom: str = 'N',
                               c_parameter_base: float = 3.6,
                               print_result: bool = True,
                               slab: float = 10.0,
                               shift_vector: list = [1.0, 1.0, 0],
                               tilt_angle: float = 5.0):

        """Creates a COF with HXL-A network.

        The HXK-A net is composed of one hexapodal and one linear building blocks.

        Parameters
        ----------
        name_bb_a : str, required
            The 4 letter code for the hexapodal Buiding Block A
        name_bb_b : str, required
            The 4 letter code for the linear Buiding Block B
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
        bond_atom : str, optional
            The atom that connect the buiding blocks (default is 'N')
        c_parameter_base : float, optional
            The base value for interlayer distance in angstroms (default is 3.6)
        print_result : bool, optional
            Parameter for the control for printing the result (default is True)
        slab : float, optional
            Default parameter for the interlayer slab (default is 10.0)
        shift_vector: list, optional
            Shift vector for the AAl and AAt stakings (defatult is [1.0,1.0,0])
        tilt_angle: float, optional
            Tilt angle for the AAt staking in degrees (default is 5.0)

        Returns
        -------
        list
            A list of strings containing:
                1. the structure name,
                2. lattice type,
                3. hall symbol of the cristaline structure,
                4. space group,
                5. number of the space group,
                6. number of operation symmetry
        """

        self.topology = 'HXL_A'
        self.dimension = 2

        bb_1 = Building_Block(name_bb_a,
                              save_dir=self.lib_path,
                              verbosity=self.verbosity)

        bb_2 = Building_Block(name_bb_b,
                              save_dir=self.lib_path,
                              verbosity=self.verbosity)

        self.name = f'{bb_1.name}-{bb_2.name}-HXL_A-{stacking}'
        self.charge = bb_1.charge + bb_2.charge
        self.chirality = bb_1.chirality or bb_2.chirality

        if self.verbosity is True:
            print(f'Starting the creation of a {self.topology} net')

        if bb_1.connectivity != 6:
            print('Building block A must present connectivity 6 insted of', len(bb_1.connectivity))
            return None
        if bb_2.connectivity != 2:
            print('Building block B must present connectivity 2 insted of', len(bb_2.connectivity))
            return None

        # Calculate the cell parameter based on the size of the building blocks
        size_a = max(bb_1.size)
        if self.verbosity:
            print('BB_A size:', size_a)

        size_b = max(bb_2.size)
        if self.verbosity:
            print('BB_B size:', size_b)

        # Defines the cell parameter a
        a = 2*(size_a + size_b) # np.cos(np.radians(30))*2*(size_a + size_b)

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
        final_pos = np.dot(bb_1.atom_pos, R.from_euler('z', 30, degrees=True).as_matrix())
        final_label = list(bb_1.atom_labels)

        # Add tbe building block 2 (C2) on [0.5, 0.5, 0.0] of the unitary cell (B1 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', 30, degrees=True).as_matrix()) + np.array([1/4, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_2.atom_labels)
        
        # Add tbe building block 2 (C2) on [0.5, 0.0, 0.0] of the unitary cell (B2 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', 90, degrees=True).as_matrix()) + np.array([0.5, 0, 0])*a))
        final_label += list(bb_2.atom_labels)
        
        # Add tbe building block 2 (C2) on [0.0, 0.5, 0.0] of the unitary cell (B3 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', -30, degrees=True).as_matrix()) + np.array([-1/4, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_2.atom_labels)

        # Changes the X atoms by the desired bond_atom
        final_label, final_pos = Tools.change_X_atoms(final_label, final_pos, bond_atom)
    
        # Cria a estrutura como entidade do pymatgen
        struct = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)
              
        # Remove os átomos duplicados
        struct.sort(reverse=True)
        struct.merge_sites(tol=.7, mode='delete')
        struct.translate_sites(range(len(struct.as_dict()['sites'])), 
                               [0, 0, 0.5], 
                               frac_coords=True, 
                               to_unit_cell=True)          
        
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

        if print_result is True:
            Tools.print_result(self.name,
                               str(self.lattice_type),
                               str(self.hall[0:2]),
                               str(self.space_group),
                               str(self.space_group_n),
                               len(symm_op))

        return [self.name, 
                str(self.lattice_type), 
                str(self.hall[0:2]), 
                str(self.space_group), 
                str(self.space_group_n), 
                len(symm_op)]


    def create_kgm_structure(self,
                             name_bb_a: str,
                             name_bb_b: str,
                             stacking: str ='AA',
                             bond_atom: str ='N',
                             c_parameter_base: float = 3.6,
                             print_result: bool = True,
                             slab: float = 10.0,
                             shift_vector: list = [1.0, 1.0, 0],
                             tilt_angle: float = 5.0):

        """Creates a COF with KGM network.

        The KGM net is composed of two tetrapodal building blocks.

        Parameters
        ----------
        name_bb_a : str, required
            The 4 letter code for the tetrapodal Buiding Block A
        name_bb_b : str, required
            The 4 letter code for the tetrapodal Buiding Block B
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
        bond_atom : str, optional
            The atom that connect the buiding blocks (default is 'N')
        c_parameter_base : float, optional
            The base value for interlayer distance in angstroms (default is 3.6)
        print_result : bool, optional
            Parameter for the control for printing the result (default is True)
        slab : float, optional
            Default parameter for the interlayer slab (default is 10.0)
        shift_vector: list, optional
            Shift vector for the AAl and AAt stakings (defatult is [1.0,1.0,0])
        tilt_angle: float, optional
            Tilt angle for the AAt staking in degrees (default is 5.0)

        Returns
        -------
        list
            A list of strings containing:
                1. the structure name,
                2. lattice type,
                3. hall symbol of the cristaline structure,
                4. space group,
                5. number of the space group,
                6. number of operation symmetry
        """

        self.topology = 'KGM'
        self.dimension = 2

        bb_1 = Building_Block(name_bb_a,
                              save_dir=self.lib_path,
                              verbosity=self.verbosity)

        bb_2 = Building_Block(name_bb_b,
                              save_dir=self.lib_path,
                              verbosity=self.verbosity)

        self.name = f'{bb_1.name}-{bb_2.name}-KGM-{stacking}'
        self.charge = bb_1.charge + bb_2.charge
        self.chirality = bb_1.chirality or bb_2.chirality

        if self.verbosity is True:
            print('Starting the creation of a KGM net')

        if bb_1.connectivity != 4:
            print('Building block A must present connectivity 4 insted of', len(bb_1.connectivity))
            return None
        if bb_2.connectivity != 4:
            print('Building block B must present connectivity 4 insted of', len(bb_2.connectivity))
            return None

        # Calculate the cell parameter based on the size of the building blocks        
        if self.verbosity:
            print('BB_A size:', bb_1.size)

        if self.verbosity:
            print('BB_B size:', bb_2.size)

        size_a = max(bb_1.size)
        size_b = max(bb_2.size)
        # Defines the cell parameter a
        a = 2*(size_a + size_b) # np.cos(np.radians(30))*2*(size_a + size_b)

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

        bb1_v1, bb1_v2, bb1_v3, bb1_v4 = bb_1.get_X_points()[1]
        bb1_lado1, bb1_lado2, bb1_lado3 = (bb1_v1 + bb1_v2) /2, (bb1_v1 + bb1_v3) /2, (bb1_v1 + bb1_v4) /2

        R_matrix = Tools.rotation_matrix_from_vectors(bb1_lado1, np.array([1, 0, 0]))

        bb_1.atom_pos = np.dot(bb_1.atom_pos, np.transpose(R_matrix))

        # Add tbe building block 1 (C4) on the center [1/2, 1/2, 0.0] of the unitary cell (A1 site) 
        final_pos = np.dot(bb_1.atom_pos, R.from_euler('z', 30, degrees=True).as_matrix()) + np.array([1/4, np.sqrt(3)/4, 0])*a
        final_label = list(bb_1.atom_labels)

        # Add tbe building block 1 (C4) on [1/2, 0.0, 0.0] of the unitary cell (A2 site) -45
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler('z', -90, degrees=True).as_matrix()) + np.array([1/2, 0, 0])*a))
        final_label += list(bb_1.atom_labels)

        # Add tbe building block 1 (C4) on [0.0, 1/2, 0.0] of the unitary cell (A3 site) 30
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler('z', -30, degrees=True).as_matrix()) + np.array([-1/4, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_1.atom_labels)

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
    
        # Simetriza a estrutura
        symm = SpacegroupAnalyzer(struct, symprec=self.symm_tol, angle_tolerance=self.angle_tol)
        struct_symm_prim = symm.get_refined_structure()

        if stacking not in self.available_stacking[self.topology]:
                raise Exception(f"""{stacking} is not in the available stack list for HCB net.
    Available options are: {self.available_stacking[self.topology]}""")
        
        # Create A stacking, a 2D isolated sheet with slab
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

        if print_result is True:
            Tools.print_result(self.name,
                               str(self.lattice_type),
                               str(self.hall[0:2]),
                               str(self.space_group),
                               str(self.space_group_n),
                               len(symm_op))

        return [self.name, 
                str(self.lattice_type), 
                str(self.hall[0:2]), 
                str(self.space_group), 
                str(self.space_group_n), 
                len(symm_op)]

    
    def create_kgm_a_structure(self,
                               name_bb_a: str,
                               name_bb_b: str,
                               stacking: str = 'AA',
                               bond_atom: str = 'N',
                               c_parameter_base: float = 3.6,
                               print_result: bool = True,
                               slab: float = 10.0,
                               shift_vector: list = [1.0, 1.0, 0],
                               tilt_angle: float = 5.0):

        """Creates a COF with KGM-A network.

        The KGM-A net is composed of one tetrapodal and one linear building blocks.

        Parameters
        ----------
        name_bb_a : str, required
            The 4 letter code for the tetrapodal Buiding Block A
        name_bb_b : str, required
            The 4 letter code for the linear Buiding Block B
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
        bond_atom : str, optional
            The atom that connect the buiding blocks (default is 'N')
        c_parameter_base : float, optional
            The base value for interlayer distance in angstroms (default is 3.6)
        print_result : bool, optional
            Parameter for the control for printing the result (default is True)
        slab : float, optional
            Default parameter for the interlayer slab (default is 10.0)
        shift_vector: list, optional
            Shift vector for the AAl and AAt stakings (defatult is [1.0,1.0,0])
        tilt_angle: float, optional
            Tilt angle for the AAt staking in degrees (default is 5.0)

        Returns
        -------
        list
            A list of strings containing:
                1. the structure name,
                2. lattice type,
                3. hall symbol of the cristaline structure,
                4. space group,
                5. number of the space group,
                6. number of operation symmetry
        """

        self.topology = 'KGM_A'
        self.dimension = 2

        bb_1 = Building_Block(name_bb_a,
                              save_dir=self.lib_path,
                              verbosity=self.verbosity)

        bb_2 = Building_Block(name_bb_b,
                              save_dir=self.lib_path,
                              verbosity=self.verbosity)

        self.name = f'{bb_1.name}-{bb_2.name}-KGM_A-{stacking}'
        self.charge = bb_1.charge + bb_2.charge
        self.chirality = bb_1.chirality or bb_2.chirality

        if self.verbosity is True:
            print('Starting the creation of a KGM_A net')

        if bb_1.connectivity != 4:
            print('Building block A must present connectivity 4 insted of', len(bb_1.connectivity))
            return None
        if bb_2.connectivity != 2:
            print('Building block B must present connectivity 2 insted of', len(bb_2.connectivity))
            return None

        # Calculate the cell parameter based on the size of the building blocks
        size_a = max(bb_1.size)
        if self.verbosity:
            print('BB_A size:', size_a)

        size_b = max(bb_2.size)
        if self.verbosity:
            print('BB_B size:', size_b)

        # Defines the cell parameter a
        a = 4*(size_a + size_b) # np.cos(np.radians(30))*2*(size_a + size_b)

        if self.verbosity:
            print('Calculated cell parameter a:', a)

        # Gets the maximum distace in the z axis to create the c parameter
        delta_a = abs(max(np.transpose(bb_1.atom_pos)[2])) + abs(min(np.transpose(bb_1.atom_pos)[2]))
        delta_b = abs(max(np.transpose(bb_2.atom_pos)[2])) + abs(min(np.transpose(bb_2.atom_pos)[2]))

        # Build the matrix of unitary hexagonal cell
        if stacking == 'A':
            c = max([slab, max([delta_a, delta_b])])
        else:
            c = c_parameter_base + max([delta_a, delta_b])

        # Define the cell lattice as an hexagonal cell
        lattice = [[a, 0, 0],
                   [-0.5*a, np.sqrt(3)/2*a, 0],
                   [0, 0, c]]

        if self.verbosity is True:
            print('Unitary cell built:', lattice)
        
        bb1_v1, bb1_v2, bb1_v3, bb1_v4 = bb_1.get_X_points()[1]
        bb1_lado1, bb1_lado2, bb1_lado3 = (bb1_v1 + bb1_v2) /2, (bb1_v1 + bb1_v3) /2, (bb1_v1 + bb1_v4) /2

        R_matrix = Tools.rotation_matrix_from_vectors(bb1_lado1, np.array([1, 0, 0]))

        bb_1.atom_pos = np.dot(bb_1.atom_pos, np.transpose(R_matrix))

        # Add tbe building block 1 (C4) on the center [1/2, 1/2, 1/2] of the unitary cell (A1 site) 
        final_pos = np.dot(bb_1.atom_pos, R.from_euler('z', 30, degrees=True).as_matrix()) + np.array([1/4, np.sqrt(3)/4, 0])*a
        final_label = list(bb_1.atom_labels)

        # Add tbe building block 1 (C4) on [1/2, 0.0, 0.0] of the unitary cell (A2 site)
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler('z', -90, degrees=True).as_matrix()) + np.array([1/2, 0, 0])*a))
        final_label += list(bb_1.atom_labels)

        # Add tbe building block 1 (C4) on [0.0, 1/2, 0.0] of the unitary cell (A3 site)
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler('z', -30, degrees=True).as_matrix()) + np.array([-1/4, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_1.atom_labels)

        '''
        Olde version of the positioning code
        # Add tbe building block 1 (C4) on the center of the unitary cell (A1 site) 
        final_pos = np.dot(bb_1.atom_pos, R.from_euler('z', -30, degrees=True).as_matrix()) + np.array([1/4, np.sqrt(3)/4, 0])*a
        final_label = list(bb_1.atom_labels)

        # Add tbe building block 1 (C4) on [0.5, 0.0, 0.0] of the unitary cell (A2 site)
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler('z', -60, degrees=True).as_matrix()) + np.array([1/2, 0, 0])*a))
        final_label += list(bb_1.atom_labels)

        # Add tbe building block 1 (C4) on [0.0, 0.5, 0.0] of the unitary cell (A3 site)
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler('z', 60, degrees=True).as_matrix()) + np.array([-1/4, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_1.atom_labels)'''

        # Add the building block 2 (C2) on [1/2, 1/4, 0.0] of the unitary cell (B1 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', -30, degrees=True).as_matrix()) + np.array([3/8, np.sqrt(3)/8, 0])*a))
        final_label += list(bb_2.atom_labels)

        # Add the building block 2 (C2) on [1/4, 3/4, 0.0] of the unitary cell (B2 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', -30, degrees=True).as_matrix()) + np.array([1/8, 3*np.sqrt(3)/8, 0])*a))
        final_label += list(bb_2.atom_labels)

        # Add tbe building block 2 (C2) on [3/4, 1/4, 0.0] of the unitary cell (B3 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', 30, degrees=True).as_matrix()) + np.array([5/8, np.sqrt(3)/8, 0])*a))
        final_label += list(bb_2.atom_labels)

        # Add tbe building block 2 (C2) on [1/4, 3/4, 0.0] of the unitary cell (B4 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', 30, degrees=True).as_matrix()) + np.array([-1/8, 3*np.sqrt(3)/8, 0])*a))
        final_label += list(bb_2.atom_labels)
        
        # Add tbe building block 2 (C2) on [3/4, 1/2, 0.0] of the unitary cell (B5 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', 90, degrees=True).as_matrix()) + np.array([4/8, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_2.atom_labels)

        # Add tbe building block 2 (C2) on [1/2, 3/4, 0.0] of the unitary cell (B6 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler('z', 90, degrees=True).as_matrix()) + np.array([0.0, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_2.atom_labels)

        # Change the X atoms by the desired bond_atom
        final_label, final_pos = Tools.change_X_atoms(final_label, final_pos, bond_atom)

        # Creates a pymatgen structure
        struct = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)

        # Remove dupolicate atoms
        struct.sort(reverse=True)
        struct.merge_sites(tol=.7, mode='delete')
        struct.translate_sites(range(len(struct.as_dict()['sites'])),
                               [0, 0, 0.5], 
                               frac_coords=True, 
                               to_unit_cell=True)

        if stacking not in self.available_stacking[self.topology]:
                raise Exception(f"""{stacking} is not in the available stack list for HCB net.
        Available options are: {self.available_stacking[self.topology]}""")

        # Create A stacking. The slab is defined by the c_cell parameter
        if stacking == 'A':
            self.stacking = 'A'
            symm = SpacegroupAnalyzer(struct,
                                      symprec=self.symm_tol,
                                      angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        # Create AA staking. By default one sheet per unitary cell is used
        if stacking == 'AA':
            self.stacking = 'AA'
            symm = SpacegroupAnalyzer(struct, 
                                      symprec=self.symm_tol,
                                      angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        if stacking == 'AB1x':
            self.stacking = 'AB1x'

            lattice = Lattice(np.array(struct.as_dict()['lattice']['matrix'])*(1, 1, 2))

            A = struct
            A_labels = np.array([i['label'] for i in A.as_dict()['sites']])
            A_pos = np.array([i['xyz'] for i in A.as_dict()['sites']])

            B = struct
            B_labels = np.array([i['label'] for i in B.as_dict()['sites']])
            B_pos = np.array([i['xyz'] for i in B.as_dict()['sites']])
            B_pos += np.array([lattice.a*1/2, 0, lattice.c*1/2])

            AB_1 = Structure(lattice, 
                             np.concatenate((A_labels, B_labels)),
                             np.concatenate((A_pos, B_pos)),
                             coords_are_cartesian=True)

            dict_structure = AB_1.as_dict()

            self.lattice = dict_structure['lattice']['matrix']

            self.atom_labels = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        if stacking == 'AB1y':
            self.stacking = 'AB1y'

            lattice = Lattice(np.array(struct.as_dict()['lattice']['matrix'])*(1, 1, 2))

            A = struct
            A_labels = np.array([i['label'] for i in A.as_dict()['sites']])
            A_pos = np.array([i['xyz'] for i in A.as_dict()['sites']])

            B = struct
            B_labels = np.array([i['label'] for i in B.as_dict()['sites']])
            B_pos = np.array([i['xyz'] for i in B.as_dict()['sites']])
            B_pos += np.array([lattice.a*1/4, lattice.b*1/4*np.sqrt(3), lattice.c*1/2])

            AB_1 = Structure(lattice, 
                             np.concatenate((A_labels, B_labels)),
                             np.concatenate((A_pos, B_pos)), 
                             coords_are_cartesian=True)

            dict_structure = AB_1.as_dict()

            self.lattice = dict_structure['lattice']['matrix']

            self.atom_labels = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        if stacking == 'AB1xy':
            self.stacking = 'AB1xy'

            lattice = Lattice(np.array(struct.as_dict()['lattice']['matrix'])*(1, 1, 2))

            A = struct
            A_labels = np.array([i['label'] for i in A.as_dict()['sites']])
            A_pos = np.array([i['xyz'] for i in A.as_dict()['sites']])

            B = struct
            B_labels = np.array([i['label'] for i in B.as_dict()['sites']])
            B_pos = np.array([i['xyz'] for i in B.as_dict()['sites']])
            B_pos += np.array([lattice.a*1/4, -lattice.b*1/4*np.sqrt(3), lattice.c*1/2])

            AB_1 = Structure(lattice, 
                             np.concatenate((A_labels, B_labels)),
                             np.concatenate((A_pos, B_pos)), 
                             coords_are_cartesian=True)

            dict_structure = AB_1.as_dict()

            self.lattice = dict_structure['lattice']['matrix']

            self.atom_labels = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        if stacking == 'AB2':
            self.stacking = 'AB2'

            lattice = Lattice(np.array(struct.as_dict()['lattice']['matrix'])*(1, 1, 2))

            A = struct
            A_labels = np.array([i['label'] for i in A.as_dict()['sites']])
            A_pos = np.array([i['xyz'] for i in A.as_dict()['sites']])

            B = struct
            B_labels = np.array([i['label'] for i in B.as_dict()['sites']])
            B_pos = np.array([i['xyz'] for i in B.as_dict()['sites']])

            B_pos += np.array([lattice.a*1/2, 
                               lattice.b*np.sqrt(3)/8, 
                               lattice.c*1/2])

            AB_1 = Structure(lattice, 
                             np.concatenate((A_labels, B_labels)),
                             np.concatenate((A_pos, B_pos)), 
                             coords_are_cartesian=True)

            dict_structure = AB_1.as_dict()

            self.lattice = dict_structure['lattice']['matrix']

            self.atom_labels = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        # Create AAl stacking.
        # Tetragonal cell with two sheets per cell shifited by the shift_vector in angstroms. 
        if stacking == 'AAl':
            self.stacking = 'AAl'

            lattice = Lattice(np.array(struct.as_dict()['lattice']['matrix'])*(1, 1, 2))

            A = struct
            A_labels = np.array([i['label'] for i in A.as_dict()['sites']])
            A_pos = np.array([i['xyz'] for i in A.as_dict()['sites']])

            B = struct
            B_labels = np.array([i['label'] for i in B.as_dict()['sites']])
            B_pos = np.array([i['xyz'] for i in B.as_dict()['sites']])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(*lattice.parameters)
            shift_vector = np.dot(r, np.array(shift_vector))

            B_pos += np.array([lattice.a, lattice.b*1/2*np.sqrt(3), 1]) * shift_vector
            B_pos += np.array([0, 0, lattice.c*1/2])

            AB_1 = Structure(lattice, 
                             np.concatenate((A_labels, B_labels)),
                             np.concatenate((A_pos, B_pos)), 
                             coords_are_cartesian=True)

            dict_structure = AB_1.as_dict()

            self.lattice = dict_structure['lattice']['matrix']

            self.atom_labels = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        # Create AA tilted stacking. Tilted tetragonal cell with two sheets per cell tilted by tilt_angle. 
        if stacking == 'AAt':
            self.stacking = 'AAt'

            a_cell = struct.as_dict()['lattice']['a']
            b_cell = struct.as_dict()['lattice']['b']
            c_cell = struct.as_dict()['lattice']['c'] * 2
            alpha = struct.as_dict()['lattice']['alpha'] - tilt_angle
            beta= struct.as_dict()['lattice']['beta'] - tilt_angle
            gamma = struct.as_dict()['lattice']['gamma']

            new_cell = Tools.cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])

            lattice = Lattice(new_cell)

            A = struct
            A_labels = np.array([i['label'] for i in A.as_dict()['sites']])
            A_pos = np.array([i['xyz'] for i in A.as_dict()['sites']])

            B = struct
            B_labels = np.array([i['label'] for i in B.as_dict()['sites']])
            B_pos = np.array([i['xyz'] for i in B.as_dict()['sites']])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(*lattice.parameters)
            shift_vector = np.dot(r, np.array(shift_vector))

            B_pos += np.array([lattice.a, lattice.b*1/2*np.sqrt(3), 1]) * shift_vector
            B_pos += np.array([0, 0, lattice.c*1/2])

            AB_1 = Structure(lattice, 
                             np.concatenate((A_labels, B_labels)),
                             np.concatenate((A_pos, B_pos)), 
                             coords_are_cartesian=True)

            dict_structure = AB_1.as_dict()

            self.lattice = dict_structure['lattice']['matrix']

            self.atom_labels = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        dict_structure = self.symm_structure.as_dict()

        self.lattice = dict_structure['lattice']['matrix']

        self.atom_labels = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.n_atoms = len(self.symm_structure)
        self.composition = self.symm_structure.formula

        if self.verbosity is True:
            print(self.symm_structure)

        # Get the simmetry information of the generated structure
        self.lattice_type =  symm.get_lattice_type()
        self.space_group =  symm.get_space_group_symbol()
        self.space_group_n =  symm.get_space_group_number()

        symm_op =  symm.get_point_group_operations()
        self.hall =  symm.get_hall()

        if print_result is True:
            Tools.print_result(self.name,
                               str(self.lattice_type),
                               str(self.hall[0:2]),
                               str(self.space_group),
                               str(self.space_group_n),
                               len(symm_op))

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]


################   Saving methods  ############################

    def save_cif(self, supercell : tuple = (1, 1, 1), path: str = None):
        """Save the structure in .cif format

        Parameters
        ----------
        supercell : tuple, optional
            List containing the supercell parameters.
            Default  = [1, 1, 1]
        path : str, optional
            Path to save the .cif file.
        """

        if path is not None:
            self.out_path = path

        os.makedirs(self.out_path, exist_ok=True)

        self.symm_structure.make_supercell(supercell)

        self.symm_structure.to(filename=os.path.join(self.out_path, self.name + '.cif'))

    def save_json(self, supercell : tuple = (1, 1, 1), path: str = None):
        """Save the structure in .json format

        Parameters
        ----------
        supercell : tuple, optional
            List containing the supercell parameters.
            Default  = [1, 1, 1]
        path : str, optional
            Path to save the .json file.
        """

        if path is not None:
            self.out_path = path

        os.makedirs(self.out_path, exist_ok=True)

        self.symm_structure.make_supercell(supercell)

        dict_sctructure = self.symm_structure.as_dict()

        a = dict_sctructure['lattice']['a']
        b = dict_sctructure['lattice']['b']
        c = dict_sctructure['lattice']['c']

        alpha = round(dict_sctructure['lattice']['alpha'], 3)
        beta = round(dict_sctructure['lattice']['beta'], 3)
        gamma = round(dict_sctructure['lattice']['gamma'], 3)

        atom_labels = [i['label'] for i in dict_sctructure['sites']]

        atom_pos = [i['xyz'] for i in dict_sctructure['sites']]

        Tools.save_json(self.out_path,
                        self.name,
                        [a, b, c, alpha, beta, gamma],
                        atom_labels, atom_pos)

    def save_xsf(self, supercell : tuple = (1, 1, 1), path: str = None):
        """Save the structure in XCrysDen .xsf format

        Parameters
        ----------
        supercell : tuple, optional
            List containing the supercell parameters.
            Default  = [1, 1, 1]
        path : str, optional
            Path to save the .xsf file.
        """

        if path is not None:
            self.out_path = path

        os.makedirs(self.out_path, exist_ok=True)

        self.symm_structure.make_supercell(supercell)

        dict_sctructure = self.symm_structure.as_dict()

        a = dict_sctructure['lattice']['a']
        b = dict_sctructure['lattice']['b']
        c = dict_sctructure['lattice']['c']

        alpha = round(dict_sctructure['lattice']['alpha'], 3)
        beta = round(dict_sctructure['lattice']['beta'], 3)
        gamma = round(dict_sctructure['lattice']['gamma'], 3)

        atom_labels = [i['label'] for i in dict_sctructure['sites']]

        atom_pos = [i['xyz'] for i in dict_sctructure['sites']]

        Tools.save_xsf(self.out_path,
                       self.name,
                       [a, b, c, alpha, beta, gamma],
                       atom_labels, atom_pos)

    def save_pdb(self, supercell : tuple = (1, 1, 1), path: str = None):
        """Save the structure in .pdb format

        Parameters
        ----------
        supercell : tuple, optional
            List containing the supercell parameters.
            Default  = [1, 1, 1]
        path : str, optional
            Path to save the .pdb file.
        """
        if path is not None:
            self.out_path = path

        os.makedirs(self.out_path, exist_ok=True)

        self.symm_structure.make_supercell(supercell)

        dict_sctructure = self.symm_structure.as_dict()

        a = dict_sctructure['lattice']['a'] 
        b = dict_sctructure['lattice']['b'] 
        c = dict_sctructure['lattice']['c']
        alpha = round(dict_sctructure['lattice']['alpha'], 3)
        beta = round(dict_sctructure['lattice']['beta'], 3)
        gamma = round(dict_sctructure['lattice']['gamma'], 3)

        atom_labels = [i['label'] for i in dict_sctructure['sites']]

        atom_pos = [i['xyz'] for i in dict_sctructure['sites']]

        Tools.save_pdb(self.out_path, self.name, [a, b, c, alpha, beta, gamma], atom_labels, atom_pos)


    def save_vasp(self, supercell : tuple = (1, 1, 1), path: str = None):
        """Save the structure in VASP POSCAR .vasp format

        Parameters
        ----------
        supercell : tuple, optional
            List containing the supercell parameters. 
            Default  = [1, 1, 1]
        path : str, optional
            Path to save the .vasp file.        
        """
        if path is not None:
            self.out_path = path

        os.makedirs(self.out_path, exist_ok=True)

        self.symm_structure.make_supercell(supercell)

        file_name = os.path.join(self.out_path, self.name + '.vasp')
        self.symm_structure.to(fmt='poscar', filename=file_name)

    def save_turbomole(self, supercell : tuple = (1, 1, 1), path: str = None):
        """Save the structure in Turbomole .coord format

        Parameters
        ----------
        supercell : tuple, optional
            List containing the supercell parameters. 
            Default  = [1, 1, 1]
        path : str, optional
            Path to save the .coord file.        
        """
        if path is not None:
            self.out_path = path

        os.makedirs(self.out_path, exist_ok=True)

        self.symm_structure.make_supercell(supercell)

        dict_sctructure = self.symm_structure.as_dict()

        a = dict_sctructure['lattice']['a'] 
        b = dict_sctructure['lattice']['b'] 
        c = dict_sctructure['lattice']['c']

        alpha = round(dict_sctructure['lattice']['alpha'], 3)
        beta = round(dict_sctructure['lattice']['beta'], 3)
        gamma = round(dict_sctructure['lattice']['gamma'], 3)

        atom_labels = [i['label'] for i in dict_sctructure['sites']]

        atom_pos = [i['xyz'] for i in dict_sctructure['sites']]

        temp_file = open(os.path.join(self.out_path ,self.name + '.coord'), 'w')
        temp_file.write('$coord angs\n')

        for i in range(len(atom_labels)):
            temp_file.write('{:>15.7f}{:>15.7f}{:>15.7f}   {:<5s}\n'.format(atom_pos[i][0], 
                                                                            atom_pos[i][1],
                                                                            atom_pos[i][2],
                                                                            atom_labels[i]))

        temp_file.write('$periodic 3\n')
        temp_file.write('$cell\n')
        temp_file.write(f'{a}  {b}  {c}  {alpha}  {beta}  {gamma}\n')
        temp_file.write('$opt\n')
        temp_file.write('   engine=inertial\n')
        temp_file.write('$end\n')

        temp_file.close()

    def save_xyz(self, supercell : tuple = (1, 1, 1), path: str = None):
        """Save the structure in .xyz format

        Parameters
        ----------
        supercell : tuple, optional
            List containing the supercell parameters. 
            Default  = [1, 1, 1]
        path : str, optional
            Path to save the .xyz file.        
        """
        if path is not None:
            self.out_path = path

        os.makedirs(self.out_path, exist_ok=True)

        self.symm_structure.make_supercell(supercell)

        dict_sctructure = self.symm_structure.as_dict()

        a = dict_sctructure['lattice']['a'] 
        b = dict_sctructure['lattice']['b'] 
        c = dict_sctructure['lattice']['c']

        alpha = round(dict_sctructure['lattice']['alpha'], 3)
        beta = round(dict_sctructure['lattice']['beta'], 3)
        gamma = round(dict_sctructure['lattice']['gamma'], 3)

        atom_labels = [i['label'] for i in dict_sctructure['sites']]

        atom_pos = [i['xyz'] for i in dict_sctructure['sites']]

        temp_file = open(os.path.join(self.out_path, self.name + '.xyz'), 'w')
        temp_file.write(f'{len(atom_labels)} \n')

        temp_file.write(f'{a}  {b}  {c}  {alpha}  {beta}  {gamma}\n')

        for i in range(len(atom_labels)):
            temp_file.write('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atom_labels[i], 
                                                                         atom_pos[i][0], 
                                                                         atom_pos[i][1], 
                                                                         atom_pos[i][2]))

        temp_file.close()

    def save_qe(self,
                supercell : tuple = (1, 1, 1),
                path: str = None,
                in_angstrom_format : bool = True,
                ecut : int = 40,
                erho : int = 360,
                k_dist : float = 0.3):
        """Save the structure in QuantumESPRESSO .in format

        Parameters
        ----------
        supercell : tuple, optional
            List containing the supercell parameters.
            Default = [1, 1, 1]
        path : str, optional
            Path to save the .in file.
        in_angstrom_format : bool, optional
            Save the atomic positions in angstroms insted of crystal coordinates
            Default = True
        ecut : int, optional
            Cutoff for wavefunctions. Default = 40
        erho : int, optional
            Cutoff for electronic density. Default = 360
        k_dist : 0.3, optional
            Distance of the k-points on the reciprocal space.
            Default = 0.3
        """
        if path is not None:
            self.out_path = path

        os.makedirs(self.out_path, exist_ok=True)

        self.symm_structure.make_supercell(supercell)

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

        if in_angstrom_format is False:
            ion_pos = ion_conv_crystal
        if in_angstrom_format is True:
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
        if in_angstrom_format is False:
            coords_type = 'crystal'
        out_file.write(f'ATOMIC_POSITIONS ({coords_type})\n')

        for atom in ion_pos:
            out_file.write('{:<5s}{:>15.9f}{:>15.9f}{:>15.9f}\n'.format(atom[0],
                                                                        atom[1],
                                                                        atom[2],
                                                                        atom[3]))

        out_file.write('K_POINTS automatic\n')
        out_file.write(f'   {kx} {ky} {kz}  1 1 1\n')

        out_file.write('EOF\n')
        out_file.write('$CMD_PW < $PREFIX.$CALC.in > $PREFIX.$CALC.out')

        out_file.close()
