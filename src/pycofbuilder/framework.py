# -*- coding: utf-8 -*-
# Copyright (c) Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

"""
This class implements definitions for Reticulum buiding
"""

import os
import numpy as np

# Import pymatgen
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from scipy.spatial.transform import Rotation as R

import pycofbuilder.tools as Tools
from pycofbuilder.building_block import Building_Block

from pycofbuilder.data.topology import TOPOLOGY_DICT


class Framework():
    """
    A class used to represent a Covalent Organic Framework as a reticular entity.

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
    from_name()
        Creates a reticulum from a name
    from_building_blocks()
        Creates a reticulum from two building blocks
    """

    def __init__(self, name=None, verbosity=False, out_dir=None):

        _ROOTDIR = os.path.abspath(os.path.dirname(__file__))

        self.verbosity = verbosity
        self.main_path = os.path.join(_ROOTDIR, 'data')

        if out_dir is None:
            self.out_path = os.path.join(os.getcwd(), 'out')
        else:
            self.out_path = out_dir

        self.lib_path = os.path.join(self.out_path, 'building_blocks')

        self.BB1_name = None
        self.BB2_name = None
        self.topology = None
        self.stacking = None
        self.smiles = None

        self.atom_labels = []
        self.atom_pos = []
        self.lattice = np.eye(3)
        self.lattice_sgs = None
        self.space_group = None
        self.space_group_n = None

        self.dimention = None
        self.n_atoms = self._get_n_atoms()
        self.mass = None
        self.composition = None
        self.charge = 0
        self.multiplicity = 1
        self.chirality = False

        self.symm_tol = 0.1
        self.angle_tol = 0.1

        # Falta adicionar: 'HXL', 'KGD_A'
        self.available_2D_topologies = ['HCB', 'HCB_A',
                                        'SQL', 'SQL_A',
                                        'KGM', 'KGM_A',
                                        'KGD',
                                        'HXL_A',
                                        'FXT_A']

        # Falta add: ['dia', 'bor', 'srs', 'pts', 'ctn', 'rra', 'fcc', 'lon', 'stp', 'acs', 'tbo', 'bcu', 'fjh', 'ceq']
        self.available_3D_topologies = ['DIA', 'BOR']  # Temporary
        self.available_topologies = self.available_2D_topologies + \
            self.available_3D_topologies

        # Define available stackings for all 2D topologies
        self.available_stacking = {'HCB': ['AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'HCB_A': ['AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'SQL': ['AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'SQL_A': ['AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'KGD': ['A', 'AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'HXL_A': ['A', 'AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'KGM': ['AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
                                   'KGM_A': ['AA', 'AB1x', 'AB1y', 'AB1xy', 'AB2', 'AAl', 'AAt'],
                                   'FXT_A': ['AA', 'AB1x', 'AB1y', 'AB1xy', 'AB2', 'AAl', 'AAt'],
                                   'DIA': [1, 2, 3, 4],  # Temporary
                                   'BOR': [5, 8, 6, 7]  # Temporary
                                   }

        if name is None:
            self.name = ""
        else:
            self.name = name
            self.from_name(self.name)

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return f'Reticulum({self.BB1_name}, {self.BB2_name}, {self.topology}, {self.stacking})'

    def _get_n_atoms(self) -> int:
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

    def check_name_concistency(self, FrameworkName) -> tuple[str, str, str, str]:
        """Checks if the name is in the correct format."""

        string_error = 'FrameworkName must be in the format: BB1_BB2_Net_Stacking'
        assert isinstance(FrameworkName, str), string_error

        name_error = 'FrameworkName must be in the format: BB1_BB2_Net_Stacking'
        assert len(FrameworkName.split('-')) == 4, name_error

        BB1_name, BB2_name, Net, Stacking = FrameworkName.split('-')

        net_error = f'Net must be one of the following: {self.available_topologies}'
        assert Net in self.available_topologies, net_error

        stacking_error = f'Stacking must be one of the following: {self.available_stacking[Net]}'
        assert Stacking in self.available_stacking[Net], stacking_error

        return BB1_name, BB2_name, Net, Stacking

    def from_name(self, FrameworkName, **kwargs) -> None:
        """Creates a COF from a given FrameworkName.

        Parameters
        ----------
        FrameworkName : str, required
            The name of the COF to be created

        Returns
        -------
        COF
            The COF object
        """
        BB1_name, BB2_name, Net, Stacking = self.check_name_concistency(FrameworkName)

        BB1 = Building_Block(name=BB1_name, save_dir=self.out_path)
        BB2 = Building_Block(name=BB2_name, save_dir=self.out_path)

        self.from_building_blocks(BB1, BB2, Net, Stacking, **kwargs)

    def from_building_blocks(self, BB1, BB2, Net, Stacking, **kwargs):
        """Creates a COF from the building blocks.

        Parameters
        ----------
        BB1 : BuildingBlock, required
            The first building block
        BB2 : BuildingBlock, required
            The second building block
        Net : str, required
            The network of the COF
        Stacking : str, required
            The stacking of the COF

        Returns
        -------
        COF
            The COF object
        """
        self.name = f'{BB1.name}-{BB2.name}-{Net}-{Stacking}'
        self.BB1_name = BB1.name
        self.BB2_name = BB2.name
        self.topology = Net
        self.stacking = Stacking

        # Check if the BB1 has the smiles attribute
        if hasattr(BB1, 'smiles') and hasattr(BB2, 'smiles'):
            self.smiles = f'{BB1.smiles}.{BB2.smiles}'
        else:
            print('WARNING: The smiles attribute is not available for the building blocks')

        net_build_dict = {
            'HCB': self._create_hcb_structure,
            'HCB_A': self._create_hcb_a_structure,
            'SQL': self._create_sql_structure,
            'SQL_A': self._create_sql_a_structure,
            'KGM': self._create_kgm_structure,
            'KGM_A': self._create_kgm_a_structure,
            'KGD': self._create_kgd_structure,
            'HXL_A': self._create_hxl_a_structure,
            # 'FXT': self._create_fxt_structure,
            'FXT_A': self._create_fxt_a_structure,
            }

        result = net_build_dict[Net](BB1, BB2, stacking=Stacking, **kwargs)

        return result

    def save(self, fmt: str = 'cif', supercell: list = [1, 1, 1], save_dir=None, **kwargs) -> None:
        '''
        Save the structure in a specif file format.

        Parameters
        ----------
        save_format : str, optional
            The file format to be saved
            Can be `json`, `cif`, `xyz`, `turbomole`, `vasp`, `xsf`, `pdb`, `pqr`, `qe`.
            Default: 'cif'
        supercell : list, optional
            The supercell to be used to save the structure.
            Default: [1,1,1]
        '''

        save_dict = {
            'json': Tools.save_json,
            'cjson': Tools.save_chemjson,
            'cif': Tools.save_cif,
            'xyz': Tools.save_xyz,
            'turbomole': Tools.save_turbomole,
            'vasp': Tools.save_vasp,
            'xsf': Tools.save_xsf,
            'pdb': Tools.save_pdb,
            'pqr': Tools.save_pqr,
            'qe': Tools.save_qe
         }

        file_format_error = f'Format must be one of the following: {save_dict.keys()}'
        assert fmt in save_dict.keys(), file_format_error

        structure = self.symm_structure
        structure.make_supercell(supercell)
        structure_dict = structure.as_dict()

        cell = structure_dict['lattice']['matrix']
        atom_pos = [site['xyz'] for site in structure_dict['sites']]
        atom_labels = [site['species'][0]['element'] for site in structure_dict['sites']]

        if save_dir is None:
            save_path = self.out_path
        else:
            save_path = save_dir

        save_dict[fmt](path=save_path,
                       file_name=self.name,
                       cell=cell,
                       atom_labels=atom_labels,
                       atom_pos=atom_pos,
                       **kwargs)

# --------------- Net creation methods -------------------------- #

    def _create_hcb_structure(self,
                              BB_T3_A,
                              BB_T3_B,
                              stacking: str = 'AA',
                              print_result: bool = True,
                              slab: float = 10.0,
                              shift_vector: list = [1.0, 1.0, 0],
                              tilt_angle: float = 5.0):
        """Creates a COF with HCB network.

        The HCB net is composed of two tripodal building blocks.

        Parameters
        ----------
        BB_T3_1 : BuildingBlock, required
            The BuildingBlock object of the tripodal Buiding Block A
        BB_T3_2 : BuildingBlock, required
            The BuildingBlock object of the tripodal Buiding Block B
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
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

        connectivity_error = 'Building block {} must present connectivity {}'
        assert BB_T3_A.connectivity == 3, connectivity_error.format('A', 3)
        assert BB_T3_B.connectivity == 3, connectivity_error.format('B', 3)

        self.topology = 'HCB'
        self.dimension = 2

        self.charge = BB_T3_A.charge + BB_T3_B.charge
        self.chirality = BB_T3_A.chirality or BB_T3_B.chirality

        self.name = f'{BB_T3_A.name}-{BB_T3_B.name}-HCB-{stacking}'

        self.charge = BB_T3_A.charge + BB_T3_B.charge
        self.chirality = BB_T3_A.chirality or BB_T3_B.chirality

        Tools.print_comand(f'Starting the creation of {self.name}',
                           self.verbosity,
                           ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = Tools.get_bond_atom(BB_T3_A.conector, BB_T3_B.conector)

        Tools.print_comand(f'Bond atom detected: {bond_atom}',
                           self.verbosity,
                           ['debug', 'high'])

        # Get the topology information
        topology_info = TOPOLOGY_DICT[self.topology]

        # Measure the base size of the building blocks
        size = BB_T3_A.size[0] + BB_T3_B.size[0]

        # Calculate the delta size to add to the c parameter
        delta_a = abs(max(np.transpose(BB_T3_A.atom_pos)[2])) + abs(min(np.transpose(BB_T3_A.atom_pos)[2]))
        delta_b = abs(max(np.transpose(BB_T3_B.atom_pos)[2])) + abs(min(np.transpose(BB_T3_B.atom_pos)[2]))

        delta_max = max([delta_a, delta_b])

        # Calculate the cell parameters
        a = topology_info['a'] * size
        b = topology_info['b'] * size
        c = topology_info['c'] + delta_max
        alpha = topology_info['alpha']
        beta = topology_info['beta']
        gamma = topology_info['gamma']

        # Create the lattice
        lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        final_label = []
        final_pos = []

        # Add the building blocks to the structure
        vertice_data = topology_info['vertices'][0]
        final_label += BB_T3_A.atom_labels
        vertice_pos = np.array(vertice_data['position'])*a

        R_Matrix = R.from_euler('z',
                                vertice_data['angle'],
                                degrees=True).as_matrix()

        rotated_pos = np.dot(BB_T3_A.atom_pos, R_Matrix) + vertice_pos
        final_pos += rotated_pos.tolist()

        vertice_data = topology_info['vertices'][1]
        final_label += BB_T3_B.atom_labels
        vertice_pos = np.array(vertice_data['position'])*a

        R_Matrix = R.from_euler('z',
                                vertice_data['angle'],
                                degrees=True).as_matrix()

        rotated_pos = np.dot(BB_T3_B.atom_pos, R_Matrix) + vertice_pos
        final_pos += rotated_pos.tolist()

        # Replace "X" on final_label with the correct bond atom
        final_label = [x.replace('X', bond_atom) for x in final_label]

        FinalFramework = Structure(lattice,
                                   final_label,
                                   final_pos,
                                   coords_are_cartesian=True)
        FinalFramework.sort(reverse=True)

        # Remove duplicated atoms
        FinalFramework.merge_sites(tol=.5, mode='delete')

        # Translates the structure to the center of the cell
        FinalFramework.translate_sites(
            range(len(FinalFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True)

        # Symmetrize the structure
        try:
            symm = SpacegroupAnalyzer(FinalFramework,
                                      symprec=self.symm_tol,
                                      angle_tolerance=self.angle_tol)
            SymmPrimFramework = symm.get_refined_structure()
        except Exception:
            print('Error in the symmetry analysis. Try to increase the symm_tol or angle_tol parameters')
            return None

        # Create A stacking, a 2D isolated sheet with slab
        if stacking == 'A':
            self.stacking = 'A'
            self.symm_structure = SymmPrimFramework

        # Create AA staking. By default one sheet per unitary cell is used
        if stacking == 'AA':
            self.stacking = 'AA'
            self.symm_structure = SymmPrimFramework

        # Create AB1 staking.
        if stacking == 'AB1':
            self.stacking = 'AB1'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['xyz'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal
            B = ion_conv_crystal

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_1 = Structure(lattice, AB_label + AB_label,
                             AB, coords_are_cartesian=True)

            B_list = np.arange(len(AB_label)) + len(AB_label)

            AB_1.translate_sites(
                B_list, [2/3, 1/3, 0.5], frac_coords=True, to_unit_cell=True)
            AB_1_symm = SpacegroupAnalyzer(AB_1,
                                           symprec=self.symm_tol,
                                           angle_tolerance=self.angle_tol)

            self.symm_structure = AB_1_symm.get_refined_structure()

        # Create AB2 stacking
        if stacking == 'AB2':
            self.stacking = 'AB2'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['xyz'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal
            B = ion_conv_crystal

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_2 = Structure(lattice, AB_label + AB_label,
                             AB, coords_are_cartesian=True)

            B_list = np.arange(len(AB_label)) + len(AB_label)
            AB_2.translate_sites(
                B_list, [1/2, 0, 0.5], frac_coords=True, to_unit_cell=True)
            AB_2_symm = SpacegroupAnalyzer(
                AB_2, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = AB_2_symm.get_refined_structure()

        # Create AAl stacking.

        # Hexagonal cell with two sheets per cell shifited by the shift_vector in angstroms.
        if stacking == 'AAl':
            self.stacking = 'AAl'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 2)

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal * np.array([1, 1, 0.5])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(
                *Tools.cell_to_cellpar(cell))
            shift_vector = np.dot(r, np.array(shift_vector))

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal * np.array([1, 1, 1.5]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            self.symm_structure = Structure(lattice,
                                            AB_label + AB_label,
                                            AB,
                                            coords_are_cartesian=False)

        # Create AA tilted stacking.
        # Tilted Hexagonal cell by tilt_angle with two sheets per cell shifited by the shift_vector
        # in angstroms.
        if stacking == 'AAt':
            self.stacking = 'AAt'

            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['xyz'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = SymmPrimFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            new_cell = Tools.cellpar_to_cell(
                [a_cell, b_cell, c_cell, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + np.array([0, 0, 0.5*c_cell]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(new_cell)
            self.symm_structure = Structure(lattice,
                                            AB_label + AB_label,
                                            AB,
                                            coords_are_cartesian=True)

        if stacking == 'ABC1':
            self.stacking = 'ABC1'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (2/3, 1/3, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (4/3, 2/3, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice,
                              ABC_label + ABC_label + ABC_label,
                              ABC,
                              coords_are_cartesian=False)

            ABC_f_symm = SpacegroupAnalyzer(ABC_f,
                                            symprec=self.symm_tol,
                                            angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_primitive_standard_structure()

        if stacking == 'ABC2':
            self.stacking = 'ABC2'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework['sites']])
            ion_conv_crystal = np.array([i['abc']
                                        for i in SymmPrimFramework['sites']])
            cell = np.array(SymmPrimFramework['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (1/3, 0, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (2/3, 0, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice,
                              ABC_label + ABC_label + ABC_label,
                              ABC,
                              coords_are_cartesian=False)

            ABC_f_symm = SpacegroupAnalyzer(ABC_f,
                                            symprec=self.symm_tol,
                                            angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_primitive_standard_structure()

        dict_structure = self.symm_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

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

    def _create_hcb_a_structure(self,
                                BB_T3: str,
                                BB_L2: str,
                                stacking: str = 'AA',
                                print_result: bool = True,
                                slab: float = 10.0,
                                shift_vector: list = [1.0, 1.0, 0],
                                tilt_angle: float = 5.0):
        """Creates a COF with HCB-A network.

        The HCB-A net is composed of one tripodal and one linear building blocks.

        Parameters
        ----------
        BB_T3 : BuildingBlock, required
            The BuildingBlock object of the tripodal Buiding Block
        BB_L2 : BuildingBlock, required
            The BuildingBlock object of the linear Buiding Block
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
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

        connectivity_error = 'Building block {} must present connectivity {}'
        assert BB_T3.connectivity == 3, connectivity_error.format('A', 3)
        assert BB_L2.connectivity == 2, connectivity_error.format('B', 2)

        self.topology = 'HCB_A'
        self.dimension = 2

        self.charge = BB_L2.charge + BB_T3.charge
        self.chirality = BB_L2.chirality or BB_T3.chirality

        self.name = f'{BB_T3.name}-{BB_L2.name}-HCB_A-{stacking}'

        Tools.print_comand(f'Starting the creation of {self.name}',
                           self.verbosity,
                           ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = Tools.get_bond_atom(BB_T3.conector, BB_L2.conector)

        Tools.print_comand(f'Bond atom detected: {bond_atom}',
                           self.verbosity,
                           ['debug', 'high'])

        # Get the topology information
        topology_info = TOPOLOGY_DICT[self.topology]

        # Measure the base size of the building blocks
        size = BB_T3.size[0] + BB_L2.size[0]

        # Calculate the delta size to add to the c parameter
        delta_a = abs(max(np.transpose(BB_T3.atom_pos)[2])) + abs(min(np.transpose(BB_T3.atom_pos)[2]))
        delta_b = abs(max(np.transpose(BB_L2.atom_pos)[2])) + abs(min(np.transpose(BB_L2.atom_pos)[2]))

        delta_max = max([delta_a, delta_b])

        # Calculate the cell parameters
        a = topology_info['a'] * size
        b = topology_info['b'] * size
        c = topology_info['c'] + delta_max
        alpha = topology_info['alpha']
        beta = topology_info['beta']
        gamma = topology_info['gamma']

        # Create the lattice
        lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        final_label = []
        final_pos = []

        # Add the building blocks to the structure
        for vertice_data in topology_info['vertices']:
            final_label += BB_T3.atom_labels
            vertice_pos = np.array(vertice_data['position'])*a

            R_Matrix = R.from_euler('z', vertice_data['angle'], degrees=True).as_matrix()

            rotated_pos = np.dot(BB_T3.atom_pos, R_Matrix) + vertice_pos
            final_pos += rotated_pos.tolist()

        # Add the building blocks to the structure
        for edge_data in topology_info['edges']:
            final_label += BB_L2.atom_labels
            edge_pos = np.array(edge_data['position'])*a

            R_Matrix = R.from_euler('z', edge_data['angle'], degrees=True).as_matrix()

            rotated_pos = np.dot(BB_L2.atom_pos, R_Matrix) + edge_pos
            final_pos += rotated_pos.tolist()

        # Replace "X" on final_label with the correct bond atom
        final_label = [x.replace('X', bond_atom) for x in final_label]

        FinalFramework = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)
        FinalFramework.sort(reverse=True)

        # Remove duplicated atoms
        FinalFramework.merge_sites(tol=.5, mode='delete')

        # Translates the structure to the center of the cell
        FinalFramework.translate_sites(
            range(len(FinalFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True)

        # Symmetrize the structure
        try:
            symm = SpacegroupAnalyzer(FinalFramework, symprec=0.3, angle_tolerance=3.0)
            SymmPrimFramework = symm.get_refined_structure()
        except Exception:
            print('Error in the symmetry analysis. Try to increase the symm_tol or angle_tol parameters')
            return None

        # Create A stacking. The slab is defined by the c_cell parameter
        if stacking == 'A':
            self.stacking = 'A'
            self.symm_structure = SymmPrimFramework

        # Create AA staking. By default one sheet per unitary cell is used
        if stacking == 'AA':
            self.stacking = 'AA'
            self.symm_structure = SymmPrimFramework

        if stacking == 'AB1':
            self.stacking = 'AB1'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']]
                )
            ion_conv_crystal = np.array(
                [i['xyz'] for i in SymmPrimFramework.as_dict()['sites']]
                )
            cell = np.array(SymmPrimFramework.as_dict()['lattice']['matrix'])*(1, 1, 2)

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
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['xyz'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal
            B = ion_conv_crystal

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_2 = Structure(lattice, AB_label + AB_label,
                             AB, coords_are_cartesian=True)

            B_list = np.arange(len(AB_label)) + len(AB_label)
            AB_2.translate_sites(
                B_list, [1/2, 0, 0.5], frac_coords=True, to_unit_cell=True)
            AB_2_symm = SpacegroupAnalyzer(
                AB_2, symprec=0.05, angle_tolerance=0.5)

            self.symm_structure = AB_2_symm.get_refined_structure()

        # Create AAl stacking.
        # Hexagonal cell with two sheets per cell shifited by the shift_vector in angstroms.
        if stacking == 'AAl':
            self.stacking = 'AAl'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 2)

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal * np.array([1, 1, 0.5])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(
                *Tools.cell_to_cellpar(cell))
            shift_vector = np.dot(r, np.array(shift_vector))

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal * np.array([1, 1, 1.5]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            self.symm_structure = Structure(lattice,
                                            AB_label + AB_label,
                                            AB,
                                            coords_are_cartesian=False)

        # Create AA tilted stacking.
        # Tilted Hexagonal cell with two sheets per cell tilted by tilt_angle.
        if stacking == 'AAt':
            self.stacking = 'AAt'

            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['xyz'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = SymmPrimFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            new_cell = Tools.cellpar_to_cell(
                [a_cell, b_cell, c_cell, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + np.array([0, 0, 0.5*c_cell]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(new_cell)
            self.symm_structure = Structure(lattice,
                                            AB_label + AB_label,
                                            AB,
                                            coords_are_cartesian=True)

        if stacking == 'ABC1':
            self.stacking = 'ABC1'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (2/3, 1/3, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (4/3, 2/3, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice,
                              ABC_label + ABC_label + ABC_label,
                              ABC,
                              coords_are_cartesian=False)

            ABC_f_symm = SpacegroupAnalyzer(
                ABC_f, symprec=0.2, angle_tolerance=2.0)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        if stacking == 'ABC2':
            self.stacking = 'ABC2'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (1/3, 0, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (2/3, 0, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice,
                              ABC_label + ABC_label + ABC_label,
                              ABC,
                              coords_are_cartesian=False)

            ABC_f_symm = SpacegroupAnalyzer(
                ABC_f, symprec=0.2, angle_tolerance=2.0)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        dict_structure = self.symm_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

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

    def _create_sql_structure(self,
                              BB_S4_A: str,
                              BB_S4_B: str,
                              stacking: str = 'AA',
                              print_result: bool = True,
                              slab: float = 10.0,
                              shift_vector: list = [1.0, 1.0, 0],
                              tilt_angle: float = 5.0):
        """Creates a COF with SQL network.

        The SQL net is composed of two tetrapodal building blocks.

        Parameters
        ----------
        BB_S4_A : BuildingBlock, required
            The BuildingBlock object of the tetrapodal Buiding Block A
        BB_S4_B : BuildingBlock, required
            The BuildingBlock object of the tetrapodal Buiding Block B
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
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

        connectivity_error = 'Building block {} must present connectivity {}'
        assert BB_S4_A.connectivity == 4, connectivity_error.format('A', 4)
        assert BB_S4_B.connectivity == 4, connectivity_error.format('B', 4)

        self.topology = 'SQL'
        self.dimension = 2

        self.charge = BB_S4_A.charge + BB_S4_B.charge
        self.chirality = BB_S4_A.chirality or BB_S4_B.chirality

        self.name = f'{BB_S4_A.name}-{BB_S4_B.name}-{self.topology}-{stacking}'

        Tools.print_comand(f'Starting the creation of {self.name}',
                           self.verbosity,
                           ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = Tools.get_bond_atom(BB_S4_A.conector, BB_S4_B.conector)

        Tools.print_comand(f'Bond atom detected: {bond_atom}',
                           self.verbosity,
                           ['debug', 'high'])

        # Get the topology information
        topology_info = TOPOLOGY_DICT[self.topology]

        # Measure the base size of the building blocks
        size = BB_S4_A.size[0] + BB_S4_B.size[0]

        # Calculate the delta size to add to the c parameter
        delta_a = abs(max(np.transpose(BB_S4_A.atom_pos)[2])) + abs(min(np.transpose(BB_S4_B.atom_pos)[2]))
        delta_b = abs(max(np.transpose(BB_S4_A.atom_pos)[2])) + abs(min(np.transpose(BB_S4_B.atom_pos)[2]))

        delta_max = max([delta_a, delta_b])

        # Calculate the cell parameters
        a = topology_info['a'] * size
        b = topology_info['b'] * size
        c = topology_info['c'] + delta_max
        alpha = topology_info['alpha']
        beta = topology_info['beta']
        gamma = topology_info['gamma']

        # Create the lattice
        lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        final_label = []
        final_pos = []

        # Add the first building block to the structure
        vertice_data = topology_info['vertices'][0]
        final_label += BB_S4_A.atom_labels
        vertice_pos = np.array(vertice_data['position'])*a

        R_Matrix = R.from_euler('z', vertice_data['angle'], degrees=True).as_matrix()

        rotated_pos = np.dot(BB_S4_A.atom_pos, R_Matrix) + vertice_pos
        final_pos += rotated_pos.tolist()

        # Add the second building block to the structure
        vertice_data = topology_info['vertices'][1]
        final_label += BB_S4_B.atom_labels
        vertice_pos = np.array(vertice_data['position'])*a

        R_Matrix = R.from_euler('z', vertice_data['angle'], degrees=True).as_matrix()

        rotated_pos = np.dot(BB_S4_B.atom_pos, R_Matrix) + vertice_pos
        final_pos += rotated_pos.tolist()

        # Replace "X" on final_label with the correct bond atom
        final_label = [x.replace('X', bond_atom) for x in final_label]

        FinalFramework = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)
        FinalFramework.sort(reverse=True)

        # Remove duplicated atoms
        FinalFramework.merge_sites(tol=.5, mode='delete')

        # Translates the structure to the center of the cell
        FinalFramework.translate_sites(
            range(len(FinalFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True)

        # Symmetrize the structure
        try:
            symm = SpacegroupAnalyzer(FinalFramework, symprec=0.3, angle_tolerance=3.0)
            SymmPrimFramework = symm.get_refined_structure()
        except Exception:
            print('Error in the symmetry analysis. Try to increase the symm_tol or angle_tol parameters')
            return None

        if stacking not in self.available_stacking[self.topology]:
            raise Exception(f"""{stacking} is not in the available stack list for HCB net.
    Available options are: {self.available_stacking[self.topology]}""")

        # Create A stacking, a 2D isolated sheet with slab
        if stacking == 'A':
            self.stacking = 'A'
            self.symm_structure = SymmPrimFramework

        if stacking == 'AA':
            self.stacking = 'AA'
            self.symm_structure = SymmPrimFramework

        if stacking == 'AB1':
            self.stacking = 'AB1'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal*(1, 1, 0.5)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1.5) + (1/4, 1/4, 0))

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_1 = Structure(lattice, AB_label + AB_label,
                             AB, coords_are_cartesian=False)
            AB_1_symm = SpacegroupAnalyzer(
                AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = AB_1_symm.get_refined_structure()

        if stacking == 'AB2':
            self.stacking = 'AB2'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal*(1, 1, 0.5)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1.5) + (1/2, 0, 0))

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
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 2)

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal * np.array([1, 1, 0.5])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(
                *Tools.cell_to_cellpar(cell))
            shift_vector = np.dot(r, np.array(shift_vector))

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal * np.array([1, 1, 1.5]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            self.symm_structure = Structure(lattice,
                                            AB_label + AB_label,
                                            AB,
                                            coords_are_cartesian=False)

        # Create AA tilted stacking.
        # Tilted tetragonal cell with two sheets per cell tilted by tilt_angle.
        if stacking == 'AAt':
            self.stacking = 'AAt'

            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['xyz'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = SymmPrimFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            new_cell = Tools.cellpar_to_cell(
                [a_cell, b_cell, c_cell, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + np.array([0, 0, 0.5*c_cell]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(new_cell)
            self.symm_structure = Structure(lattice,
                                            AB_label + AB_label,
                                            AB,
                                            coords_are_cartesian=True)

        if stacking == 'ABC1':
            self.stacking = 'ABC1'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (2/3, 1/3, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (4/3, 2/3, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice,
                              ABC_label + ABC_label + ABC_label,
                              ABC,
                              coords_are_cartesian=False)

            ABC_f_symm = SpacegroupAnalyzer(ABC_f,
                                            symprec=self.symm_tol,
                                            angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        if stacking == 'ABC2':
            self.stacking = 'ABC2'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (1/3, 0, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (2/3, 0, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice,
                              ABC_label + ABC_label + ABC_label,
                              ABC,
                              coords_are_cartesian=False)

            ABC_f_symm = SpacegroupAnalyzer(ABC_f,
                                            symprec=self.symm_tol,
                                            angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        dict_structure = self.symm_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

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

    def _create_sql_a_structure(self,
                                BB_S4: str,
                                BB_L2: str,
                                stacking: str = 'AA',
                                c_parameter_base: float = 3.6,
                                print_result: bool = True,
                                slab: float = 10.0,
                                shift_vector: list = [1.0, 1.0, 0],
                                tilt_angle: float = 5.0):
        """Creates a COF with SQL-A network.

        The SQL-A net is composed of one tetrapodal and one linear building blocks.

        Parameters
        ----------
        BB_S4 : BuildingBlock, required
            The BuildingBlock object of the tetrapodal Buiding Block
        BB_L2 : BuildingBlock, required
            The BuildingBlock object of the bipodal Buiding Block
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
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

        connectivity_error = 'Building block {} must present connectivity {}'
        assert BB_S4.connectivity == 4, connectivity_error.format('A', 4)
        assert BB_L2.connectivity == 2, connectivity_error.format('B', 2)

        self.topology = 'SQL_A'
        self.dimension = 2

        self.charge = BB_S4.charge + BB_L2.charge
        self.chirality = BB_S4.chirality or BB_L2.chirality

        self.name = f'{BB_S4.name}-{BB_L2.name}-{self.topology}-{stacking}'

        Tools.print_comand(f'Starting the creation of {self.name}',
                           self.verbosity,
                           ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = Tools.get_bond_atom(BB_S4.conector, BB_L2.conector)

        Tools.print_comand(f'Bond atom detected: {bond_atom}',
                           self.verbosity,
                           ['debug', 'high'])

        # Get the topology information
        topology_info = TOPOLOGY_DICT[self.topology]

        # Measure the base size of the building blocks
        size = BB_S4.size[0] + BB_L2.size[0]

        # Calculate the delta size to add to the c parameter
        delta_a = abs(max(np.transpose(BB_S4.atom_pos)[2])) + abs(min(np.transpose(BB_S4.atom_pos)[2]))
        delta_b = abs(max(np.transpose(BB_L2.atom_pos)[2])) + abs(min(np.transpose(BB_L2.atom_pos)[2]))

        delta_max = max([delta_a, delta_b])

        # Calculate the cell parameters
        a = topology_info['a'] * size
        b = topology_info['b'] * size
        c = topology_info['c'] + delta_max
        alpha = topology_info['alpha']
        beta = topology_info['beta']
        gamma = topology_info['gamma']

        # Create the lattice
        lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        final_label = []
        final_pos = []

        # Add the building blocks to the structure
        for vertice_data in topology_info['vertices']:
            final_label += BB_S4.atom_labels
            vertice_pos = np.array(vertice_data['position'])*a

            R_Matrix = R.from_euler('z', vertice_data['angle'], degrees=True).as_matrix()

            rotated_pos = np.dot(BB_S4.atom_pos, R_Matrix) + vertice_pos
            final_pos += rotated_pos.tolist()

        # Add the building blocks to the structure
        for edge_data in topology_info['edges']:
            final_label += BB_L2.atom_labels
            edge_pos = np.array(edge_data['position'])*a

            R_Matrix = R.from_euler('z', edge_data['angle'], degrees=True).as_matrix()

            rotated_pos = np.dot(BB_L2.atom_pos, R_Matrix) + edge_pos
            final_pos += rotated_pos.tolist()

        # Replace "X" on final_label with the correct bond atom
        bond_atom = Tools.get_bond_atom(BB_S4.conector, BB_L2.conector)

        final_label = [x.replace('X', bond_atom) for x in final_label]

        FinalFramework = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)
        FinalFramework.sort(reverse=True)

        # Remove duplicated atoms
        FinalFramework.merge_sites(tol=.5, mode='delete')

        # Translates the structure to the center of the cell
        FinalFramework.translate_sites(
            range(len(FinalFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True)

        # Symmetrize the structure
        try:
            symm = SpacegroupAnalyzer(FinalFramework, symprec=0.3, angle_tolerance=3.0)
            SymmPrimFramework = symm.get_refined_structure()
        except Exception:
            print('Error in the symmetry analysis. Try to increase the symm_tol or angle_tol parameters')
            return None

        if stacking not in self.available_stacking[self.topology]:
            raise Exception(f"""{stacking} is not in the available stack list for HCB net.
    Available options are: {self.available_stacking[self.topology]}""")

        # Create A stacking, a 2D isolated sheet with slab
        if stacking == 'A':
            self.symm_structure = SymmPrimFramework

        if stacking == 'AA':
            self.symm_structure = SymmPrimFramework

        if stacking == 'AB1':
            labels_conv_crystal = [[i['label']] for i in SymmPrimFramework.as_dict()['sites']]
            labels_conv_crystal = np.array(labels_conv_crystal)

            ion_conv_crystal = [i['abc'] for i in SymmPrimFramework.as_dict()['sites']]
            ion_conv_crystal = np.array(ion_conv_crystal).astype(float)

            cell = np.array(SymmPrimFramework.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal*(1, 1, 0.5)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1.5) + (1/4, 1/4, 0))

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_1 = Structure(lattice, AB_label + AB_label, AB, coords_are_cartesian=False)
            AB_1_symm = SpacegroupAnalyzer(AB_1,
                                           symprec=self.symm_tol,
                                           angle_tolerance=self.angle_tol)

            self.symm_structure = AB_1_symm.get_refined_structure()

        if stacking == 'AB2':
            labels_conv_crystal = [[i['label']] for i in SymmPrimFramework.as_dict()['sites']]
            labels_conv_crystal = np.array(labels_conv_crystal)

            ion_conv_crystal = [i['abc'] for i in SymmPrimFramework.as_dict()['sites']]
            ion_conv_crystal = np.array(ion_conv_crystal).astype(float)

            cell = np.array(SymmPrimFramework.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal*(1, 1, 0.5)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1.5) + (1/2, 0, 0))

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_2 = Structure(lattice, AB_label + AB_label, AB, coords_are_cartesian=False)
            AB_2_symm = SpacegroupAnalyzer(AB_2,
                                           symprec=self.symm_tol,
                                           angle_tolerance=self.angle_tol)

            self.symm_structure = AB_2_symm.get_refined_structure()

        # Create AAl stacking.
        # Tetragonal cell with two sheets per cell shifited by the shift_vector in angstroms.
        if stacking == 'AAl':
            labels_conv_crystal = [[i['label']] for i in SymmPrimFramework.as_dict()['sites']]
            labels_conv_crystal = np.array(labels_conv_crystal)

            ion_conv_crystal = [i['abc'] for i in SymmPrimFramework.as_dict()['sites']]
            ion_conv_crystal = np.array(ion_conv_crystal).astype(float)

            cell = np.array(SymmPrimFramework.as_dict()['lattice']['matrix'])*(1, 1, 2)

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal * np.array([1, 1, 0.5])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(
                *Tools.cell_to_cellpar(cell))
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
            labels_conv_crystal = [[i['label']] for i in SymmPrimFramework.as_dict()['sites']]
            labels_conv_crystal = np.array(labels_conv_crystal)

            ion_conv_crystal = [i['abc'] for i in SymmPrimFramework.as_dict()['sites']]
            ion_conv_crystal = np.array(ion_conv_crystal).astype(float)

            cell = SymmPrimFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_parameter_base = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            new_cell = Tools.cellpar_to_cell(
                [a_cell, b_cell, c_parameter_base, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + \
                np.array([0, 0, 0.5*c_parameter_base]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(new_cell)
            self.symm_structure = Structure(
                lattice, AB_label+AB_label, AB, coords_are_cartesian=True)

        if stacking == 'ABC1':
            labels_conv_crystal = [[i['label']] for i in SymmPrimFramework.as_dict()['sites']]
            labels_conv_crystal = np.array(labels_conv_crystal)

            ion_conv_crystal = [i['abc'] for i in SymmPrimFramework.as_dict()['sites']]
            ion_conv_crystal = np.array(ion_conv_crystal).astype(float)

            cell = np.array(SymmPrimFramework.as_dict()['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (2/3, 1/3, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (4/3, 2/3, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice,
                              ABC_label + ABC_label + ABC_label,
                              ABC,
                              coords_are_cartesian=False)

            ABC_f_symm = SpacegroupAnalyzer(ABC_f,
                                            symprec=self.symm_tol,
                                            angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        if stacking == 'ABC2':
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (1/3, 0, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (2/3, 0, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label +
                              ABC_label, ABC, coords_are_cartesian=False)
            ABC_f_symm = SpacegroupAnalyzer(
                ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        dict_structure = self.symm_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

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

    def _create_kgd_structure(self,
                              BB_H6: str,
                              BB_T3: str,
                              stacking: str = 'AA',
                              print_result: bool = True,
                              slab: float = 10.0,
                              shift_vector: list = [1.0, 1.0, 0],
                              tilt_angle: float = 5.0):
        """Creates a COF with KGD network.

        The KGD net is composed of one hexapodal and one tripodal building blocks.

        Parameters
        ----------
        BB_H6 : BuildingBlock, required
            The BuildingBlock object of the hexapodal Buiding Block
        BB_T3 : BuildingBlock, required
            The BuildingBlock object of the tripodal Buiding Block
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
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
        connectivity_error = 'Building block {} must present connectivity {}'
        assert BB_H6.connectivity == 6, connectivity_error.format('A', 6)
        assert BB_T3.connectivity == 3, connectivity_error.format('B', 3)

        self.topology = 'KGD'
        self.dimension = 2

        self.charge = BB_H6.charge + BB_T3.charge
        self.chirality = BB_H6.chirality or BB_T3.chirality

        self.name = f'{BB_H6.name}-{BB_T3.name}-{self.topology}-{stacking}'

        Tools.print_comand(f'Starting the creation of {self.name}',
                           self.verbosity,
                           ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = Tools.get_bond_atom(BB_H6.conector, BB_T3.conector)

        Tools.print_comand(f'Bond atom detected: {bond_atom}',
                           self.verbosity,
                           ['debug', 'high'])

        # Get the topology information
        topology_info = TOPOLOGY_DICT[self.topology]

        # Measure the base size of the building blocks
        size = BB_H6.size[0] + BB_T3.size[0]

        # Calculate the delta size to add to the c parameter
        delta_a = abs(max(np.transpose(BB_H6.atom_pos)[2])) + abs(min(np.transpose(BB_H6.atom_pos)[2]))
        delta_b = abs(max(np.transpose(BB_T3.atom_pos)[2])) + abs(min(np.transpose(BB_T3.atom_pos)[2]))

        delta_max = max([delta_a, delta_b])

        # Calculate the cell parameters
        a = topology_info['a'] * size
        b = topology_info['b'] * size
        c = topology_info['c'] + delta_max
        alpha = topology_info['alpha']
        beta = topology_info['beta']
        gamma = topology_info['gamma']

        # Create the lattice
        lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        final_label = []
        final_pos = []

        # Add the building blocks to the structure
        for vertice_data in topology_info['vertices']:
            final_label += BB_H6.atom_labels
            vertice_pos = np.array(vertice_data['position'])*a

            R_Matrix = R.from_euler('z',
                                    vertice_data['angle'],
                                    degrees=True).as_matrix()

            rotated_pos = np.dot(BB_H6.atom_pos, R_Matrix) + vertice_pos
            final_pos += rotated_pos.tolist()

        # Add the building blocks to the structure
        for edge_data in topology_info['edges']:
            final_label += BB_T3.atom_labels
            edge_pos = np.array(edge_data['position'])*a

            R_Matrix = R.from_euler('z',
                                    edge_data['angle'],
                                    degrees=True).as_matrix()

            rotated_pos = np.dot(BB_T3.atom_pos, R_Matrix) + edge_pos
            final_pos += rotated_pos.tolist()

        # Replace "X" on final_label with the correct bond atom
        bond_atom = Tools.get_bond_atom(BB_H6.conector, BB_T3.conector)

        final_label = [x.replace('X', bond_atom) for x in final_label]

        FinalFramework = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)
        FinalFramework.sort(reverse=True)

        # Remove duplicated atoms
        FinalFramework.merge_sites(tol=.5, mode='delete')

        # Translates the structure to the center of the cell
        FinalFramework.translate_sites(
            range(len(FinalFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True)

        # Symmetrize the structure
        try:
            symm = SpacegroupAnalyzer(FinalFramework, symprec=0.3, angle_tolerance=3.0)
            SymmPrimFramework = symm.get_refined_structure()
        except Exception:
            print('Error in the symmetry analysis. Try to increase the symm_tol or angle_tol parameters')
            return None

        if stacking not in self.available_stacking[self.topology]:
            raise Exception(f"""{stacking} is not in the available stack list for HCB net.
    Available options are: {self.available_stacking[self.topology]}""")

        # Create A stacking, a 2D isolated sheet with slab
        if stacking == 'A':
            self.stacking = 'A'
            self.symm_structure = SymmPrimFramework

        if stacking == 'AA':
            self.stacking = 'AA'
            self.symm_structure = SymmPrimFramework

        # Create AB1 staking.
        if stacking == 'AB1':
            self.stacking = 'AB1'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz']
                                        for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal
            B = ion_conv_crystal

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_1 = Structure(lattice, AB_label + AB_label,
                             AB, coords_are_cartesian=True)

            B_list = np.arange(len(AB_label)) + len(AB_label)

            AB_1.translate_sites(
                B_list, [2/3, 1/3, 0.5], frac_coords=True, to_unit_cell=True)
            symm = SpacegroupAnalyzer(
                AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        # Create AB2 stacking
        if stacking == 'AB2':
            self.stacking = 'AB2'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz']
                                        for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal
            B = ion_conv_crystal

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_2 = Structure(lattice, AB_label + AB_label,
                             AB, coords_are_cartesian=True)

            B_list = np.arange(len(AB_label)) + len(AB_label)
            AB_2.translate_sites(
                B_list, [1/2, 0, 0.5], frac_coords=True, to_unit_cell=True)
            symm = SpacegroupAnalyzer(
                AB_2, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        # Create AAl stacking.
        # Hexagonal cell with two sheets per cell shifited by the shift_vector in angstroms.
        if stacking == 'AAl':
            self.stacking = 'AAl'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc']
                                        for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()['lattice']['matrix'])*(1, 1, 2)

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal * np.array([1, 1, 0.5])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(
                *Tools.cell_to_cellpar(cell))
            shift_vector = np.dot(r, np.array(shift_vector))

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal * np.array([1, 1, 1.5]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            symm = Structure(lattice, AB_label+AB_label,
                             AB, coords_are_cartesian=False)

            self.symm_structure = symm.get_refined_structure()

        # Create AA tilted stacking.
        # Tilted Hexagonal cell by tilt_angle with two sheets per cell shifited by
        # the shift_vector in angstroms.
        if stacking == 'AAt':
            self.stacking = 'AAt'

            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz']
                                        for i in SymmPrimFramework.as_dict()['sites']])
            cell = SymmPrimFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            new_cell = Tools.cellpar_to_cell(
                [a_cell, b_cell, c_cell, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + np.array([0, 0, 0.5*c_cell]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(new_cell)
            symm = Structure(lattice, AB_label+AB_label,
                             AB, coords_are_cartesian=True)

            self.symm_structure = symm.get_refined_structure()

        if stacking == 'ABC1':
            self.stacking = 'ABC1'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc']
                                        for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (2/3, 1/3, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (4/3, 2/3, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label +
                              ABC_label, ABC, coords_are_cartesian=False)
            symm = SpacegroupAnalyzer(
                ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_primitive_standard_structure()

        if stacking == 'ABC2':
            self.stacking = 'ABC2'
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in SymmPrimFramework['sites']])
            cell = np.array(SymmPrimFramework['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (1/3, 0, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (2/3, 0, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label +
                              ABC_label, ABC, coords_are_cartesian=False)
            symm = SpacegroupAnalyzer(
                ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_primitive_standard_structure()

        dict_structure = self.symm_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

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

    def _create_hxl_a_structure(self,
                                BB_H6: str,
                                BB_L2: str,
                                stacking: str = 'AA',
                                print_result: bool = True,
                                slab: float = 10.0,
                                shift_vector: list = [1.0, 1.0, 0],
                                tilt_angle: float = 5.0):
        """Creates a COF with HXL-A network.

        The HXK-A net is composed of one hexapodal and one linear building blocks.

        Parameters
        ----------
        BB_H6 : BuildingBlock, required
            The BuildingBlock object of the tetrapodal Buiding Block A
        BB_L2 : BuildingBlock, required
            The BuildingBlock object of the tetrapodal Buiding Block B
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
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

        connectivity_error = 'Building block {} must present connectivity {}'
        assert BB_H6.connectivity == 6, connectivity_error.format('A', 6)
        assert BB_L2.connectivity == 2, connectivity_error.format('B', 2)

        self.topology = 'HXL_A'
        self.dimension = 2

        self.charge = BB_H6.charge + BB_L2.charge
        self.chirality = BB_H6.chirality or BB_L2.chirality

        self.name = f'{BB_H6.name}-{BB_L2.name}-{self.topology}-{stacking}'

        Tools.print_comand(f'Starting the creation of {self.name}',
                           self.verbosity,
                           ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = Tools.get_bond_atom(BB_H6.conector, BB_L2.conector)

        Tools.print_comand(f'Bond atom detected: {bond_atom}',
                           self.verbosity,
                           ['debug', 'high'])

        # Get the topology information
        topology_info = TOPOLOGY_DICT[self.topology]

        # Measure the base size of the building blocks
        size = BB_H6.size[0] + BB_L2.size[0]

        # Calculate the delta size to add to the c parameter
        delta_a = abs(max(np.transpose(BB_H6.atom_pos)[2])) + abs(min(np.transpose(BB_H6.atom_pos)[2]))
        delta_b = abs(max(np.transpose(BB_L2.atom_pos)[2])) + abs(min(np.transpose(BB_L2.atom_pos)[2]))

        delta_max = max([delta_a, delta_b])

        # Calculate the cell parameters
        a = topology_info['a'] * size
        b = topology_info['b'] * size
        c = topology_info['c'] + delta_max
        alpha = topology_info['alpha']
        beta = topology_info['beta']
        gamma = topology_info['gamma']

        # Create the lattice
        lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        final_label = []
        final_pos = []

        # Add the building blocks to the structure
        for vertice_data in topology_info['vertices']:
            final_label += BB_H6.atom_labels
            vertice_pos = np.array(vertice_data['position'])*a

            R_Matrix = R.from_euler('z',
                                    vertice_data['angle'],
                                    degrees=True).as_matrix()

            rotated_pos = np.dot(BB_H6.atom_pos, R_Matrix) + vertice_pos
            final_pos += rotated_pos.tolist()

        # Add the building blocks to the structure
        for edge_data in topology_info['edges']:
            final_label += BB_L2.atom_labels
            edge_pos = np.array(edge_data['position'])*a

            R_Matrix = R.from_euler('z',
                                    edge_data['angle'],
                                    degrees=True).as_matrix()

            rotated_pos = np.dot(BB_L2.atom_pos, R_Matrix) + edge_pos
            final_pos += rotated_pos.tolist()

        # Replace "X" on final_label with the correct bond atom
        bond_atom = Tools.get_bond_atom(BB_H6.conector, BB_L2.conector)

        final_label = [x.replace('X', bond_atom) for x in final_label]

        FinalFramework = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)
        FinalFramework.sort(reverse=True)

        # Remove duplicated atoms
        FinalFramework.merge_sites(tol=.5, mode='delete')

        # Translates the structure to the center of the cell
        FinalFramework.translate_sites(
            range(len(FinalFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True)

        # Symmetrize the structure
        try:
            symm = SpacegroupAnalyzer(FinalFramework, symprec=0.3, angle_tolerance=3.0)
            SymmPrimFramework = symm.get_refined_structure()
        except Exception:
            print('Error in the symmetry analysis. Try to increase the symm_tol or angle_tol parameters')
            return None

        if stacking not in self.available_stacking[self.topology]:
            raise Exception(f"""{stacking} is not in the available stack list for HCB net.
    Available options are: {self.available_stacking[self.topology]}""")

        # Create A stacking, a 2D isolated sheet with slab
        if stacking == 'A':
            self.symm_structure = SymmPrimFramework

        if stacking == 'AA':
            self.symm_structure = SymmPrimFramework

        # Create AB1 staking.
        if stacking == 'AB1':
            self.stacking = 'AB1'
            labels_conv_crystal = np.array(
                [[i['label']] for i in FinalFramework.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz']
                                        for i in FinalFramework.as_dict()['sites']])
            cell = np.array(FinalFramework.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal
            B = ion_conv_crystal

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_1 = Structure(lattice, AB_label + AB_label,
                             AB, coords_are_cartesian=True)

            B_list = np.arange(len(AB_label)) + len(AB_label)

            AB_1.translate_sites(
                B_list, [2/3, 1/3, 0.5], frac_coords=True, to_unit_cell=True)
            symm = SpacegroupAnalyzer(
                AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        # Create AB2 stacking
        if stacking == 'AB2':
            self.stacking = 'AB2'
            labels_conv_crystal = np.array(
                [[i['label']] for i in FinalFramework.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz']
                                        for i in FinalFramework.as_dict()['sites']])
            cell = np.array(FinalFramework.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal
            B = ion_conv_crystal

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_2 = Structure(lattice, AB_label + AB_label,
                             AB, coords_are_cartesian=True)

            B_list = np.arange(len(AB_label)) + len(AB_label)
            AB_2.translate_sites(
                B_list, [1/2, 0, 0.5], frac_coords=True, to_unit_cell=True)
            symm = SpacegroupAnalyzer(
                AB_2, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        # Create AAl stacking.
        # Hexagonal cell with two sheets per cell shifited by the shift_vector in angstroms.
        if stacking == 'AAl':
            self.stacking = 'AAl'
            labels_conv_crystal = np.array(
                [[i['label']] for i in FinalFramework.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc']
                                        for i in FinalFramework.as_dict()['sites']])
            cell = np.array(FinalFramework.as_dict()['lattice']['matrix'])*(1, 1, 2)

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal * np.array([1, 1, 0.5])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(
                *Tools.cell_to_cellpar(cell))
            shift_vector = np.dot(r, np.array(shift_vector))

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal * np.array([1, 1, 1.5]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            symm = Structure(lattice, AB_label+AB_label,
                             AB, coords_are_cartesian=False)

            self.symm_structure = symm.get_refined_structure()

        # Create AA tilted stacking.
        # Tilted Hexagonal cell by tilt_angle with two sheets per cell
        # shifited by the shift_vector in angstroms.
        if stacking == 'AAt':
            self.stacking = 'AAt'

            labels_conv_crystal = np.array(
                [[i['label']] for i in FinalFramework.as_dict()['sites']])
            ion_conv_crystal = np.array([i['xyz']
                                        for i in FinalFramework.as_dict()['sites']])
            cell = FinalFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            new_cell = Tools.cellpar_to_cell(
                [a_cell, b_cell, c_cell, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + np.array([0, 0, 0.5*c_cell]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(new_cell)
            symm = Structure(lattice, AB_label+AB_label,
                             AB, coords_are_cartesian=True)

            self.symm_structure = symm.get_refined_structure()

        if stacking == 'ABC1':
            self.stacking = 'ABC1'
            labels_conv_crystal = np.array(
                [[i['label']] for i in FinalFramework.as_dict()['sites']])
            ion_conv_crystal = np.array([i['abc']
                                        for i in FinalFramework.as_dict()['sites']])
            cell = np.array(FinalFramework.as_dict()['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (2/3, 1/3, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (4/3, 2/3, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label +
                              ABC_label, ABC, coords_are_cartesian=False)
            symm = SpacegroupAnalyzer(
                ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_primitive_standard_structure()

        if stacking == 'ABC2':
            self.stacking = 'ABC2'
            labels_conv_crystal = np.array(
                [[i['label']] for i in FinalFramework['sites']])
            ion_conv_crystal = np.array([i['abc'] for i in FinalFramework['sites']])
            cell = np.array(FinalFramework['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (1/3, 0, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (2/3, 0, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label +
                              ABC_label, ABC, coords_are_cartesian=False)
            symm = SpacegroupAnalyzer(
                ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_primitive_standard_structure()

        dict_structure = self.symm_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

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

    def _create_fxt_a_structure(self,
                                BB_S4: str,
                                BB_L2: str,
                                stacking: str = 'AA',
                                c_parameter_base: float = 3.6,
                                print_result: bool = True,
                                slab: float = 10.0,
                                shift_vector: list = [1.0, 1.0, 0],
                                tilt_angle: float = 5.0):
        """Creates a COF with FXT-A network.

        The SQL-A net is composed of one tetrapodal and one linear building blocks.

        Parameters
        ----------
        BB_S4 : BuildingBlock, required
            The BuildingBlock object of the tetrapodal Buiding Block
        BB_L2 : BuildingBlock, required
            The BuildingBlock object of the bipodal Buiding Block
        stacking : str, optional
            The stacking pattern of the COF layers (default is 'AA')
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

        connectivity_error = 'Building block {} must present connectivity {}'
        assert BB_S4.connectivity == 4, connectivity_error.format('A', 4)
        assert BB_L2.connectivity == 2, connectivity_error.format('B', 2)

        self.topology = 'FXT_A'
        self.dimension = 2

        self.charge = BB_S4.charge + BB_L2.charge
        self.chirality = BB_S4.chirality or BB_L2.chirality

        self.name = f'{BB_S4.name}-{BB_L2.name}-{self.topology}-{stacking}'

        Tools.print_comand(f'Starting the creation of {self.name}',
                           self.verbosity,
                           ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = Tools.get_bond_atom(BB_S4.conector, BB_L2.conector)

        Tools.print_comand(f'Bond atom detected: {bond_atom}',
                           self.verbosity,
                           ['debug', 'high'])

        # Get the topology information
        topology_info = TOPOLOGY_DICT[self.topology]

        # Measure the base size of the building blocks
        size = BB_S4.size[0] + BB_L2.size[0]

        # Calculate the delta size to add to the c parameter
        delta_a = abs(max(np.transpose(BB_S4.atom_pos)[2])) + abs(min(np.transpose(BB_S4.atom_pos)[2]))
        delta_b = abs(max(np.transpose(BB_L2.atom_pos)[2])) + abs(min(np.transpose(BB_L2.atom_pos)[2]))

        delta_max = max([delta_a, delta_b])

        # Calculate the cell parameters
        a = topology_info['a'] * size
        b = topology_info['b'] * size
        c = topology_info['c'] + delta_max
        alpha = topology_info['alpha']
        beta = topology_info['beta']
        gamma = topology_info['gamma']

        # Create the lattice
        lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        final_label = []
        final_pos = []

        # Add the building blocks to the structure
        for vertice_data in topology_info['vertices']:
            final_label += BB_S4.atom_labels
            vertice_pos = np.array(vertice_data['position'])*a

            R_Matrix = R.from_euler('z', vertice_data['angle'], degrees=True).as_matrix()

            rotated_pos = np.dot(BB_S4.atom_pos, R_Matrix) + vertice_pos
            final_pos += rotated_pos.tolist()

        # Add the building blocks to the structure
        for edge_data in topology_info['edges']:
            final_label += BB_L2.atom_labels
            edge_pos = np.array(edge_data['position'])*a

            R_Matrix = R.from_euler('z', edge_data['angle'], degrees=True).as_matrix()

            rotated_pos = np.dot(BB_L2.atom_pos, R_Matrix) + edge_pos
            final_pos += rotated_pos.tolist()

        # Replace "X" on final_label with the correct bond atom
        bond_atom = Tools.get_bond_atom(BB_S4.conector, BB_L2.conector)

        final_label = [x.replace('X', bond_atom) for x in final_label]

        FinalFramework = Structure(lattice, final_label, final_pos, coords_are_cartesian=True)
        FinalFramework.sort(reverse=True)

        # Remove duplicated atoms
        FinalFramework.merge_sites(tol=.5, mode='delete')

        # Translates the structure to the center of the cell
        FinalFramework.translate_sites(
            range(len(FinalFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True)

        # Symmetrize the structure
        try:
            symm = SpacegroupAnalyzer(FinalFramework, symprec=0.3, angle_tolerance=3.0)
            SymmPrimFramework = symm.get_refined_structure()
        except Exception:
            print('Error in the symmetry analysis. Try to increase the symm_tol or angle_tol parameters')
            return None

        if stacking not in self.available_stacking[self.topology]:
            raise Exception(f"""{stacking} is not in the available stack list for HCB net.
    Available options are: {self.available_stacking[self.topology]}""")

        # Create A stacking, a 2D isolated sheet with slab
        if stacking == 'A':
            self.symm_structure = SymmPrimFramework

        if stacking == 'AA':
            self.symm_structure = SymmPrimFramework

        if stacking == 'AB1':
            labels_conv_crystal = [[i['label']] for i in SymmPrimFramework.as_dict()['sites']]
            labels_conv_crystal = np.array(labels_conv_crystal)

            ion_conv_crystal = [i['abc'] for i in SymmPrimFramework.as_dict()['sites']]
            ion_conv_crystal = np.array(ion_conv_crystal).astype(float)

            cell = np.array(SymmPrimFramework.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal*(1, 1, 0.5)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1.5) + (1/4, 1/4, 0))

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_1 = Structure(lattice, AB_label + AB_label, AB, coords_are_cartesian=False)
            AB_1_symm = SpacegroupAnalyzer(AB_1,
                                           symprec=self.symm_tol,
                                           angle_tolerance=self.angle_tol)

            self.symm_structure = AB_1_symm.get_refined_structure()

        if stacking == 'AB2':
            labels_conv_crystal = [[i['label']] for i in SymmPrimFramework.as_dict()['sites']]
            labels_conv_crystal = np.array(labels_conv_crystal)

            ion_conv_crystal = [i['abc'] for i in SymmPrimFramework.as_dict()['sites']]
            ion_conv_crystal = np.array(ion_conv_crystal).astype(float)

            cell = np.array(SymmPrimFramework.as_dict()['lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal*(1, 1, 0.5)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1.5) + (1/2, 0, 0))

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_2 = Structure(lattice, AB_label + AB_label, AB, coords_are_cartesian=False)
            AB_2_symm = SpacegroupAnalyzer(AB_2,
                                           symprec=self.symm_tol,
                                           angle_tolerance=self.angle_tol)

            self.symm_structure = AB_2_symm.get_refined_structure()

        # Create AAl stacking.
        # Tetragonal cell with two sheets per cell shifited by the shift_vector in angstroms.
        if stacking == 'AAl':
            labels_conv_crystal = [[i['label']] for i in SymmPrimFramework.as_dict()['sites']]
            labels_conv_crystal = np.array(labels_conv_crystal)

            ion_conv_crystal = [i['abc'] for i in SymmPrimFramework.as_dict()['sites']]
            ion_conv_crystal = np.array(ion_conv_crystal).astype(float)

            cell = np.array(SymmPrimFramework.as_dict()['lattice']['matrix'])*(1, 1, 2)

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal * np.array([1, 1, 0.5])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(
                *Tools.cell_to_cellpar(cell))
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
            labels_conv_crystal = [[i['label']] for i in SymmPrimFramework.as_dict()['sites']]
            labels_conv_crystal = np.array(labels_conv_crystal)

            ion_conv_crystal = [i['abc'] for i in SymmPrimFramework.as_dict()['sites']]
            ion_conv_crystal = np.array(ion_conv_crystal).astype(float)

            cell = SymmPrimFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_parameter_base = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            new_cell = Tools.cellpar_to_cell(
                [a_cell, b_cell, c_parameter_base, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + \
                np.array([0, 0, 0.5*c_parameter_base]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(new_cell)
            self.symm_structure = Structure(
                lattice, AB_label+AB_label, AB, coords_are_cartesian=True)

        if stacking == 'ABC1':
            labels_conv_crystal = [[i['label']] for i in SymmPrimFramework.as_dict()['sites']]
            labels_conv_crystal = np.array(labels_conv_crystal)

            ion_conv_crystal = [i['abc'] for i in SymmPrimFramework.as_dict()['sites']]
            ion_conv_crystal = np.array(ion_conv_crystal).astype(float)

            cell = np.array(SymmPrimFramework.as_dict()['lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (2/3, 1/3, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (4/3, 2/3, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice,
                              ABC_label + ABC_label + ABC_label,
                              ABC,
                              coords_are_cartesian=False)

            ABC_f_symm = SpacegroupAnalyzer(ABC_f,
                                            symprec=self.symm_tol,
                                            angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        if stacking == 'ABC2':
            labels_conv_crystal = np.array(
                [[i['label']] for i in SymmPrimFramework.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in SymmPrimFramework.as_dict()['sites']])
            cell = np.array(SymmPrimFramework.as_dict()[
                            'lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (1/3, 0, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (2/3, 0, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label +
                              ABC_label, ABC, coords_are_cartesian=False)
            ABC_f_symm = SpacegroupAnalyzer(
                ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        dict_structure = self.symm_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

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

    def _create_kgm_structure(self,
                              name_bb_a: str,
                              name_bb_b: str,
                              stacking: str = 'AA',
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
            print('Building block A must present connectivity 4 insted of',
                  len(bb_1.connectivity))
            return None
        if bb_2.connectivity != 4:
            print('Building block B must present connectivity 4 insted of',
                  len(bb_2.connectivity))
            return None

        # Detect the bond atom from the connection groups type
        bond_atom = Tools.get_bond_atom(bb_1.conector, bb_2.conector)

        Tools.print_comand(f'Bond atom detected: {bond_atom}',
                           self.verbosity,
                           ['debug', 'high'])

        # Calculate the cell parameter based on the size of the building blocks
        if self.verbosity:
            print('BB_A size:', bb_1.size)

        if self.verbosity:
            print('BB_B size:', bb_2.size)

        size_a = max(bb_1.size)
        size_b = max(bb_2.size)
        # Defines the cell parameter a
        a = 2*(size_a + size_b)  # np.cos(np.radians(30))*2*(size_a + size_b)

        if self.verbosity:
            print('Calculated cell parameter a:', a)

        # Gets the maximum distace in the z axis to create the c parameter
        delta_a = abs(max(np.transpose(bb_1.atom_pos)[
                      2])) + abs(min(np.transpose(bb_1.atom_pos)[2]))
        delta_b = abs(max(np.transpose(bb_2.atom_pos)[
                      2])) + abs(min(np.transpose(bb_2.atom_pos)[2]))

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

        bb1_v1, bb1_v2, bb1_v3, bb1_v4 = bb_1._get_X_points()[1]
        bb1_lado1, _, _ = (
            bb1_v1 + bb1_v2) / 2, (bb1_v1 + bb1_v3) / 2, (bb1_v1 + bb1_v4) / 2

        R_matrix = Tools.rotation_matrix_from_vectors(
            bb1_lado1, np.array([1, 0, 0]))

        bb_1.atom_pos = np.dot(bb_1.atom_pos, np.transpose(R_matrix))

        # Add tbe building block 1 (C4) on the center [1/2, 1/2, 0.0] of the unitary cell (A1 site)
        final_pos = np.dot(bb_1.atom_pos, R.from_euler(
            'z', 30, degrees=True).as_matrix()) + np.array([1/4, np.sqrt(3)/4, 0])*a
        final_label = list(bb_1.atom_labels)

        # Add tbe building block 1 (C4) on [1/2, 0.0, 0.0] of the unitary cell (A2 site) -45
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler(
            'z', -90, degrees=True).as_matrix()) + np.array([1/2, 0, 0])*a))
        final_label += list(bb_1.atom_labels)

        # Add tbe building block 1 (C4) on [0.0, 1/2, 0.0] of the unitary cell (A3 site) 30
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler(
            'z', -30, degrees=True).as_matrix()) + np.array([-1/4, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_1.atom_labels)

        # Changes the X atoms by the desired bond_atom
        final_label, final_pos = Tools.change_X_atoms(
            final_label, final_pos, bond_atom)

        # Cria a estrutura como entidade do pymatgen
        struct = Structure(lattice, final_label, final_pos,
                           coords_are_cartesian=True)

        # Remove os átomos duplicados
        struct.sort(reverse=True)
        struct.merge_sites(tol=.7, mode='delete')
        struct.translate_sites(range(len(struct.as_dict()['sites'])), [
                               0, 0, 0.5], frac_coords=True, to_unit_cell=True)

        # atom_labels = np.array([[i['label']]
        #                       for i in struct.as_dict()['sites']]).flatten()
        # atom_pos = np.array([i['xyz'] for i in struct.as_dict()['sites']])
        cell = np.array(struct.as_dict()['lattice']['matrix'])

        # Change the X atoms by the desired bond_atom
        final_label, final_pos = Tools.change_X_atoms(
            final_label, final_pos, bond_atom)

        # Cria a estrutura como entidade do pymatgen
        struct = Structure(lattice, final_label, final_pos,
                           coords_are_cartesian=True)

        # Remove os átomos duplicados
        struct.sort(reverse=True)
        struct.merge_sites(tol=.7, mode='delete')
        struct.translate_sites(range(len(struct.as_dict()['sites'])), [
                               0, 0, 0.5], frac_coords=True, to_unit_cell=True)

        # atom_labels = np.array([[i['label']]
        #                       for i in struct.as_dict()['sites']]).flatten()
        # atom_pos = np.array([i['xyz'] for i in struct.as_dict()['sites']])
        cell = np.array(struct.as_dict()['lattice']['matrix'])

        # Simetriza a estrutura
        symm = SpacegroupAnalyzer(
            struct, symprec=self.symm_tol, angle_tolerance=self.angle_tol)
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
            labels_conv_crystal = np.array(
                [[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()[
                            'lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal*(1, 1, 0.5)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1.5) + (1/4, 1/4, 0))

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_1 = Structure(lattice, AB_label + AB_label,
                             AB, coords_are_cartesian=False)
            AB_1_symm = SpacegroupAnalyzer(
                AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = AB_1_symm.get_refined_structure()

        if stacking == 'AB2':
            self.stacking = 'AB2'
            labels_conv_crystal = np.array(
                [[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()[
                            'lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal*(1, 1, 0.5)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1.5) + (1/2, 0, 0))

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_2 = Structure(lattice, AB_label + AB_label,
                             AB, coords_are_cartesian=False)
            AB_2_symm = SpacegroupAnalyzer(
                AB_2, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = AB_2_symm.get_refined_structure()

        # Create AAl stacking.
        # Tetragonal cell with two sheets per cell shifited by the shift_vector in angstroms.
        if stacking == 'AAl':
            self.stacking = 'AAl'
            labels_conv_crystal = np.array(
                [[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()[
                            'lattice']['matrix'])*(1, 1, 2)

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal * np.array([1, 1, 0.5])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(
                *Tools.cell_to_cellpar(cell))
            shift_vector = np.dot(r, np.array(shift_vector))

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal * np.array([1, 1, 1.5]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            self.symm_structure = Structure(
                lattice, AB_label+AB_label, AB, coords_are_cartesian=False)

        # Create AA tilted stacking.
        # Tilted tetragonal cell with two sheets per cell tilted by tilt_angle.
        if stacking == 'AAt':
            self.stacking = 'AAt'

            labels_conv_crystal = np.array(
                [[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['xyz'] for i in struct_symm_prim.as_dict()['sites']])
            cell = struct_symm_prim.as_dict()['lattice']

            # Shift the cell by the tilt angle
            print(cell)
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            new_cell = Tools.cellpar_to_cell(
                [a_cell, b_cell, c_cell, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + np.array([0, 0, 0.5*c_cell]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(new_cell)
            self.symm_structure = Structure(
                lattice, AB_label+AB_label, AB, coords_are_cartesian=True)

        if stacking == 'ABC1':
            self.stacking = 'ABC1'
            labels_conv_crystal = np.array(
                [[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()[
                            'lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (2/3, 1/3, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (4/3, 2/3, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label +
                              ABC_label, ABC, coords_are_cartesian=False)
            ABC_f_symm = SpacegroupAnalyzer(
                ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        if stacking == 'ABC2':
            self.stacking = 'ABC2'
            labels_conv_crystal = np.array(
                [[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()[
                            'lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1) + (1/3, 0, 0))
            C = Tools.translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (2/3, 0, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label +
                              ABC_label, ABC, coords_are_cartesian=False)
            ABC_f_symm = SpacegroupAnalyzer(
                ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = ABC_f_symm.get_refined_structure()

        dict_structure = self.symm_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

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

    def _create_kgm_a_structure(self,
                                name_bb_a: str,
                                name_bb_b: str,
                                stacking: str = 'AA',
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
            print('Building block A must present connectivity 4 insted of',
                  len(bb_1.connectivity))
            return None
        if bb_2.connectivity != 2:
            print('Building block B must present connectivity 2 insted of',
                  len(bb_2.connectivity))
            return None

        # Detect the bond atom from the connection groups type
        bond_atom = Tools.get_bond_atom(bb_1.conector, bb_2.conector)

        Tools.print_comand(f'Bond atom detected: {bond_atom}',
                           self.verbosity,
                           ['debug', 'high'])

        # Calculate the cell parameter based on the size of the building blocks
        size_a = max(bb_1.size)
        if self.verbosity:
            print('BB_A size:', size_a)

        size_b = max(bb_2.size)
        if self.verbosity:
            print('BB_B size:', size_b)

        # Defines the cell parameter a
        a = 4*(size_a + size_b)  # np.cos(np.radians(30))*2*(size_a + size_b)

        if self.verbosity:
            print('Calculated cell parameter a:', a)

        # Gets the maximum distace in the z axis to create the c parameter
        delta_a = abs(max(np.transpose(bb_1.atom_pos)[
                      2])) + abs(min(np.transpose(bb_1.atom_pos)[2]))
        delta_b = abs(max(np.transpose(bb_2.atom_pos)[
                      2])) + abs(min(np.transpose(bb_2.atom_pos)[2]))

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

        bb1_v1, bb1_v2, bb1_v3, bb1_v4 = bb_1._get_X_points()[1]
        bb1_lado1, _, _ = (
            bb1_v1 + bb1_v2) / 2, (bb1_v1 + bb1_v3) / 2, (bb1_v1 + bb1_v4) / 2

        R_matrix = Tools.rotation_matrix_from_vectors(
            bb1_lado1, np.array([1, 0, 0]))

        bb_1.atom_pos = np.dot(bb_1.atom_pos, np.transpose(R_matrix))

        # Add tbe building block 1 (C4) on the center [1/2, 1/2, 1/2] of the unitary cell (A1 site)
        final_pos = np.dot(bb_1.atom_pos, R.from_euler(
            'z', 30, degrees=True).as_matrix()) + np.array([1/4, np.sqrt(3)/4, 0])*a
        final_label = list(bb_1.atom_labels)

        # Add tbe building block 1 (C4) on [1/2, 0.0, 0.0] of the unitary cell (A2 site)
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler(
            'z', -90, degrees=True).as_matrix()) + np.array([1/2, 0, 0])*a))
        final_label += list(bb_1.atom_labels)

        # Add tbe building block 1 (C4) on [0.0, 1/2, 0.0] of the unitary cell (A3 site)
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler(
            'z', -30, degrees=True).as_matrix()) + np.array([-1/4, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_1.atom_labels)

        '''
        Older version of the positioning code
        # Add tbe building block 1 (C4) on the center of the unitary cell (A1 site)
        final_pos = np.dot(bb_1.atom_pos, R.from_euler('z', -30, degrees=True).as_matrix())
         + np.array([1/4, np.sqrt(3)/4, 0])*a
        final_label = list(bb_1.atom_labels)

        # Add tbe building block 1 (C4) on [0.5, 0.0, 0.0] of the unitary cell (A2 site)
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos,
        R.from_euler('z', -60, degrees=True).as_matrix()) + np.array([1/2, 0, 0])*a))
        final_label += list(bb_1.atom_labels)

        # Add tbe building block 1 (C4) on [0.0, 0.5, 0.0] of the unitary cell (A3 site)
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos,
        R.from_euler('z', 60, degrees=True).as_matrix()) + np.array([-1/4, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_1.atom_labels)'''

        # Add the building block 2 (C2) on [1/2, 1/4, 0.0] of the unitary cell (B1 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler(
            'z', -30, degrees=True).as_matrix()) + np.array([3/8, np.sqrt(3)/8, 0])*a))
        final_label += list(bb_2.atom_labels)

        # Add the building block 2 (C2) on [1/4, 3/4, 0.0] of the unitary cell (B2 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler(
            'z', -30, degrees=True).as_matrix()) + np.array([1/8, 3*np.sqrt(3)/8, 0])*a))
        final_label += list(bb_2.atom_labels)

        # Add tbe building block 2 (C2) on [3/4, 1/4, 0.0] of the unitary cell (B3 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler(
            'z', 30, degrees=True).as_matrix()) + np.array([5/8, np.sqrt(3)/8, 0])*a))
        final_label += list(bb_2.atom_labels)

        # Add tbe building block 2 (C2) on [1/4, 3/4, 0.0] of the unitary cell (B4 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler(
            'z', 30, degrees=True).as_matrix()) + np.array([-1/8, 3*np.sqrt(3)/8, 0])*a))
        final_label += list(bb_2.atom_labels)

        # Add tbe building block 2 (C2) on [3/4, 1/2, 0.0] of the unitary cell (B5 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler(
            'z', 90, degrees=True).as_matrix()) + np.array([4/8, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_2.atom_labels)

        # Add tbe building block 2 (C2) on [1/2, 3/4, 0.0] of the unitary cell (B6 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler(
            'z', 90, degrees=True).as_matrix()) + np.array([0.0, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_2.atom_labels)

        # Change the X atoms by the desired bond_atom
        final_label, final_pos = Tools.change_X_atoms(
            final_label, final_pos, bond_atom)

        # Creates a pymatgen structure
        struct = Structure(lattice, final_label, final_pos,
                           coords_are_cartesian=True)

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

            lattice = Lattice(
                np.array(struct.as_dict()['lattice']['matrix'])*(1, 1, 2))

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

            self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

            self.atom_labels = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(
                AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        if stacking == 'AB1y':
            self.stacking = 'AB1y'

            lattice = Lattice(
                np.array(struct.as_dict()['lattice']['matrix'])*(1, 1, 2))

            A = struct
            A_labels = np.array([i['label'] for i in A.as_dict()['sites']])
            A_pos = np.array([i['xyz'] for i in A.as_dict()['sites']])

            B = struct
            B_labels = np.array([i['label'] for i in B.as_dict()['sites']])
            B_pos = np.array([i['xyz'] for i in B.as_dict()['sites']])
            B_pos += np.array([lattice.a*1/4, lattice.b*1 /
                              4*np.sqrt(3), lattice.c*1/2])

            AB_1 = Structure(lattice,
                             np.concatenate((A_labels, B_labels)),
                             np.concatenate((A_pos, B_pos)),
                             coords_are_cartesian=True)

            dict_structure = AB_1.as_dict()

            self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

            self.atom_labels = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(
                AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        if stacking == 'AB1xy':
            self.stacking = 'AB1xy'

            lattice = Lattice(
                np.array(struct.as_dict()['lattice']['matrix'])*(1, 1, 2))

            A = struct
            A_labels = np.array([i['label'] for i in A.as_dict()['sites']])
            A_pos = np.array([i['xyz'] for i in A.as_dict()['sites']])

            B = struct
            B_labels = np.array([i['label'] for i in B.as_dict()['sites']])
            B_pos = np.array([i['xyz'] for i in B.as_dict()['sites']])
            B_pos += np.array([lattice.a*1/4, -lattice.b *
                              1/4*np.sqrt(3), lattice.c*1/2])

            AB_1 = Structure(lattice,
                             np.concatenate((A_labels, B_labels)),
                             np.concatenate((A_pos, B_pos)),
                             coords_are_cartesian=True)

            dict_structure = AB_1.as_dict()

            self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

            self.atom_labels = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(
                AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        if stacking == 'AB2':
            self.stacking = 'AB2'

            lattice = Lattice(
                np.array(struct.as_dict()['lattice']['matrix'])*(1, 1, 2))

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

            self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

            self.atom_labels = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(
                AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        # Create AAl stacking.
        # Tetragonal cell with two sheets per cell shifited by the shift_vector in angstroms.
        if stacking == 'AAl':
            self.stacking = 'AAl'

            lattice = Lattice(
                np.array(struct.as_dict()['lattice']['matrix'])*(1, 1, 2))

            A = struct
            A_labels = np.array([i['label'] for i in A.as_dict()['sites']])
            A_pos = np.array([i['xyz'] for i in A.as_dict()['sites']])

            B = struct
            B_labels = np.array([i['label'] for i in B.as_dict()['sites']])
            B_pos = np.array([i['xyz'] for i in B.as_dict()['sites']])

            # Calculates the shift vector in crystal units
            r = Tools.get_cartesian_to_fractional_matrix(*lattice.parameters)
            shift_vector = np.dot(r, np.array(shift_vector))

            B_pos += np.array([lattice.a, lattice.b*1/2 *
                              np.sqrt(3), 1]) * shift_vector
            B_pos += np.array([0, 0, lattice.c*1/2])

            AB_1 = Structure(lattice,
                             np.concatenate((A_labels, B_labels)),
                             np.concatenate((A_pos, B_pos)),
                             coords_are_cartesian=True)

            dict_structure = AB_1.as_dict()

            self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

            self.atom_labels = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(
                AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        # Create AA tilted stacking.
        # Tilted tetragonal cell with two sheets per cell tilted by tilt_angle.
        if stacking == 'AAt':
            self.stacking = 'AAt'

            a_cell = struct.as_dict()['lattice']['a']
            b_cell = struct.as_dict()['lattice']['b']
            c_cell = struct.as_dict()['lattice']['c'] * 2
            alpha = struct.as_dict()['lattice']['alpha'] - tilt_angle
            beta = struct.as_dict()['lattice']['beta'] - tilt_angle
            gamma = struct.as_dict()['lattice']['gamma']

            new_cell = Tools.cellpar_to_cell(
                [a_cell, b_cell, c_cell, alpha, beta, gamma])

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

            B_pos += np.array([lattice.a, lattice.b*1/2 *
                              np.sqrt(3), 1]) * shift_vector
            B_pos += np.array([0, 0, lattice.c*1/2])

            AB_1 = Structure(lattice,
                             np.concatenate((A_labels, B_labels)),
                             np.concatenate((A_pos, B_pos)),
                             coords_are_cartesian=True)

            dict_structure = AB_1.as_dict()

            self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

            self.atom_labels = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(
                AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.symm_structure = symm.get_refined_structure()

        dict_structure = self.symm_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

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
