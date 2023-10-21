# -*- coding: utf-8 -*-
# Copyright (c) Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

"""
This class implements definitions for a Framework buiding
"""

import os
import numpy as np

# Import pymatgen
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from scipy.spatial.transform import Rotation as R

from pycofbuilder.tools import (print_command,
                                get_bond_atom,
                                translate_inside,
                                get_cartesian_to_fractional_matrix,
                                cell_to_cellpar,
                                cellpar_to_cell,
                                rotation_matrix_from_vectors,
                                change_X_atoms,
                                print_framework_name)

from pycofbuilder.io_tools import (save_json,
                                   save_chemjson,
                                   save_cif,
                                   save_xyz,
                                   save_turbomole,
                                   save_vasp,
                                   save_xsf,
                                   save_pdb,
                                   save_pqr,
                                   save_qe)

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
    atom_types : list = []
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

    def __init__(self, name=None, verbosity=False, out_dir=None, save_bb=False):

        _ROOTDIR = os.path.abspath(os.path.dirname(__file__))

        self.verbosity = verbosity
        self.main_path = os.path.join(_ROOTDIR, 'data')

        if out_dir is None:
            self.out_path = os.path.join(os.getcwd(), 'out')
        else:
            self.out_path = out_dir

        self.save_bb = save_bb

        self.lib_path = os.path.join(self.out_path, 'building_blocks')

        self.BB1_name = None
        self.BB2_name = None
        self.topology = None
        self.stacking = None
        self.smiles = None

        self.atom_types = []
        self.atom_pos = []
        self.atom_labels = []
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
        self.available_2D_top = ['HCB', 'HCB_A',
                                 'SQL', 'SQL_A',
                                 'KGD',
                                 'HXL', 'HXL_A',
                                 'FXT', 'FXT_A']

        # Falta add: ['dia', 'bor', 'srs', 'pts', 'ctn', 'rra', 'fcc', 'lon', 'stp', 'acs', 'tbo', 'bcu', 'fjh', 'ceq']
        self.available_3D_top = ['DIA', 'BOR']  # Temporary
        self.available_topologies = self.available_2D_top + self.available_3D_top

        # Define available stackings for all 2D topologies
        self.available_stacking = {
            'HCB': ['A', 'AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
            'HCB_A': ['A', 'AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
            'SQL': ['A', 'AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
            'SQL_A': ['A', 'AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
            'KGD': ['A', 'AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
            'HXL_A': ['A', 'AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
            'FXT': ['A', 'AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
            'FXT_A': ['A', 'AA', 'AB1', 'AB2', 'AAl', 'AAt', 'ABC1', 'ABC2'],
            'DIA': [0, 1, 2, 3, 4],  # Temporary
            'BOR': [0, 5, 8, 6, 7]  # Temporary
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
        return len(self.atom_types)

    def get_available_topologies(self, dimensionality='all'):

        if dimensionality == 'all' or dimensionality == '2D':
            print('Available 2D Topologies:')
            for i in self.available_2D_top:
                print(i.upper())

        if dimensionality == 'all' or dimensionality == '3D':
            print('Available 3D Topologies:')
            for i in self.available_3D_top:
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

        BB1 = Building_Block(name=BB1_name, save_dir=self.out_path, save_bb=self.save_bb)
        BB2 = Building_Block(name=BB2_name, save_dir=self.out_path, save_bb=self.save_bb)

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
            'KGD': self._create_kgd_structure,
            # 'HXL': self._create_hxl_structure,
            'HXL_A': self._create_hxl_a_structure,
            'FXT': self._create_fxt_structure,
            'FXT_A': self._create_fxt_a_structure,
            }

        result = net_build_dict[Net](BB1, BB2, stacking=Stacking, **kwargs)

        return result

    def save(self, fmt: str = 'cif', supercell: list = [1, 1, 1], save_dir=None, primitive=False) -> None:
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
            'json': save_json,
            'cjson': save_chemjson,
            'cif': save_cif,
            'xyz': save_xyz,
            'turbomole': save_turbomole,
            'vasp': save_vasp,
            'xsf': save_xsf,
            'pdb': save_pdb,
            'pqr': save_pqr,
            'qe': save_qe
         }

        file_format_error = f'Format must be one of the following: {save_dict.keys()}'
        assert fmt in save_dict.keys(), file_format_error

        if primitive:
            structure = self.prim_structure

        else:
            structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

        final_structure = structure.make_supercell(supercell, in_place=False)

        structure_dict = final_structure.as_dict()

        cell = structure_dict['lattice']['matrix']

        atom_types = [site['species'][0]['element'] for site in structure_dict['sites']]
        atom_labels = [site['properties']['source'] for site in structure_dict['sites']]
        atom_pos = [site['xyz'] for site in structure_dict['sites']]

        if save_dir is None:
            save_path = self.out_path
        else:
            save_path = save_dir

        save_dict[fmt](path=save_path,
                       file_name=self.name,
                       cell=cell,
                       atom_types=atom_types,
                       atom_labels=atom_labels,
                       atom_pos=atom_pos)

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

        self.name = f'{BB_T3_A.name}-{BB_T3_B.name}-HCB-{stacking}'
        self.topology = 'HCB'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_T3_A.charge + BB_T3_B.charge
        self.chirality = BB_T3_A.chirality or BB_T3_B.chirality

        print_command(f'Starting the creation of {self.name}', self.verbosity, ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_T3_A.conector, BB_T3_B.conector)

        # Replace "X" the building block
        BB_T3_A.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_T3_A.remove_X()
        BB_T3_B.remove_X()

        print_command(f'Bond atom detected: {bond_atom}', self.verbosity, ['debug', 'high'])

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

        if self.stacking == 'A':
            c = slab

        # Create the lattice
        self.lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        self.atom_types = []
        self.atom_labels = []
        self.atom_pos = []

        # Add the A1 building blocks to the structure
        vertice_data = topology_info['vertices'][0]
        self.atom_types += BB_T3_A.atom_types
        vertice_pos = np.array(vertice_data['position'])*a

        R_Matrix = R.from_euler('z',
                                vertice_data['angle'],
                                degrees=True).as_matrix()

        rotated_pos = np.dot(BB_T3_A.atom_pos, R_Matrix) + vertice_pos
        self.atom_pos += rotated_pos.tolist()

        self.atom_labels += ['C1' if i == 'C' else i for i in BB_T3_A.atom_labels]

        # Add the A2 building block to the structure
        vertice_data = topology_info['vertices'][1]
        self.atom_types += BB_T3_B.atom_types
        vertice_pos = np.array(vertice_data['position'])*a

        R_Matrix = R.from_euler('z',
                                vertice_data['angle'],
                                degrees=True).as_matrix()

        rotated_pos = np.dot(BB_T3_B.atom_pos, R_Matrix) + vertice_pos
        self.atom_pos += rotated_pos.tolist()

        self.atom_labels += ['C2' if i == 'C' else i for i in BB_T3_B.atom_labels]

        # Creates a pymatgen structure
        StartingFramework = Structure(
            self.lattice,
            self.atom_types,
            self.atom_pos,
            coords_are_cartesian=True,
            site_properties={'source': self.atom_labels}
        ).get_sorted_structure()

        # Translates the structure to the center of the cell
        StartingFramework.translate_sites(
            range(len(StartingFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True)

        dict_structure = StartingFramework.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [2/3, 1/3, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'AB2':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [1/2, 0, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [1/2, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'ABC1':
            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (2/3, 1/3, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (4/3, 2/3, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        if stacking == 'ABC2':
            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (1/3, 0, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (2/3, 0, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        if stacking == 'AAl':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        # Create AA tilted stacking.
        if stacking == 'AAt':
            cell = StartingFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            self.lattice = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        dict_structure = stacked_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure, symprec=0.1, angle_tolerance=.5)

        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            print_command(self.prim_structure, self.verbosity, ['debug'])

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()
        except Exception as e:
            print_command(e, self.verbosity, ['debug'])

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        if print_result is True:
            print_framework_name(self.name,
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

        self.name = f'{BB_T3.name}-{BB_L2.name}-HCB_A-{stacking}'
        self.topology = 'HCB_A'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_L2.charge + BB_T3.charge
        self.chirality = BB_L2.chirality or BB_T3.chirality

        print_command(f'Starting the creation of {self.name}', self.verbosity, ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_T3.conector, BB_L2.conector)

        # Replace "X" the building block
        BB_L2.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_T3.remove_X()
        BB_L2.remove_X()

        print_command(f'Bond atom detected: {bond_atom}', self.verbosity, ['debug', 'high'])

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

        if self.stacking == 'A':
            c = slab

        # Create the lattice
        self.lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        self.atom_types = []
        self.atom_labels = []
        self.atom_pos = []

        # Add the building blocks to the structure
        for vertice_data in topology_info['vertices']:
            self.atom_types += BB_T3.atom_types
            vertice_pos = np.array(vertice_data['position'])*a

            R_Matrix = R.from_euler('z', vertice_data['angle'], degrees=True).as_matrix()

            rotated_pos = np.dot(BB_T3.atom_pos, R_Matrix) + vertice_pos
            self.atom_pos += rotated_pos.tolist()

            self.atom_labels += ['C1' if i == 'C' else i for i in BB_T3.atom_labels]

        # Add the building blocks to the structure
        for edge_data in topology_info['edges']:
            self.atom_types += BB_L2.atom_types

            R_Matrix = R.from_euler('z', edge_data['angle'], degrees=True).as_matrix()

            rotated_pos = np.dot(BB_L2.atom_pos, R_Matrix) + np.array(edge_data['position'])*a

            self.atom_pos += rotated_pos.tolist()

            self.atom_labels += ['C2' if i == 'C' else i for i in BB_L2.atom_labels]

        StartingFramework = Structure(
            self.lattice,
            self.atom_types,
            self.atom_pos,
            coords_are_cartesian=True,
            site_properties={'source': self.atom_labels}
        ).get_sorted_structure()

        # Translates the structure to the center of the cell
        StartingFramework.translate_sites(
            range(len(StartingFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True
        )

        dict_structure = StartingFramework.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [2/3, 1/3, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'AB2':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [1/2, 0, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [1/2, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'ABC1':
            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (2/3, 1/3, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (4/3, 2/3, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        if stacking == 'ABC2':
            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (1/3, 0, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (2/3, 0, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        if stacking == 'AAl':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        # Create AA tilted stacking.
        if stacking == 'AAt':
            cell = StartingFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            self.lattice = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        dict_structure = stacked_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure, symprec=0.05, angle_tolerance=.5)
        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            print_command(self.prim_structure, self.verbosity, ['debug'])

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()

        except Exception as e:
            print_command(e, self.verbosity, ['debug'])

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        if print_result is True:
            print_framework_name(self.name,
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

        self.name = f'{BB_S4_A.name}-{BB_S4_B.name}-SQL-{stacking}'
        self.topology = 'SQL'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_S4_A.charge + BB_S4_B.charge
        self.chirality = BB_S4_A.chirality or BB_S4_B.chirality

        print_command(f'Starting the creation of {self.name}', self.verbosity, ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_S4_A.conector, BB_S4_B.conector)

        # Replace "X" the building block
        BB_S4_A.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_S4_A.remove_X()
        BB_S4_B.remove_X()

        print_command(f'Bond atom detected: {bond_atom}', self.verbosity, ['debug', 'high'])

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

        if self.stacking == 'A':
            c = slab

        # Create the lattice
        self.lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        self.atom_types = []
        self.atom_labels = []
        self.atom_pos = []

        # Add the first building block to the structure
        vertice_data = topology_info['vertices'][0]
        self.atom_types += BB_S4_A.atom_types
        vertice_pos = np.array(vertice_data['position'])*a

        R_Matrix = R.from_euler('z', vertice_data['angle'], degrees=True).as_matrix()

        rotated_pos = np.dot(BB_S4_A.atom_pos, R_Matrix) + vertice_pos
        self.atom_pos += rotated_pos.tolist()

        self.atom_labels += ['C1' if i == 'C' else i for i in BB_S4_A.atom_labels]

        # Add the second building block to the structure
        vertice_data = topology_info['vertices'][1]
        self.atom_types += BB_S4_B.atom_types
        vertice_pos = np.array(vertice_data['position'])*a

        R_Matrix = R.from_euler('z', vertice_data['angle'], degrees=True).as_matrix()

        rotated_pos = np.dot(BB_S4_B.atom_pos, R_Matrix) + vertice_pos
        self.atom_pos += rotated_pos.tolist()

        self.atom_labels += ['C2' if i == 'C' else i for i in BB_S4_B.atom_labels]

        StartingFramework = Structure(
            self.lattice,
            self.atom_types,
            self.atom_pos,
            coords_are_cartesian=True,
            site_properties={'source': self.atom_labels}
        ).get_sorted_structure()

        # Translates the structure to the center of the cell
        StartingFramework.translate_sites(
            range(len(StartingFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True
        )

        dict_structure = StartingFramework.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':

            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [1/4, 1/4, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [1/4, 1/4, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'AB2':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [1/2, 0, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [1/2, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'ABC1':

            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (1/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (1/3, 1/3, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 2/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (2/3, 2/3, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        if stacking == 'ABC2':
            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (1/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (1/3, 0, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 2/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (2/3, 0, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        # Create AAl stacking. Tetragonal cell with two sheets
        # per cell shifited by the shift_vector in angstroms.
        if stacking == 'AAl':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        # Create AA tilted stacking.
        # Tilted tetragonal cell with two sheets per cell tilted by tilt_angle.
        if stacking == 'AAt':
            cell = StartingFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            self.lattice = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        dict_structure = stacked_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure, symprec=0.05, angle_tolerance=.5)
        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            print_command(self.prim_structure, self.verbosity, ['debug'])

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()

        except Exception as e:
            print_command(e, self.verbosity, ['debug'])

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        if print_result is True:
            print_framework_name(self.name,
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

        self.name = f'{BB_S4.name}-{BB_L2.name}-SQL_A-{stacking}'
        self.topology = 'SQL_A'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_S4.charge + BB_L2.charge
        self.chirality = BB_S4.chirality or BB_L2.chirality

        print_command(f'Starting the creation of {self.name}', self.verbosity, ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_S4.conector, BB_L2.conector)

        # Replace "X" the building block
        BB_L2.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_S4.remove_X()
        BB_L2.remove_X()

        print_command(f'Bond atom detected: {bond_atom}', self.verbosity, ['debug', 'high'])

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

        if self.stacking == 'A':
            c = slab

        # Create the lattice
        self.lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        self.atom_types = []
        self.atom_labels = []
        self.atom_pos = []

        # Add the building blocks to the structure
        for vertice_data in topology_info['vertices']:
            self.atom_types += BB_S4.atom_types
            vertice_pos = np.array(vertice_data['position'])*a

            R_Matrix = R.from_euler('z', vertice_data['angle'], degrees=True).as_matrix()

            rotated_pos = np.dot(BB_S4.atom_pos, R_Matrix) + vertice_pos
            self.atom_pos += rotated_pos.tolist()

            self.atom_labels += ['C1' if i == 'C' else i for i in BB_S4.atom_labels]

        # Add the building blocks to the structure
        for edge_data in topology_info['edges']:
            self.atom_types += BB_L2.atom_types
            edge_pos = np.array(edge_data['position'])*a

            R_Matrix = R.from_euler('z', edge_data['angle'], degrees=True).as_matrix()

            rotated_pos = np.dot(BB_L2.atom_pos, R_Matrix) + edge_pos
            self.atom_pos += rotated_pos.tolist()

            self.atom_labels += ['C2' if i == 'C' else i for i in BB_L2.atom_labels]

        StartingFramework = Structure(
            self.lattice,
            self.atom_types,
            self.atom_pos,
            coords_are_cartesian=True,
            site_properties={'source': self.atom_labels}
        ).get_sorted_structure()

        # Translates the structure to the center of the cell
        StartingFramework.translate_sites(
            range(len(StartingFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True
        )

        dict_structure = StartingFramework.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':

            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [1/4, 1/4, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [1/4, 1/4, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'AB2':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [1/2, 0, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [1/2, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'ABC1':

            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (1/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (1/3, 1/3, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 2/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (2/3, 2/3, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        if stacking == 'ABC2':
            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (1/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (1/3, 0, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 2/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (2/3, 0, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        # Create AAl stacking. Tetragonal cell with two sheets
        # per cell shifited by the shift_vector in angstroms.
        if stacking == 'AAl':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        # Create AA tilted stacking.
        # Tilted tetragonal cell with two sheets per cell tilted by tilt_angle.
        if stacking == 'AAt':
            cell = StartingFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            self.lattice = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        dict_structure = stacked_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure, symprec=0.05, angle_tolerance=.5)
        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            print_command(self.prim_structure, self.verbosity, ['debug'])

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()

        except Exception as e:
            print_command(e, self.verbosity, ['debug'])

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        if print_result is True:
            print_framework_name(self.name,
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

        self.name = f'{BB_H6.name}-{BB_T3.name}-KGD-{stacking}'
        self.topology = 'KGD'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_H6.charge + BB_T3.charge
        self.chirality = BB_H6.chirality or BB_T3.chirality

        print_command(f'Starting the creation of {self.name}', self.verbosity, ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_H6.conector, BB_T3.conector)

        # Replace "X" the building block
        BB_H6.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_H6.remove_X()
        BB_T3.remove_X()

        print_command(f'Bond atom detected: {bond_atom}', self.verbosity, ['debug', 'high'])

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

        if self.stacking == 'A':
            c = slab

        # Create the lattice
        self.lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        self.atom_types = []
        self.atom_labels = []
        self.atom_pos = []

        # Add the building blocks to the structure
        for vertice_data in topology_info['vertices']:
            self.atom_types += BB_H6.atom_types
            vertice_pos = np.array(vertice_data['position'])*a

            R_Matrix = R.from_euler('z',
                                    vertice_data['angle'],
                                    degrees=True).as_matrix()

            rotated_pos = np.dot(BB_H6.atom_pos, R_Matrix) + vertice_pos
            self.atom_pos += rotated_pos.tolist()

            self.atom_labels += ['C1' if i == 'C' else i for i in BB_H6.atom_labels]

        # Add the building blocks to the structure
        for edge_data in topology_info['edges']:
            self.atom_types += BB_T3.atom_types
            edge_pos = np.array(edge_data['position'])*a

            R_Matrix = R.from_euler('z',
                                    edge_data['angle'],
                                    degrees=True).as_matrix()

            rotated_pos = np.dot(BB_T3.atom_pos, R_Matrix) + edge_pos
            self.atom_pos += rotated_pos.tolist()

            self.atom_labels += ['C2' if i == 'C' else i for i in BB_T3.atom_labels]

        StartingFramework = Structure(
            self.lattice,
            self.atom_types,
            self.atom_pos,
            coords_are_cartesian=True,
            site_properties={'source': self.atom_labels}
        ).get_sorted_structure()

        # Translates the structure to the center of the cell
        StartingFramework.translate_sites(
            range(len(StartingFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True
        )

        dict_structure = StartingFramework.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [2/3, 1/3, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'AB2':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [1/2, 0, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [1/2, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'ABC1':
            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (2/3, 1/3, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (4/3, 2/3, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        if stacking == 'ABC2':
            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (1/3, 0, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (2/3, 0, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        if stacking == 'AAl':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        # Create AA tilted stacking.
        if stacking == 'AAt':
            cell = StartingFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            self.lattice = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        dict_structure = stacked_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure, symprec=0.1, angle_tolerance=.5)

        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            print_command(self.prim_structure, self.verbosity, ['debug'])

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()
        except Exception as e:
            print_command(e, self.verbosity, ['debug'])

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        if print_result is True:
            print_framework_name(self.name,
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

        self.name = f'{BB_H6.name}-{BB_L2.name}-HXL_A-{stacking}'
        self.topology = 'HXL_A'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_H6.charge + BB_L2.charge
        self.chirality = BB_H6.chirality or BB_L2.chirality

        print_command(f'Starting the creation of {self.name}', self.verbosity, ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_H6.conector, BB_L2.conector)

        # Replace "X" the building block
        BB_L2.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_H6.remove_X()
        BB_L2.remove_X()

        print_command(f'Bond atom detected: {bond_atom}', self.verbosity, ['debug', 'high'])

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

        if self.stacking == 'A':
            c = slab

        # Create the lattice
        self.lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        self.atom_types = []
        self.atom_labels = []
        self.atom_pos = []

        # Add the building blocks to the structure
        for vertice_data in topology_info['vertices']:
            self.atom_types += BB_H6.atom_types
            vertice_pos = np.array(vertice_data['position'])*a

            R_Matrix = R.from_euler('z',
                                    vertice_data['angle'],
                                    degrees=True).as_matrix()

            rotated_pos = np.dot(BB_H6.atom_pos, R_Matrix) + vertice_pos
            self.atom_pos += rotated_pos.tolist()

            self.atom_labels += ['C1' if i == 'C' else i for i in BB_H6.atom_labels]

        # Add the building blocks to the structure
        for edge_data in topology_info['edges']:
            self.atom_types += BB_L2.atom_types
            edge_pos = np.array(edge_data['position'])*a

            R_Matrix = R.from_euler('z',
                                    edge_data['angle'],
                                    degrees=True).as_matrix()

            rotated_pos = np.dot(BB_L2.atom_pos, R_Matrix) + edge_pos
            self.atom_pos += rotated_pos.tolist()

            self.atom_labels += ['C2' if i == 'C' else i for i in BB_L2.atom_labels]

        StartingFramework = Structure(
            self.lattice,
            self.atom_types,
            self.atom_pos,
            coords_are_cartesian=True,
            site_properties={'source': self.atom_labels}
        ).get_sorted_structure()

        # Translates the structure to the center of the cell
        StartingFramework.translate_sites(
            range(len(StartingFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True
        )

        dict_structure = StartingFramework.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [2/3, 1/3, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'AB2':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [1/2, 0, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [1/2, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'ABC1':
            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (2/3, 1/3, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (4/3, 2/3, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        if stacking == 'ABC2':
            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (1/3, 0, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (2/3, 0, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        if stacking == 'AAl':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        # Create AA tilted stacking.
        if stacking == 'AAt':
            cell = StartingFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            self.lattice = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        dict_structure = stacked_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure, symprec=0.1, angle_tolerance=.5)

        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            print_command(self.prim_structure, self.verbosity, ['debug'])

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()
        except Exception as e:
            print_command(e, self.verbosity, ['debug'])

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        if print_result is True:
            print_framework_name(self.name,
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

    def _create_fxt_structure(self,
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

        self.name = f'{BB_S4_A.name}-{BB_S4_B.name}-SQL-{stacking}'
        self.topology = 'FXT'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_S4_A.charge + BB_S4_B.charge
        self.chirality = BB_S4_A.chirality or BB_S4_B.chirality

        print_command(f'Starting the creation of {self.name}', self.verbosity, ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_S4_A.conector, BB_S4_B.conector)

        # Replace "X" the building block
        BB_S4_A.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_S4_A.remove_X()
        BB_S4_B.remove_X()

        print_command(f'Bond atom detected: {bond_atom}', self.verbosity, ['debug', 'high'])

        # Get the topology information
        topology_info = TOPOLOGY_DICT[self.topology]

        # Measure the base size of the building blocks
        size = 2 * (BB_S4_A.size[0] + BB_S4_B.size[0])

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

        if self.stacking == 'A':
            c = slab

        # Create the lattice
        self.lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        self.atom_types = []
        self.atom_labels = []
        self.atom_pos = []

        # Add the first building block to the structure
        vertice_data = topology_info['vertices'][1]
        self.atom_types += BB_S4_A.atom_types
        vertice_pos = np.array(vertice_data['position'])*a

        R_Matrix = R.from_euler('z', vertice_data['angle'], degrees=True).as_matrix()

        rotated_pos = np.dot(BB_S4_A.atom_pos, R_Matrix) + vertice_pos
        self.atom_pos += rotated_pos.tolist()

        self.atom_labels += ['C1' if i == 'C' else i for i in BB_S4_A.atom_labels]

        # Add the second building block to the structure
        for n in [0, 2]:
            vertice_data = topology_info['vertices'][n]
            self.atom_types += BB_S4_B.atom_types
            vertice_pos = np.array(vertice_data['position'])*a

            R_Matrix = R.from_euler('z', vertice_data['angle'], degrees=True).as_matrix()

            rotated_pos = np.dot(BB_S4_B.atom_pos, R_Matrix) + vertice_pos
            self.atom_pos += rotated_pos.tolist()

            self.atom_labels += ['C2' if i == 'C' else i for i in BB_S4_B.atom_labels]

        StartingFramework = Structure(
            self.lattice,
            self.atom_types,
            self.atom_pos,
            coords_are_cartesian=True,
            site_properties={'source': self.atom_labels}
        ).get_sorted_structure()

        # Translates the structure to the center of the cell
        StartingFramework.translate_sites(
            range(len(StartingFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True
        )

        dict_structure = StartingFramework.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':

            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [1/4, 1/4, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [1/4, 1/4, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'AB2':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [1/2, 0, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [1/2, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'ABC1':

            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (1/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (1/3, 1/3, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 2/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (2/3, 2/3, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        if stacking == 'ABC2':
            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (1/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (1/3, 0, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 2/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (2/3, 0, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        # Create AAl stacking. Tetragonal cell with two sheets
        # per cell shifited by the shift_vector in angstroms.
        if stacking == 'AAl':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        # Create AA tilted stacking.
        # Tilted tetragonal cell with two sheets per cell tilted by tilt_angle.
        if stacking == 'AAt':
            cell = StartingFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            self.lattice = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        dict_structure = stacked_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure, symprec=0.05, angle_tolerance=.5)
        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            print_command(self.prim_structure, self.verbosity, ['debug'])

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()

        except Exception as e:
            print_command(e, self.verbosity, ['debug'])

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        if print_result is True:
            print_framework_name(self.name,
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

        The FXT-A net is composed of one tetrapodal and one linear building blocks.

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

        self.name = f'{BB_S4.name}-{BB_L2.name}-SQL_A-{stacking}'
        self.topology = 'FXT_A'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_S4.charge + BB_L2.charge
        self.chirality = BB_S4.chirality or BB_L2.chirality

        print_command(f'Starting the creation of {self.name}', self.verbosity, ['debug', 'high'])

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_S4.conector, BB_L2.conector)

        # Replace "X" the building block
        BB_L2.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_S4.remove_X()
        BB_L2.remove_X()

        print_command(f'Bond atom detected: {bond_atom}', self.verbosity, ['debug', 'high'])

        # Get the topology information
        topology_info = TOPOLOGY_DICT[self.topology]

        # Measure the base size of the building blocks
        size = 2 * (BB_S4.size[0] + BB_L2.size[0])

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

        if self.stacking == 'A':
            c = slab

        # Create the lattice
        self.lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

        # Create the structure
        self.atom_types = []
        self.atom_labels = []
        self.atom_pos = []

        # Add the building blocks to the structure
        for vertice_data in topology_info['vertices']:
            self.atom_types += BB_S4.atom_types
            vertice_pos = np.array(vertice_data['position'])*a

            R_Matrix = R.from_euler('z', vertice_data['angle'], degrees=True).as_matrix()

            rotated_pos = np.dot(BB_S4.atom_pos, R_Matrix) + vertice_pos
            self.atom_pos += rotated_pos.tolist()

            self.atom_labels += ['C1' if i == 'C' else i for i in BB_S4.atom_labels]

        # Add the building blocks to the structure
        for edge_data in topology_info['edges']:
            self.atom_types += BB_L2.atom_types
            edge_pos = np.array(edge_data['position'])*a

            R_Matrix = R.from_euler('z', edge_data['angle'], degrees=True).as_matrix()

            rotated_pos = np.dot(BB_L2.atom_pos, R_Matrix) + edge_pos
            self.atom_pos += rotated_pos.tolist()

            self.atom_labels += ['C2' if i == 'C' else i for i in BB_L2.atom_labels]

        StartingFramework = Structure(
            self.lattice,
            self.atom_types,
            self.atom_pos,
            coords_are_cartesian=True,
            site_properties={'source': self.atom_labels}
        ).get_sorted_structure()

        # Translates the structure to the center of the cell
        StartingFramework.translate_sites(
            range(len(StartingFramework.as_dict()['sites'])),
            [0, 0, 0.5],
            frac_coords=True,
            to_unit_cell=True
        )

        dict_structure = StartingFramework.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [2/3, 1/3, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'AB2':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [1/2, 0, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [1/2, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        if stacking == 'ABC1':
            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (2/3, 1/3, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (4/3, 2/3, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        if stacking == 'ABC2':
            self.lattice *= (1, 1, 3)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            _, B_list, C_list = np.split(np.arange(len(self.atom_types)), 3)

            # Translate the second sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                (1/3, 0, 1/3),
                frac_coords=True,
                to_unit_cell=True
                )

            # Translate the third sheet by the vector (2/3, 1/3, 0) to generate the B positions
            stacked_structure.translate_sites(
                C_list,
                (2/3, 0, 2/3),
                frac_coords=True,
                to_unit_cell=True
            )

        if stacking == 'AAl':
            self.lattice *= (1, 1, 2)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        # Create AA tilted stacking.
        if stacking == 'AAt':
            cell = StartingFramework.as_dict()['lattice']

            # Shift the cell by the tilt angle
            a_cell = cell['a']
            b_cell = cell['b']
            c_cell = cell['c'] * 2
            alpha = cell['alpha'] - tilt_angle
            beta = cell['beta'] - tilt_angle
            gamma = cell['gamma']

            self.lattice = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.lattice,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

            # Get the index of the atoms in the second sheet
            B_list = np.split(np.arange(len(self.atom_types)), 2)[1]

            # Translate the second sheet by the vector [2/3, 1/3, 0.5] to generate the B positions
            stacked_structure.translate_sites(
                B_list,
                [0, 0, 0.5],
                frac_coords=True,
                to_unit_cell=True
                )

        dict_structure = stacked_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure, symprec=0.05, angle_tolerance=.5)
        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            print_command(self.prim_structure, self.verbosity, ['debug'])

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()

        except Exception as e:
            print_command(e, self.verbosity, ['debug'])

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        if print_result is True:
            print_framework_name(self.name,
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
        bond_atom = get_bond_atom(bb_1.conector, bb_2.conector)

        print_command(f'Bond atom detected: {bond_atom}', self.verbosity, ['debug', 'high'])

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

        R_matrix = rotation_matrix_from_vectors(bb1_lado1, np.array([1, 0, 0]))

        bb_1.atom_pos = np.dot(bb_1.atom_pos, np.transpose(R_matrix))

        # Add tbe building block 1 (C4) on the center [1/2, 1/2, 0.0] of the unitary cell (A1 site)
        final_pos = np.dot(bb_1.atom_pos, R.from_euler(
            'z', 30, degrees=True).as_matrix()) + np.array([1/4, np.sqrt(3)/4, 0])*a
        final_label = list(bb_1.atom_types)

        # Add tbe building block 1 (C4) on [1/2, 0.0, 0.0] of the unitary cell (A2 site) -45
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler(
            'z', -90, degrees=True).as_matrix()) + np.array([1/2, 0, 0])*a))
        final_label += list(bb_1.atom_types)

        # Add tbe building block 1 (C4) on [0.0, 1/2, 0.0] of the unitary cell (A3 site) 30
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler(
            'z', -30, degrees=True).as_matrix()) + np.array([-1/4, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_1.atom_types)

        # Changes the X atoms by the desired bond_atom
        final_label, final_pos = change_X_atoms(final_label, final_pos, bond_atom)

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
        final_label, final_pos = change_X_atoms(final_label, final_pos, bond_atom)

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
            self.prim_structure = struct_symm_prim

        # Create AA staking. By default one sheet per unitary cell is used
        if stacking == 'AA':
            self.stacking = 'AA'
            self.prim_structure = struct_symm_prim

        if stacking == 'AB1':
            self.stacking = 'AB1'
            labels_conv_crystal = np.array(
                [[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()[
                            'lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal*(1, 1, 0.5)
            B = translate_inside(
                ion_conv_crystal*(1, 1, 1.5) + (1/4, 1/4, 0))

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_1 = Structure(lattice, AB_label + AB_label,
                             AB, coords_are_cartesian=False)
            AB_1_symm = SpacegroupAnalyzer(
                AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.prim_structure = AB_1_symm.get_refined_structure()

        if stacking == 'AB2':
            self.stacking = 'AB2'
            labels_conv_crystal = np.array(
                [[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()[
                            'lattice']['matrix'])*(1, 1, 2)

            A = ion_conv_crystal*(1, 1, 0.5)
            B = translate_inside(
                ion_conv_crystal*(1, 1, 1.5) + (1/2, 0, 0))

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            AB_2 = Structure(lattice, AB_label + AB_label,
                             AB, coords_are_cartesian=False)
            AB_2_symm = SpacegroupAnalyzer(
                AB_2, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.prim_structure = AB_2_symm.get_refined_structure()

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
            r = get_cartesian_to_fractional_matrix(
                *cell_to_cellpar(cell))
            shift_vector = np.dot(r, np.array(shift_vector))

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal * np.array([1, 1, 1.5]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)
            self.prim_structure = Structure(
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

            new_cell = cellpar_to_cell(
                [a_cell, b_cell, c_cell, alpha, beta, gamma])

            # Shift the first sheet to be at 0.25 * c
            A = ion_conv_crystal

            # Shift the first sheet to be at 0.75 * c and translate by the shift_vector
            B = ion_conv_crystal + np.array([0, 0, 0.5*c_cell]) + shift_vector

            AB = np.concatenate((A, B))
            AB_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(new_cell)
            self.prim_structure = Structure(
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
            B = translate_inside(
                ion_conv_crystal*(1, 1, 1) + (2/3, 1/3, 0))
            C = translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (4/3, 2/3, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label +
                              ABC_label, ABC, coords_are_cartesian=False)
            ABC_f_symm = SpacegroupAnalyzer(
                ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.prim_structure = ABC_f_symm.get_refined_structure()

        if stacking == 'ABC2':
            self.stacking = 'ABC2'
            labels_conv_crystal = np.array(
                [[i['label']] for i in struct_symm_prim.as_dict()['sites']])
            ion_conv_crystal = np.array(
                [i['abc'] for i in struct_symm_prim.as_dict()['sites']])
            cell = np.array(struct_symm_prim.as_dict()[
                            'lattice']['matrix'])*(1, 1, 3)

            A = ion_conv_crystal*(1, 1, 5/3)
            B = translate_inside(
                ion_conv_crystal*(1, 1, 1) + (1/3, 0, 0))
            C = translate_inside(
                ion_conv_crystal*(1, 1, 1/3) + (2/3, 0, 0))

            ABC = np.concatenate((A, B, C))
            ABC_label = [i[0] for i in labels_conv_crystal]

            lattice = Lattice(cell)

            ABC_f = Structure(lattice, ABC_label+ABC_label +
                              ABC_label, ABC, coords_are_cartesian=False)
            ABC_f_symm = SpacegroupAnalyzer(
                ABC_f, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.prim_structure = ABC_f_symm.get_refined_structure()

        dict_structure = self.prim_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.n_atoms = len(self.prim_structure)
        self.composition = self.prim_structure.formula

        if self.verbosity is True:
            print(self.prim_structure)

        # Get the simmetry information of the generated structure
        self.lattice_type = symm.get_lattice_type()
        self.space_group = symm.get_space_group_symbol()
        self.space_group_n = symm.get_space_group_number()

        symm_op = symm.get_point_group_operations()
        self.hall = symm.get_hall()

        if print_result is True:
            print_framework_name(self.name,
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
        bond_atom = get_bond_atom(bb_1.conector, bb_2.conector)

        print_command(f'Bond atom detected: {bond_atom}', self.verbosity, ['debug', 'high'])

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

        R_matrix = rotation_matrix_from_vectors(
            bb1_lado1, np.array([1, 0, 0]))

        bb_1.atom_pos = np.dot(bb_1.atom_pos, np.transpose(R_matrix))

        # Add tbe building block 1 (C4) on the center [1/2, 1/2, 1/2] of the unitary cell (A1 site)
        final_pos = np.dot(bb_1.atom_pos, R.from_euler(
            'z', 30, degrees=True).as_matrix()) + np.array([1/4, np.sqrt(3)/4, 0])*a
        final_label = list(bb_1.atom_types)

        # Add tbe building block 1 (C4) on [1/2, 0.0, 0.0] of the unitary cell (A2 site)
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler(
            'z', -90, degrees=True).as_matrix()) + np.array([1/2, 0, 0])*a))
        final_label += list(bb_1.atom_types)

        # Add tbe building block 1 (C4) on [0.0, 1/2, 0.0] of the unitary cell (A3 site)
        final_pos = np.vstack((final_pos, np.dot(bb_1.atom_pos, R.from_euler(
            'z', -30, degrees=True).as_matrix()) + np.array([-1/4, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_1.atom_types)

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
        final_label += list(bb_2.atom_types)

        # Add the building block 2 (C2) on [1/4, 3/4, 0.0] of the unitary cell (B2 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler(
            'z', -30, degrees=True).as_matrix()) + np.array([1/8, 3*np.sqrt(3)/8, 0])*a))
        final_label += list(bb_2.atom_types)

        # Add tbe building block 2 (C2) on [3/4, 1/4, 0.0] of the unitary cell (B3 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler(
            'z', 30, degrees=True).as_matrix()) + np.array([5/8, np.sqrt(3)/8, 0])*a))
        final_label += list(bb_2.atom_types)

        # Add tbe building block 2 (C2) on [1/4, 3/4, 0.0] of the unitary cell (B4 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler(
            'z', 30, degrees=True).as_matrix()) + np.array([-1/8, 3*np.sqrt(3)/8, 0])*a))
        final_label += list(bb_2.atom_types)

        # Add tbe building block 2 (C2) on [3/4, 1/2, 0.0] of the unitary cell (B5 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler(
            'z', 90, degrees=True).as_matrix()) + np.array([4/8, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_2.atom_types)

        # Add tbe building block 2 (C2) on [1/2, 3/4, 0.0] of the unitary cell (B6 site)
        final_pos = np.vstack((final_pos, np.dot(bb_2.atom_pos, R.from_euler(
            'z', 90, degrees=True).as_matrix()) + np.array([0.0, np.sqrt(3)/4, 0])*a))
        final_label += list(bb_2.atom_types)

        # Change the X atoms by the desired bond_atom
        final_label, final_pos = change_X_atoms(final_label, final_pos, bond_atom)

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

            self.prim_structure = symm.get_refined_structure()

        # Create AA staking. By default one sheet per unitary cell is used
        if stacking == 'AA':
            self.stacking = 'AA'
            symm = SpacegroupAnalyzer(struct,
                                      symprec=self.symm_tol,
                                      angle_tolerance=self.angle_tol)

            self.prim_structure = symm.get_refined_structure()

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

            self.atom_types = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.prim_structure = symm.get_refined_structure()

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

            self.atom_types = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.prim_structure = symm.get_refined_structure()

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

            self.atom_types = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.prim_structure = symm.get_refined_structure()

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

            self.atom_types = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.prim_structure = symm.get_refined_structure()

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
            r = get_cartesian_to_fractional_matrix(*lattice.parameters)
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

            self.atom_types = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.prim_structure = symm.get_refined_structure()

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

            new_cell = cellpar_to_cell(
                [a_cell, b_cell, c_cell, alpha, beta, gamma])

            lattice = Lattice(new_cell)

            A = struct
            A_labels = np.array([i['label'] for i in A.as_dict()['sites']])
            A_pos = np.array([i['xyz'] for i in A.as_dict()['sites']])

            B = struct
            B_labels = np.array([i['label'] for i in B.as_dict()['sites']])
            B_pos = np.array([i['xyz'] for i in B.as_dict()['sites']])

            # Calculates the shift vector in crystal units
            r = get_cartesian_to_fractional_matrix(*lattice.parameters)
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

            self.atom_types = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]

            symm = SpacegroupAnalyzer(AB_1, symprec=self.symm_tol, angle_tolerance=self.angle_tol)

            self.prim_structure = symm.get_refined_structure()

        dict_structure = self.prim_structure.as_dict()

        self.lattice = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.n_atoms = len(self.prim_structure)
        self.composition = self.prim_structure.formula

        if self.verbosity is True:
            print(self.prim_structure)

        # Get the simmetry information of the generated structure
        self.lattice_type = symm.get_lattice_type()
        self.space_group = symm.get_space_group_symbol()
        self.space_group_n = symm.get_space_group_number()

        symm_op = symm.get_point_group_operations()
        self.hall = symm.get_hall()

        if print_result is True:
            print_framework_name(self.name,
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
