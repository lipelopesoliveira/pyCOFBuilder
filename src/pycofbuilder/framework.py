# -*- coding: utf-8 -*-
# Created by Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

"""
The Framework class implements definitions and methods for a Framework buiding
"""

import os
import copy
import numpy as np

# Import pymatgen
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from scipy.spatial.transform import Rotation as R

# Import pycofbuilder exceptions
from pycofbuilder.exceptions import (BondLenghError,
                                     BBConnectivityError)

# Import pycofbuilder building_block
from pycofbuilder.building_block import BuildingBlock

# Import pycofbuilder topology data
from pycofbuilder.data.topology import TOPOLOGY_DICT

# Import pycofbuilder tools
from pycofbuilder.tools import (get_bond_atom,
                                cell_to_cellpar,
                                cellpar_to_cell,
                                rotation_matrix_from_vectors,
                                unit_vector,
                                angle,
                                get_framework_symm_text,
                                get_bonds)

# Import pycofbuilder io_tools
from pycofbuilder.io_tools import (save_chemjson,
                                   save_cif,
                                   save_xyz,
                                   save_turbomole,
                                   save_vasp,
                                   save_xsf,
                                   save_pdb,
                                   save_pqr,
                                   save_qe,
                                   save_gjf)

from pycofbuilder.logger import create_logger


class Framework():
    """
    A class used to represent a Covalent Organic Framework as a reticular entity.

    ...

    Attributes
    ----------
    name : str
        Name of the material
    out_dir : str
        Path to save the results.
        If not defined, a `out` folder will be created in the current directory.
    verbosity : str
        Control the printing options. Can be 'none', 'normal', or 'debug'.
        Default: 'normal'
    save_bb : bool
        Control the saving of the building blocks.
        Default: True
    lib_path : str
        Path for saving the building block files.
        If not defined, a `building_blocks` folder will be created in the current directory.
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
    dist_threshold : float = 0.8
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
    """

    def __init__(self, name: str = None, **kwargs):

        self.name: str = name

        self.out_path: str = kwargs.get('out_dir', os.path.join(os.getcwd(), 'out'))
        self.save_bb: bool = kwargs.get('save_bb', True)
        self.bb_out_path: str = kwargs.get('bb_out_path', os.path.join(self.out_path, 'building_blocks'))

        self.logger = create_logger(level=kwargs.get('log_level', 'info'),
                                    format=kwargs.get('log_format', 'simple'),
                                    save_to_file=kwargs.get('save_to_file', False),
                                    log_filename=kwargs.get('log_filename', 'pycofbuilder.log'))

        self.symm_tol = kwargs.get('symm_tol', 0.1)
        self.angle_tol = kwargs.get('angle_tol', 0.5)
        self.dist_threshold = kwargs.get('dist_threshold', 0.8)
        self.bond_threshold = kwargs.get('bond_threshold', 1.3)

        self.bb1_name = None
        self.bb2_name = None
        self.topology = None
        self.stacking = None
        self.smiles = None

        self.atom_types = []
        self.atom_pos = []
        self.atom_labels = []
        self.cellMatrix = np.eye(3)
        self.cellParameters = np.array([1, 1, 1, 90, 90, 90]).astype(float)
        self.bonds = []

        self.lattice_sgs = None
        self.space_group = None
        self.space_group_n = None

        self.dimention = None
        self.n_atoms = self.get_n_atoms()
        self.mass = None
        self.composition = None
        self.charge = 0
        self.multiplicity = 1
        self.chirality = False

        self.available_2D_top = ['HCB', 'HCB_A',
                                 'SQL', 'SQL_A',
                                 'KGD',
                                 'HXL', 'HXL_A',
                                 'FXT', 'FXT_A']

        # To add: ['dia', 'bor', 'srs', 'pts', 'ctn', 'rra', 'fcc', 'lon', 'stp', 'acs', 'tbo', 'bcu', 'fjh', 'ceq']
        self.available_3D_top = ['DIA', 'DIA_A', 'BOR']  # Temporary
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
            'DIA': [str(i + 1) for i in range(15)],
            'DIA_A': [str(i + 1) for i in range(15)],
            'BOR': [str(i + 1) for i in range(15)]
        }

        if self.name is not None:
            self.from_name(self.name)

    def __str__(self) -> str:
        return self.as_string()

    def __repr__(self) -> str:
        return f'Framework({self.bb1_name}, {self.bb2_name}, {self.topology}, {self.stacking})'

    def as_string(self) -> str:
        """
        Returns a string with the Framework information.
        """

        fram_str = f'Name: {self.name}\n'

        # Get the formula of the framework

        if self.composition is not None:
            fram_str += f'Full Formula ({self.composition})\n'
            fram_str += f'Reduced Formula: ({self.composition})\n'
        else:
            fram_str += 'Full Formula ()\n'
            fram_str += 'Reduced Formula: \n'

        fram_str += 'abc   :  {:11.6f}  {:11.6f}  {:11.6f}\n'.format(*self.cellParameters[:3])
        fram_str += 'angles:  {:11.6f}  {:11.6f}  {:11.6f}\n'.format(*self.cellParameters[3:])
        fram_str += 'A: {:11.6f}  {:11.6f}   {:11.6f}\n'.format(*self.cellMatrix[0])
        fram_str += 'B: {:11.6f}  {:11.6f}   {:11.6f}\n'.format(*self.cellMatrix[1])
        fram_str += 'C: {:11.6f}  {:11.6f}   {:11.6f}\n'.format(*self.cellMatrix[2])

        fram_str += f'Cartesian Sites ({self.n_atoms})\n'
        fram_str += '  #  Type         a         b         c    label\n'
        fram_str += '---  ----  --------  --------  --------  -------\n'

        for i in range(len(self.atom_types)):
            fram_str += '{:3d}  {:4s}  {:8.5f}  {:8.5f}  {:8.5f}  {:>7}\n'.format(i,
                                                                                  self.atom_types[i],
                                                                                  self.atom_pos[i][0],
                                                                                  self.atom_pos[i][1],
                                                                                  self.atom_pos[i][2],
                                                                                  self.atom_labels[i])

        return fram_str

    def get_n_atoms(self) -> int:
        ''' Returns the number of atoms in the unitary cell'''
        return len(self.atom_types)

    def get_available_topologies(self, dimensionality: str = 'all', print_result: bool = True):
        """
        Get the available topologies implemented in the class.

        Parameters
        ----------

        dimensionality : str, optional
            The dimensionality of the topologies to be printed. Can be 'all', '2D' or '3D'.
            Default: 'all'
        print_result: bool, optional
            If True, the available topologies are printed.

        Returns
        -------
        dimensionality_list: list
            A list with the available topologies.
        """

        dimensionality_error = f'Dimensionality must be one of the following: all, 2D, 3D, not {dimensionality}'
        assert dimensionality in ['all', '2D', '3D'], dimensionality_error

        dimensionality_list = []

        if dimensionality == 'all' or dimensionality == '2D':
            if print_result:
                print('Available 2D Topologies:')
            for i in self.available_2D_top:
                if print_result:
                    print(i.upper())
                dimensionality_list.append(i)

        if dimensionality == 'all' or dimensionality == '3D':
            if print_result:
                print('Available 3D Topologies:')
            for i in self.available_3D_top:
                if print_result:
                    print(i.upper())
                dimensionality_list.append(i)

        return dimensionality_list

    def check_name_concistency(self, FrameworkName) -> tuple[str, str, str, str]:
        """
        Checks if the name is in the correct format and returns a
        tuple with the building blocks names, the net and the stacking.

        In case the name is not in the correct format, an error is raised.

        Parameters
        ----------
        FrameworkName : str, required
            The name of the COF to be created

        Returns
        -------
        tuple[str, str, str, str]
            A tuple with the building blocks names, the net and the stacking.
        """

        string_error = 'FrameworkName must be in the format: BB1_BB2_Net_Stacking'
        assert isinstance(FrameworkName, str), string_error

        name_error = 'FrameworkName must be in the format: BB1_BB2_Net_Stacking'
        assert len(FrameworkName.split('-')) == 4, name_error

        bb1_name, bb2_name, Net, Stacking = FrameworkName.split('-')

        net_error = f'{Net} not in the available list: {self.available_topologies}'
        assert Net in self.available_topologies, net_error

        stacking_error = f'{Stacking} not in the available list: {self.available_stacking[Net]}'
        assert Stacking in self.available_stacking[Net], stacking_error

        return bb1_name, bb2_name, Net, Stacking

    def from_name(self, FrameworkName, **kwargs) -> None:
        """Creates a COF from a given FrameworkName.

        Parameters
        ----------
        FrameworkName : str, required
            The name of the COF to be created

        Returns
        -------
        COF : Framework
            The COF object
        """
        bb1_name, bb2_name, Net, Stacking = self.check_name_concistency(FrameworkName)

        bb1 = BuildingBlock(name=bb1_name, bb_out_path=self.bb_out_path, save_bb=self.save_bb)
        bb2 = BuildingBlock(name=bb2_name, bb_out_path=self.bb_out_path, save_bb=self.save_bb)

        self.from_building_blocks(bb1, bb2, Net, Stacking, **kwargs)

    def from_building_blocks(self,
                             bb1: BuildingBlock,
                             bb2: BuildingBlock,
                             net: str,
                             stacking: str,
                             **kwargs):
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
        COF : Framework
            The COF object
        """
        self.name = f'{bb1.name}-{bb2.name}-{net}-{stacking}'
        self.bb1_name = bb1.name
        self.bb2_name = bb2.name
        self.topology = net
        self.stacking = stacking

        # Check if the BB1 has the smiles attribute
        if hasattr(bb1, 'smiles') and hasattr(bb2, 'smiles'):
            self.smiles = f'{bb1.smiles}.{bb2.smiles}'
        else:
            print('WARNING: The smiles attribute is not available for the building blocks')

        net_build_dict = {
            'HCB': self.create_hcb_structure,
            'HCB_A': self.create_hcb_a_structure,
            'SQL': self.create_sql_structure,
            'SQL_A': self.create_sql_a_structure,
            'KGD': self.create_kgd_structure,
            # 'HXL': self.create_hxl_structure,
            'HXL_A': self.create_hxl_a_structure,
            'FXT': self.create_fxt_structure,
            'FXT_A': self.create_fxt_a_structure,
            'DIA': self.create_dia_structure,
            'DIA_A': self.create_dia_a_structure,
            'BOR': self.create_bor_structure
            }

        result = net_build_dict[net](bb1, bb2, stacking, **kwargs)

        structure = Structure(
            self.cellMatrix,
            self.atom_types,
            self.atom_pos,
            coords_are_cartesian=True,
            site_properties={'source': self.atom_labels}
        )

        self.bonds = get_bonds(structure, self.bond_threshold)

        return result

    def save(self, 
             fmt: str = 'cif',
             supercell: list = [1, 1, 1],
             save_dir=None,
             primitive=False,
             save_bonds=True) -> None:
        '''
        Save the structure in a specif file format.

        Parameters
        ----------
        fmt : str, optional
            The file format to be saved
            Can be `json`, `cif`, `xyz`, `turbomole`, `vasp`, `xsf`, `pdb`, `pqr`, `qe`.
            Default: 'cif'
        supercell : list, optional
            The supercell to be used to save the structure.
            Default: [1,1,1]
        save_dir : str, optional
            The path to save the structure. By default, the structure is saved in a
            `out` folder created in the current directory.
        primitive : bool, optional
            If True, the primitive cell is saved. Otherwise, the conventional cell is saved.
            Default: False
        '''

        save_dict = {
            'cjson': save_chemjson,
            'cif': save_cif,
            'xyz': save_xyz,
            'turbomole': save_turbomole,
            'vasp': save_vasp,
            'xsf': save_xsf,
            'pdb': save_pdb,
            'pqr': save_pqr,
            'qe': save_qe,
            'gjf': save_gjf
         }

        file_format_error = f'Format must be one of the following: {save_dict.keys()}'
        assert fmt in save_dict.keys(), file_format_error

        if primitive:
            structure = self.prim_structure

        else:
            structure = Structure(
                self.cellMatrix,
                self.atom_types,
                self.atom_pos,
                coords_are_cartesian=True,
                site_properties={'source': self.atom_labels}
            )

        final_structure = structure.make_supercell(supercell, in_place=False)

        if save_bonds:
            bonds = get_bonds(final_structure, self.bond_threshold)
        else:
            bonds = []

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
                       atom_pos=atom_pos,
                       bonds=bonds)

# --------------- Net creation methods -------------------------- #

    def create_hcb_structure(self,
                             BB_T3_A,
                             BB_T3_B,
                             stacking: str = 'AA',
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

        connectivity_error = 'Building block {} must present connectivity {} not {}'

        if BB_T3_A.connectivity != 3:
            self.logger.error(connectivity_error.format('A', 3, BB_T3_A.connectivity))
            raise BBConnectivityError(3, BB_T3_A.connectivity)
        if BB_T3_B.connectivity != 3:
            self.logger.error(connectivity_error.format('B', 3, BB_T3_B.connectivity))
            raise BBConnectivityError(3, BB_T3_B.connectivity)

        self.name = f'{BB_T3_A.name}-{BB_T3_B.name}-HCB-{stacking}'
        self.topology = 'HCB'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_T3_A.charge + BB_T3_B.charge
        self.chirality = BB_T3_A.chirality or BB_T3_B.chirality

        self.logger.debug(f'Starting the creation of {self.name}')

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_T3_A.conector, BB_T3_B.conector)

        self.logger.debug('{} detected as bond atom for groups {} and {}'.format(bond_atom,
                                                                                 BB_T3_A.conector,
                                                                                 BB_T3_B.conector))

        # Replace "X" the building block
        BB_T3_A.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_T3_A.remove_X()
        BB_T3_B.remove_X()

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
        self.cellMatrix = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        self.cellParameters = np.array([a, b, c, alpha, beta, gamma]).astype(float)

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
            self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)
        self.cellParameters = cell_to_cellpar(self.cellMatrix)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

            self.cellMatrix = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])
            self.cellParameters = np.array([a_cell, b_cell, c_cell, alpha, beta, gamma]).astype(float)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)
        self.cellParameters = cell_to_cellpar(self.cellMatrix)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        dist_matrix = StartingFramework.distance_matrix

        # Check if there are any atoms closer than 0.8 A
        for i in range(len(dist_matrix)):
            for j in range(i+1, len(dist_matrix)):
                if dist_matrix[i][j] < self.dist_threshold:
                    raise BondLenghError(i, j, dist_matrix[i][j], self.dist_threshold)

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure,
                                  symprec=self.symm_tol,
                                  angle_tolerance=self.angle_tol)

        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            self.logger.debug(self.prim_structure)

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()
        except Exception as e:
            self.logger.exception(e)

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        symm_text = get_framework_symm_text(self.name,
                                            str(self.lattice_type),
                                            str(self.hall[0:2]),
                                            str(self.space_group),
                                            str(self.space_group_n),
                                            len(symm_op))

        self.logger.info(symm_text)

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]

    def create_hcb_a_structure(self,
                               BB_T3: str,
                               BB_L2: str,
                               stacking: str = 'AA',
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

        connectivity_error = 'Building block {} must present connectivity {} not {}'
        if BB_T3.connectivity != 3:
            self.logger.error(connectivity_error.format('A', 3, BB_T3.connectivity))
            raise BBConnectivityError(3, BB_T3.connectivity)
        if BB_L2.connectivity != 2:
            self.logger.error(connectivity_error.format('B', 3, BB_L2.connectivity))
            raise BBConnectivityError(2, BB_L2.connectivity)

        self.name = f'{BB_T3.name}-{BB_L2.name}-HCB_A-{stacking}'
        self.topology = 'HCB_A'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_L2.charge + BB_T3.charge
        self.chirality = BB_L2.chirality or BB_T3.chirality

        self.logger.debug(f'Starting the creation of {self.name}')

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_T3.conector, BB_L2.conector)

        self.logger.debug('{} detected as bond atom for groups {} and {}'.format(bond_atom,
                                                                                 BB_T3.conector,
                                                                                 BB_L2.conector))

        # Replace "X" the building block
        BB_L2.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_T3.remove_X()
        BB_L2.remove_X()

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
        self.cellMatrix = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        self.cellParameters = np.array([a, b, c, alpha, beta, gamma]).astype(float)

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
            self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

            self.cellMatrix = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])
            self.cellParameters = np.array([a_cell, b_cell, c_cell, alpha, beta, gamma]).astype(float)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        dist_matrix = stacked_structure.distance_matrix

        # Check if there are any atoms closer than 0.8 A
        for i in range(len(dist_matrix)):
            for j in range(i+1, len(dist_matrix)):
                if dist_matrix[i][j] < self.dist_threshold:
                    raise BondLenghError(i, j, dist_matrix[i][j], self.dist_threshold)

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure,
                                  symprec=self.symm_tol,
                                  angle_tolerance=self.angle_tol)

        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            self.logger.debug(self.prim_structure)

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()

        except Exception as e:
            self.logger.exception(e)

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        symm_text = get_framework_symm_text(self.name,
                                            str(self.lattice_type),
                                            str(self.hall[0:2]),
                                            str(self.space_group),
                                            str(self.space_group_n),
                                            len(symm_op))

        self.logger.info(symm_text)

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]

    def create_sql_structure(self,
                             BB_S4_A: str,
                             BB_S4_B: str,
                             stacking: str = 'AA',
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

        connectivity_error = 'Building block {} must present connectivity {} not {}'
        if BB_S4_A.connectivity != 4:
            self.logger.error(connectivity_error.format('A', 4, BB_S4_A.connectivity))
            raise BBConnectivityError(4, BB_S4_A.connectivity)
        if BB_S4_B.connectivity != 4:
            self.logger.error(connectivity_error.format('B', 4, BB_S4_B.connectivity))
            raise BBConnectivityError(4, BB_S4_B.connectivity)

        self.name = f'{BB_S4_A.name}-{BB_S4_B.name}-SQL-{stacking}'
        self.topology = 'SQL'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_S4_A.charge + BB_S4_B.charge
        self.chirality = BB_S4_A.chirality or BB_S4_B.chirality

        self.logger.debug(f'Starting the creation of {self.name}')

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_S4_A.conector, BB_S4_B.conector)

        self.logger.debug('{} detected as bond atom for groups {} and {}'.format(bond_atom,
                                                                                 BB_S4_A.conector,
                                                                                 BB_S4_B.conector))

        # Replace "X" the building block
        BB_S4_A.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_S4_A.remove_X()
        BB_S4_B.remove_X()

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
        self.cellMatrix = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        self.cellParameters = np.array([a, b, c, alpha, beta, gamma]).astype(float)

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
            self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':

            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

            self.cellMatrix = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])
            self.cellParameters = np.array([a_cell, b_cell, c_cell, alpha, beta, gamma]).astype(float)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        dist_matrix = stacked_structure.distance_matrix

        # Check if there are any atoms closer than 0.8 A
        for i in range(len(dist_matrix)):
            for j in range(i+1, len(dist_matrix)):
                if dist_matrix[i][j] < self.dist_threshold:
                    raise BondLenghError(i, j, dist_matrix[i][j], self.dist_threshold)

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure,
                                  symprec=self.symm_tol,
                                  angle_tolerance=self.angle_tol)

        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            self.logger.debug(self.prim_structure)

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()

        except Exception as e:
            self.logger.exception(e)

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        symm_text = get_framework_symm_text(self.name,
                                            str(self.lattice_type),
                                            str(self.hall[0:2]),
                                            str(self.space_group),
                                            str(self.space_group_n),
                                            len(symm_op))

        self.logger.info(symm_text)

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]

    def create_sql_a_structure(self,
                               BB_S4: str,
                               BB_L2: str,
                               stacking: str = 'AA',
                               c_parameter_base: float = 3.6,
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

        connectivity_error = 'Building block {} must present connectivity {} not {}'
        if BB_S4.connectivity != 4:
            self.logger.error(connectivity_error.format('A', 4, BB_S4.connectivity))
            raise BBConnectivityError(4, BB_S4.connectivity)
        if BB_L2.connectivity != 2:
            self.logger.error(connectivity_error.format('B', 3, BB_L2.connectivity))
            raise BBConnectivityError(2, BB_L2.connectivity)

        self.name = f'{BB_S4.name}-{BB_L2.name}-SQL_A-{stacking}'
        self.topology = 'SQL_A'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_S4.charge + BB_L2.charge
        self.chirality = BB_S4.chirality or BB_L2.chirality

        self.logger.debug(f'Starting the creation of {self.name}')

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_S4.conector, BB_L2.conector)

        self.logger.debug('{} detected as bond atom for groups {} and {}'.format(bond_atom,
                                                                                 BB_S4.conector,
                                                                                 BB_L2.conector))

        # Replace "X" the building block
        BB_L2.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_S4.remove_X()
        BB_L2.remove_X()

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
        self.cellMatrix = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        self.cellParameters = np.array([a, b, c, alpha, beta, gamma]).astype(float)

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
            self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':

            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

            self.cellMatrix = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])
            self.cellParameters = np.array([a_cell, b_cell, c_cell, alpha, beta, gamma]).astype(float)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        dist_matrix = stacked_structure.distance_matrix

        # Check if there are any atoms closer than 0.8 A
        for i in range(len(dist_matrix)):
            for j in range(i+1, len(dist_matrix)):
                if dist_matrix[i][j] < self.dist_threshold:
                    raise BondLenghError(i, j, dist_matrix[i][j], self.dist_threshold)

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure,
                                  symprec=self.symm_tol,
                                  angle_tolerance=self.angle_tol)

        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            self.logger.debug(self.prim_structure)

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()

        except Exception as e:
            self.logger.exception(e)

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        symm_text = get_framework_symm_text(self.name,
                                            str(self.lattice_type),
                                            str(self.hall[0:2]),
                                            str(self.space_group),
                                            str(self.space_group_n),
                                            len(symm_op))

        self.logger.info(symm_text)

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]

    def create_kgd_structure(self,
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
        connectivity_error = 'Building block {} must present connectivity {} not {}'
        if BB_H6.connectivity != 6:
            self.logger.error(connectivity_error.format('A', 6, BB_H6.connectivity))
            raise BBConnectivityError(6, BB_H6.connectivity)
        if BB_T3.connectivity != 3:
            self.logger.error(connectivity_error.format('B', 3, BB_T3.connectivity))
            raise BBConnectivityError(3, BB_T3.connectivity)

        self.name = f'{BB_H6.name}-{BB_T3.name}-KGD-{stacking}'
        self.topology = 'KGD'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_H6.charge + BB_T3.charge
        self.chirality = BB_H6.chirality or BB_T3.chirality

        self.logger.debug(f'Starting the creation of {self.name}')

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_H6.conector, BB_T3.conector)

        self.logger.debug('{} detected as bond atom for groups {} and {}'.format(bond_atom,
                                                                                 BB_H6.conector,
                                                                                 BB_T3.conector))

        # Replace "X" the building block
        BB_H6.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_H6.remove_X()
        BB_T3.remove_X()

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
        self.cellMatrix = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        self.cellParameters = np.array([a, b, c, alpha, beta, gamma]).astype(float)

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
            self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

            self.cellMatrix = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])
            self.cellParameters = np.array([a_cell, b_cell, c_cell, alpha, beta, gamma]).astype(float)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        dist_matrix = stacked_structure.distance_matrix

        # Check if there are any atoms closer than 0.8 A
        for i in range(len(dist_matrix)):
            for j in range(i+1, len(dist_matrix)):
                if dist_matrix[i][j] < self.dist_threshold:
                    raise BondLenghError(i, j, dist_matrix[i][j], self.dist_threshold)

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure,
                                  symprec=self.symm_tol,
                                  angle_tolerance=self.angle_tol)

        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            self.logger.debug(self.prim_structure)

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()
        except Exception as e:
            self.logger.exception(e)

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        symm_text = get_framework_symm_text(self.name,
                                            str(self.lattice_type),
                                            str(self.hall[0:2]),
                                            str(self.space_group),
                                            str(self.space_group_n),
                                            len(symm_op))

        self.logger.info(symm_text)

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]

    def create_hxl_a_structure(self,
                               BB_H6: str,
                               BB_L2: str,
                               stacking: str = 'AA',
                               print_result: bool = True,
                               slab: float = 10.0,
                               shift_vector: list = [1.0, 1.0, 0],
                               tilt_angle: float = 5.0):
        """Creates a COF with HXL-A network.

        The HXL-A net is composed of one hexapodal and one linear building blocks.

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

        connectivity_error = 'Building block {} must present connectivity {} not {}'
        if BB_H6.connectivity != 6:
            self.logger.error(connectivity_error.format('A', 6, BB_H6.connectivity))
            raise BBConnectivityError(6, BB_H6.connectivity)
        if BB_L2.connectivity != 2:
            self.logger.error(connectivity_error.format('B', 3, BB_L2.connectivity))
            raise BBConnectivityError(2, BB_L2.connectivity)

        self.name = f'{BB_H6.name}-{BB_L2.name}-HXL_A-{stacking}'
        self.topology = 'HXL_A'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_H6.charge + BB_L2.charge
        self.chirality = BB_H6.chirality or BB_L2.chirality

        self.logger.debug(f'Starting the creation of {self.name}')

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_H6.conector, BB_L2.conector)

        self.logger.debug('{} detected as bond atom for groups {} and {}'.format(bond_atom,
                                                                                 BB_H6.conector,
                                                                                 BB_L2.conector))

        # Replace "X" the building block
        BB_L2.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_H6.remove_X()
        BB_L2.remove_X()

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
        self.cellMatrix = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        self.cellParameters = np.array([a, b, c, alpha, beta, gamma]).astype(float)

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
            self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

            self.cellMatrix = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])
            self.cellParameters = np.array([a_cell, b_cell, c_cell, alpha, beta, gamma]).astype(float)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        dist_matrix = stacked_structure.distance_matrix

        # Check if there are any atoms closer than 0.8 A
        for i in range(len(dist_matrix)):
            for j in range(i+1, len(dist_matrix)):
                if dist_matrix[i][j] < self.dist_threshold:
                    raise BondLenghError(i, j, dist_matrix[i][j], self.dist_threshold)

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure,
                                  symprec=self.symm_tol,
                                  angle_tolerance=self.angle_tol)

        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            self.logger.debug(self.prim_structure)

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()
        except Exception as e:
            self.logger.exception(e)

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        symm_text = get_framework_symm_text(self.name,
                                            str(self.lattice_type),
                                            str(self.hall[0:2]),
                                            str(self.space_group),
                                            str(self.space_group_n),
                                            len(symm_op))

        self.logger.info(symm_text)

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]

    def create_fxt_structure(self,
                             BB_S4_A: str,
                             BB_S4_B: str,
                             stacking: str = 'AA',
                             print_result: bool = True,
                             slab: float = 10.0,
                             shift_vector: list = [1.0, 1.0, 0],
                             tilt_angle: float = 5.0):
        """Creates a COF with FXT network.

        The FXT net is composed of two tetrapodal building blocks.

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

        connectivity_error = 'Building block {} must present connectivity {} not {}'
        if BB_S4_A.connectivity != 4:
            self.logger.error(connectivity_error.format('A', 4, BB_S4_A.connectivity))
            raise BBConnectivityError(4, BB_S4_A.connectivity)
        if BB_S4_B.connectivity != 4:
            self.logger.error(connectivity_error.format('B', 4, BB_S4_B.connectivity))
            raise BBConnectivityError(4, BB_S4_B.connectivity)

        self.name = f'{BB_S4_A.name}-{BB_S4_B.name}-FXT-{stacking}'
        self.topology = 'FXT'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_S4_A.charge + BB_S4_B.charge
        self.chirality = BB_S4_A.chirality or BB_S4_B.chirality

        self.logger.debug(f'Starting the creation of {self.name}')

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_S4_A.conector, BB_S4_B.conector)

        self.logger.debug('{} detected as bond atom for groups {} and {}'.format(bond_atom,
                                                                                 BB_S4_A.conector,
                                                                                 BB_S4_B.conector))

        # Replace "X" the building block
        BB_S4_A.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_S4_A.remove_X()
        BB_S4_B.remove_X()

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
        self.cellMatrix = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        self.cellParameters = np.array([a, b, c, alpha, beta, gamma]).astype(float)

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
        for vertice_data in topology_info['vertices'][1:]:
            self.atom_types += BB_S4_B.atom_types
            vertice_pos = np.array(vertice_data['position'])*a

            R_Matrix = R.from_euler('z', vertice_data['angle'], degrees=True).as_matrix()

            rotated_pos = np.dot(BB_S4_B.atom_pos, R_Matrix) + vertice_pos
            self.atom_pos += rotated_pos.tolist()

            self.atom_labels += ['C2' if i == 'C' else i for i in BB_S4_B.atom_labels]

        StartingFramework = Structure(
            self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':

            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

            self.cellMatrix = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])
            self.cellParameters = np.array([a_cell, b_cell, c_cell, alpha, beta, gamma]).astype(float)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        dist_matrix = stacked_structure.distance_matrix

        # Check if there are any atoms closer than 0.8 A
        for i in range(len(dist_matrix)):
            for j in range(i+1, len(dist_matrix)):
                if dist_matrix[i][j] < self.dist_threshold:
                    raise BondLenghError(i, j, dist_matrix[i][j], self.dist_threshold)

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure,
                                  symprec=self.symm_tol,
                                  angle_tolerance=self.angle_tol)

        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            self.logger.debug(self.prim_structure)

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()

        except Exception as e:
            self.logger.exception(e)

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        symm_text = get_framework_symm_text(self.name,
                                            str(self.lattice_type),
                                            str(self.hall[0:2]),
                                            str(self.space_group),
                                            str(self.space_group_n),
                                            len(symm_op))

        self.logger.info(symm_text)

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]

    def create_fxt_a_structure(self,
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

        connectivity_error = 'Building block {} must present connectivity {} not {}'
        if BB_S4.connectivity != 4:
            self.logger.error(connectivity_error.format('A', 4, BB_S4.connectivity))
            raise BBConnectivityError(4, BB_S4.connectivity)
        if BB_L2.connectivity != 2:
            self.logger.error(connectivity_error.format('B', 3, BB_L2.connectivity))
            raise BBConnectivityError(2, BB_L2.connectivity)

        self.name = f'{BB_S4.name}-{BB_L2.name}-FXT_A-{stacking}'
        self.topology = 'FXT_A'
        self.staking = stacking
        self.dimension = 2

        self.charge = BB_S4.charge + BB_L2.charge
        self.chirality = BB_S4.chirality or BB_L2.chirality

        self.logger.debug(f'Starting the creation of {self.name}')

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_S4.conector, BB_L2.conector)

        self.logger.debug('{} detected as bond atom for groups {} and {}'.format(bond_atom,
                                                                                 BB_S4.conector,
                                                                                 BB_L2.conector))

        # Replace "X" the building block
        BB_L2.replace_X(bond_atom)

        # Remove the "X" atoms from the the building block
        BB_S4.remove_X()
        BB_L2.remove_X()

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
        self.cellMatrix = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
        self.cellParameters = np.array([a, b, c, alpha, beta, gamma]).astype(float)

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
            self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]

        if stacking == 'A' or stacking == 'AA':
            stacked_structure = StartingFramework

        if stacking == 'AB1':
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 3)
            self.cellParameters *= (1, 1, 3, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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
            self.cellMatrix *= (1, 1, 2)
            self.cellParameters *= (1, 1, 2, 1, 1, 1)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            sv = np.array(shift_vector)
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos + sv))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

            self.cellMatrix = cellpar_to_cell([a_cell, b_cell, c_cell, alpha, beta, gamma])
            self.cellParameters = np.array([a_cell, b_cell, c_cell, alpha, beta, gamma]).astype(float)

            self.atom_types = np.concatenate((self.atom_types, self.atom_types))
            self.atom_pos = np.concatenate((self.atom_pos, self.atom_pos))
            self.atom_labels = np.concatenate((self.atom_labels, self.atom_labels))

            stacked_structure = Structure(
                self.cellMatrix,
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

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = stacked_structure.formula

        dist_matrix = stacked_structure.distance_matrix

        # Check if there are any atoms closer than 0.8 A
        for i in range(len(dist_matrix)):
            for j in range(i+1, len(dist_matrix)):
                if dist_matrix[i][j] < self.dist_threshold:
                    raise BondLenghError(i, j, dist_matrix[i][j], self.dist_threshold)

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(stacked_structure,
                                  symprec=self.symm_tol,
                                  angle_tolerance=self.angle_tol)

        try:
            self.prim_structure = symm.get_refined_structure(keep_site_properties=True)

            self.logger.debug(self.prim_structure)

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()

        except Exception as e:
            self.logger.exception(e)

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        symm_text = get_framework_symm_text(self.name,
                                            str(self.lattice_type),
                                            str(self.hall[0:2]),
                                            str(self.space_group),
                                            str(self.space_group_n),
                                            len(symm_op))

        self.logger.info(symm_text)

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]

    def create_dia_structure(self,
                             BB_D41: str,
                             BB_D42: str,
                             interp_dg: str = '1',
                             d_param_base: float = 7.2,
                             print_result: bool = True,
                             **kwargs):
        """Creates a COF with DIA network.

        The DIA net is composed of two tetrapodal tetrahedical building blocks.

        Parameters
        ----------
        BB_D41 : BuildingBlock, required
            The BuildingBlock object of the tetrapodal tetrahedical Buiding Block
        BB_D42 : BuildingBlock, required
            The BuildingBlock object of the tetrapodal tetrahedical Buiding Block
        interp_dg : str, optional
            The degree of interpenetration of the framework (default is '1')
        d_param_base : float, optional
            The base value for interlayer distance in angstroms (default is 7.2)
        print_result : bool, optional
            Parameter for the control for printing the result (default is True)

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
        connectivity_error = 'Building block {} must present connectivity {} not {}'
        if BB_D41.connectivity != 4:
            self.logger.error(connectivity_error.format('A', 4, BB_D41.connectivity))
            raise BBConnectivityError(4, BB_D41.connectivity)
        if BB_D42.connectivity != 4:
            self.logger.error(connectivity_error.format('B', 4, BB_D42.connectivity))
            raise BBConnectivityError(4, BB_D42.connectivity)

        self.name = f'{BB_D41.name}-{BB_D42.name}-DIA-{interp_dg}'
        self.topology = 'DIA'
        self.staking = interp_dg
        self.dimension = 3

        self.charge = BB_D41.charge + BB_D42.charge
        self.chirality = BB_D41.chirality or BB_D42.chirality

        self.logger.debug(f'Starting the creation of {self.name}')

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_D41.conector, BB_D42.conector)

        self.logger.debug('{} detected as bond atom for groups {} and {}'.format(bond_atom,
                                                                                 BB_D41.conector,
                                                                                 BB_D42.conector))

        # Get the topology information
        topology_info = TOPOLOGY_DICT[self.topology]

        # Measure the base size of the building blocks
        size = np.average(BB_D41.size) + np.average(BB_D42.size)

        # Calculate the primitive cell vector assuming tetrahedical building blocks
        a_prim = np.sqrt(2)*size*np.sqrt((1 - np.cos(1.9106316646041868)))
        a_conv = np.sqrt(2)*a_prim

        # Create the primitive lattice
        self.cellMatrix = Lattice(a_conv/2*np.array(topology_info['lattice']))
        self.cellParameters = np.array([a_prim, a_prim, a_prim, 60, 60, 60]).astype(float)

        # Create the structure
        self.atom_types = []
        self.atom_labels = []
        self.atom_pos = []

        # Align and rotate the building block 1 to their respective positions
        BB_D41.align_to(topology_info['vertices'][0]['align_v'])

        # Determine the angle that alings the X[1] to one of the vertices of the tetrahedron
        vertice_pos = unit_vector(np.array([1, 0, 1]))
        Q_vertice_pos = BB_D41.get_X_points()[1][1]

        rotated_list = [
             R.from_rotvec(
                 angle * unit_vector(topology_info['vertices'][0]['align_v']), degrees=False
                 ).apply(Q_vertice_pos)
             for angle in np.linspace(0, 2*np.pi, 360)
             ]

        # Calculate the angle between the vertice_pos and the elements of rotated_list
        angle_list = [angle(vertice_pos, i) for i in rotated_list]

        rot_angle = np.linspace(0, 360, 360)[np.argmax(angle_list)]

        BB_D41.rotate_around(rotation_axis=np.array(topology_info['vertices'][0]['align_v']),
                             angle=rot_angle,
                             degree=True)

        BB_D41.shift(np.array(topology_info['vertices'][0]['position'])*a_conv)
        BB_D41.remove_X()

        # Add the building block 1 to the structure
        self.atom_types += BB_D41.atom_types
        self.atom_pos += BB_D41.atom_pos.tolist()
        self.atom_labels += ['C1' if i == 'C' else i for i in BB_D41.atom_labels]

        # Align and rotate the building block 1 to their respective positions
        BB_D42.align_to(topology_info['vertices'][0]['align_v'])

        # Determine the angle that alings the X[1] to one of the vertices of the tetrahedron
        vertice_pos = unit_vector(np.array([1, 0, 1]))
        Q_vertice_pos = BB_D42.get_X_points()[1][1]

        rotated_list = [
             R.from_rotvec(
                 angle * unit_vector(topology_info['vertices'][0]['align_v']), degrees=False
                 ).apply(Q_vertice_pos)
             for angle in np.linspace(0, 2*np.pi, 360)
             ]

        # Calculate the angle between the vertice_pos and the elements of rotated_list
        angle_list = [angle(vertice_pos, i) for i in rotated_list]

        rot_angle = np.linspace(0, 360, 360)[np.argmax(angle_list)]

        BB_D42.rotate_around(rotation_axis=np.array(topology_info['vertices'][0]['align_v']),
                             angle=rot_angle,
                             degree=True)

        BB_D42.atom_pos = -BB_D42.atom_pos

        BB_D42.shift(np.array(topology_info['vertices'][1]['position'])*a_conv)

        BB_D42.replace_X(bond_atom)
        BB_D42.remove_X()

        # Add the building block 2 to the structure
        self.atom_types += BB_D42.atom_types
        self.atom_pos += BB_D42.atom_pos.tolist()
        self.atom_labels += ['C2' if i == 'C' else i for i in BB_D42.atom_labels]

        atom_types, atom_labels, atom_pos = [], [], []
        for n_int in range(int(self.stacking)):
            int_direction = np.array([0, 1, 0]) * d_param_base * n_int

            atom_types += self.atom_types
            atom_pos += (np.array(self.atom_pos) + int_direction).tolist()
            atom_labels += self.atom_labels

        self.atom_types = atom_types
        self.atom_pos = atom_pos
        self.atom_labels = atom_labels

        StartingFramework = Structure(
            self.cellMatrix,
            self.atom_types,
            self.atom_pos,
            coords_are_cartesian=True,
            site_properties={'source': self.atom_labels}
        ).get_sorted_structure()

        StartingFramework.to(os.path.join(os.getcwd(), 'TESTE_DIA.cif'), fmt='cif')

        dict_structure = StartingFramework.as_dict()

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = StartingFramework.formula

        dist_matrix = StartingFramework.distance_matrix

        # Check if there are any atoms closer than 0.8 A
        for i in range(len(dist_matrix)):
            for j in range(i+1, len(dist_matrix)):
                if dist_matrix[i][j] < self.dist_threshold:
                    raise BondLenghError(i, j, dist_matrix[i][j], self.dist_threshold)

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(StartingFramework,
                                  symprec=self.symm_tol,
                                  angle_tolerance=self.angle_tol)

        try:
            self.prim_structure = symm.get_primitive_standard_structure(keep_site_properties=True)

            dict_structure = symm.get_refined_structure(keep_site_properties=True).as_dict()

            self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)
            self.cellParameters = np.array([dict_structure['lattice']['a'],
                                            dict_structure['lattice']['b'],
                                            dict_structure['lattice']['c'],
                                            dict_structure['lattice']['alpha'],
                                            dict_structure['lattice']['beta'],
                                            dict_structure['lattice']['gamma']]).astype(float)

            self.atom_types = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
            self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
            self.n_atoms = len(dict_structure['sites'])
            self.composition = self.prim_structure.formula

            self.logger.debug(self.prim_structure)

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()

        except Exception as e:
            self.logger.exception(e)

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        symm_text = get_framework_symm_text(self.name,
                                            str(self.lattice_type),
                                            str(self.hall[0:2]),
                                            str(self.space_group),
                                            str(self.space_group_n),
                                            len(symm_op))

        self.logger.info(symm_text)

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]

    def create_dia_a_structure(self,
                               BB_D4: str,
                               BB_L2: str,
                               interp_dg: str = '1',
                               d_param_base: float = 7.2,
                               print_result: bool = True,
                               **kwargs):
        """Creates a COF with DIA-A network.

        The DIA net is composed of two tetrapodal tetrahedical building blocks.

        Parameters
        ----------
        BB_D4 : BuildingBlock, required
            The BuildingBlock object of the tetrapodal tetrahedical Buiding Block
        BB_L2 : BuildingBlock, required
            The BuildingBlock object of the dipodal linear Buiding Block
        interp_dg : str, optional
            The degree of interpenetration of the framework (default is '1')
        d_param_base : float, optional
            The base value for interlayer distance in angstroms (default is 7.2)
        print_result : bool, optional
            Parameter for the control for printing the result (default is True)

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
        connectivity_error = 'Building block {} must present connectivity {} not {}'
        if BB_D4.connectivity != 4:
            self.logger.error(connectivity_error.format('A', 4, BB_D4.connectivity))
            raise BBConnectivityError(4, BB_D4.connectivity)
        if BB_L2.connectivity != 2:
            self.logger.error(connectivity_error.format('B', 2, BB_L2.connectivity))
            raise BBConnectivityError(2, BB_L2.connectivity)

        self.name = f'{BB_D4.name}-{BB_L2.name}-DIA_A-{interp_dg}'
        self.topology = 'DIA_A'
        self.staking = interp_dg
        self.dimension = 3

        self.charge = BB_D4.charge + BB_L2.charge
        self.chirality = BB_D4.chirality or BB_L2.chirality

        self.logger.debug(f'Starting the creation of {self.name}')

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_D4.conector, BB_L2.conector)

        self.logger.debug('{} detected as bond atom for groups {} and {}'.format(bond_atom,
                                                                                 BB_D4.conector,
                                                                                 BB_L2.conector))

        # Get the topology information
        topology_info = TOPOLOGY_DICT[self.topology]

        # Measure the base size of the building blocks
        size = 2 * (np.average(BB_D4.size) + np.average(BB_L2.size))

        # Calculate the primitive cell vector assuming tetrahedical building blocks
        a_prim = np.sqrt(2)*size*np.sqrt((1 - np.cos(1.9106316646041868)))
        a_conv = np.sqrt(2)*a_prim

        # Create the primitive lattice
        self.cellMatrix = Lattice(a_conv/2*np.array(topology_info['lattice']))
        self.cellParameters = np.array([a_prim, a_prim, a_prim, 60, 60, 60]).astype(float)

        # Create the structure
        self.atom_types = []
        self.atom_labels = []
        self.atom_pos = []

        # Align and rotate the building block 1 to their respective positions
        BB_D4.align_to(topology_info['vertices'][0]['align_v'])

        # Determine the angle that alings the X[1] to one of the vertices of the tetrahedron
        vertice_pos = unit_vector(np.array([1, 0, 1]))
        Q_vertice_pos = BB_D4.get_X_points()[1][1]

        rotated_list = [
             R.from_rotvec(
                 angle * unit_vector(topology_info['vertices'][0]['align_v']), degrees=False
                 ).apply(Q_vertice_pos)
             for angle in np.linspace(0, 2*np.pi, 360)
             ]

        # Calculate the angle between the vertice_pos and the elements of rotated_list
        angle_list = [angle(vertice_pos, i) for i in rotated_list]

        rot_angle = np.linspace(0, 360, 360)[np.argmax(angle_list)]

        BB_D4.rotate_around(rotation_axis=np.array(topology_info['vertices'][0]['align_v']),
                            angle=rot_angle,
                            degree=True)

        BB_D4.shift(np.array(topology_info['vertices'][0]['position'])*a_conv)
        BB_D4.remove_X()

        # Add the building block 1 to the structure
        self.atom_types += BB_D4.atom_types
        self.atom_pos += BB_D4.atom_pos.tolist()
        self.atom_labels += ['C1' if i == 'C' else i for i in BB_D4.atom_labels]

        # Add the building block 1 to the structure
        self.atom_types += BB_D4.atom_types
        self.atom_pos += list(-np.array(BB_D4.atom_pos) + np.array(topology_info['vertices'][1]['position'])*a_conv)
        self.atom_labels += ['C1' if i == 'C' else i for i in BB_D4.atom_labels]

        # Add the building blocks to the structure
        for edge_data in topology_info['edges']:
            # Copy the building block 2 object
            BB = copy.deepcopy(BB_L2)

            # Align, rotate and shift the building block 2 to their respective positions
            BB.align_to(edge_data['align_v'])
            BB.rotate_around(rotation_axis=edge_data['align_v'],
                             angle=edge_data['angle'])
            BB.shift(np.array(edge_data['position']) * a_conv)

            # Replace "X" the building block with the correct atom dicated by the connection group
            BB.replace_X(bond_atom)
            BB.remove_X()

            # Update the structure
            self.atom_types += BB.atom_types
            self.atom_pos += BB.atom_pos.tolist()
            self.atom_labels += ['C2' if i == 'C' else i for i in BB.atom_labels]

        atom_types, atom_labels, atom_pos = [], [], []
        for n_int in range(int(self.stacking)):
            int_direction = np.array([0, 1, 0]) * d_param_base * n_int

            atom_types += self.atom_types
            atom_pos += (np.array(self.atom_pos) + int_direction).tolist()
            atom_labels += self.atom_labels

        self.atom_types = atom_types
        self.atom_pos = atom_pos
        self.atom_labels = atom_labels

        StartingFramework = Structure(
            self.cellMatrix,
            self.atom_types,
            self.atom_pos,
            coords_are_cartesian=True,
            site_properties={'source': self.atom_labels}
        ).get_sorted_structure()

        StartingFramework.translate_sites(
            np.ones(len(self.atom_types)).astype(int).tolist(),
            [0, 0, 0],
            frac_coords=True,
            to_unit_cell=True
            )

        dict_structure = StartingFramework.as_dict()

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)
        self.cellParameters = np.array([dict_structure['lattice']['a'],
                                        dict_structure['lattice']['b'],
                                        dict_structure['lattice']['c'],
                                        dict_structure['lattice']['alpha'],
                                        dict_structure['lattice']['beta'],
                                        dict_structure['lattice']['gamma']]).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = StartingFramework.formula

        StartingFramework.to(os.path.join(os.getcwd(), 'TESTE_DIA-A.cif'), fmt='cif')

        dist_matrix = StartingFramework.distance_matrix

        # Check if there are any atoms closer than 0.8 A
        for i in range(len(dist_matrix)):
            for j in range(i+1, len(dist_matrix)):
                if dist_matrix[i][j] < self.dist_threshold:
                    raise BondLenghError(i, j, dist_matrix[i][j], self.dist_threshold)

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(StartingFramework,
                                  symprec=self.symm_tol,
                                  angle_tolerance=self.angle_tol)

        try:
            self.prim_structure = symm.get_primitive_standard_structure(keep_site_properties=True)

            dict_structure = symm.get_refined_structure(keep_site_properties=True).as_dict()

            self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)
            self.cellParameters = np.array([dict_structure['lattice']['a'],
                                            dict_structure['lattice']['b'],
                                            dict_structure['lattice']['c'],
                                            dict_structure['lattice']['alpha'],
                                            dict_structure['lattice']['beta'],
                                            dict_structure['lattice']['gamma']]).astype(float)

            self.atom_types = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
            self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
            self.n_atoms = len(dict_structure['sites'])
            self.composition = self.prim_structure.formula

            self.logger.debug(self.prim_structure)

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()

        except Exception as e:
            self.logger.exception(e)

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        symm_text = get_framework_symm_text(self.name,
                                            str(self.lattice_type),
                                            str(self.hall[0:2]),
                                            str(self.space_group),
                                            str(self.space_group_n),
                                            len(symm_op))

        self.logger.info(symm_text)

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]

    def create_bor_structure(self,
                             BB_D4: str,
                             BB_T3: str,
                             interp_dg: str = '1',
                             d_param_base: float = 7.2,
                             print_result: bool = True,
                             **kwargs):
        """Creates a COF with BOR network.

        The DIA net is composed of one tetrapodal tetrahedical building block and
        one tripodal triangular building block.

        Parameters
        ----------
        BB_D4 : BuildingBlock, required
            The BuildingBlock object of the tetrapodal tetrahedical Buiding Block
        BB_T3 : BuildingBlock, required
            The BuildingBlock object of the tripodal triangular Buiding Block
        interp_dg : str, optional
            The degree of interpenetration of the framework (default is '1')
        d_param_base : float, optional
            The base value for interlayer distance in angstroms (default is 7.2)
        print_result : bool, optional
            Parameter for the control for printing the result (default is True)

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
        connectivity_error = 'Building block {} must present connectivity {} not {}'
        if BB_D4.connectivity != 4:
            self.logger.error(connectivity_error.format('A', 4, BB_D4.connectivity))
            raise BBConnectivityError(4, BB_D4.connectivity)
        if BB_T3.connectivity != 3:
            self.logger.error(connectivity_error.format('B', 3, BB_T3.connectivity))
            raise BBConnectivityError(3, BB_T3.connectivity)

        # Get the topology information
        topology_info = TOPOLOGY_DICT[self.topology]

        self.name = f'{BB_D4.name}-{BB_T3.name}-BOR-{interp_dg}'
        self.topology = 'BOR'
        self.staking = interp_dg
        self.dimension = 3

        self.charge = BB_D4.charge + BB_T3.charge
        self.chirality = BB_D4.chirality or BB_T3.chirality

        self.logger.debug(f'Starting the creation of {self.name}')

        # Detect the bond atom from the connection groups type
        bond_atom = get_bond_atom(BB_D4.conector, BB_T3.conector)

        self.logger.debug('{} detected as bond atom for groups {} and {}'.format(bond_atom,
                                                                                 BB_D4.conector,
                                                                                 BB_T3.conector))

        # Get the topology information
        topology_info = TOPOLOGY_DICT[self.topology]

        # Measure the base size of the building blocks
        d_size = (np.array(BB_D4.size).mean() + np.array(BB_T3.size).mean())

        # Calculate the primitive cell vector assuming tetrahedical building blocks
        a_conv = np.sqrt(6) * d_size

        # Create the primitive lattice
        self.cellMatrix = Lattice(a_conv * np.array(topology_info['lattice']))
        self.cellParameters = np.array([a_conv, a_conv, a_conv, 90, 90, 90]).astype(float)

        # Create the structure
        atom_types = []
        atom_labels = []
        atom_pos = []

        for D_site in topology_info['vertices']:
            D4 = BB_D4.copy()
            D4.align_to(
                np.array(D_site['align_v'])
                )

            D4.rotate_around(
                rotation_axis=D_site['align_v'],
                angle=D_site['angle'])

            D4.shift(np.array(D_site['position'])*a_conv)

            atom_types += D4.atom_types
            atom_pos += D4.atom_pos.tolist()
            atom_labels += D4.atom_labels.tolist()

        # Translate all atoms to inside the cell
        for i, pos in enumerate(atom_pos):
            for j, coord in enumerate(pos):
                if coord < 0:
                    atom_pos[i][j] += a_conv

        X_pos = [atom_pos[i] for i in np.where(np.array(atom_types) == 'X')[0]]

        T_site = topology_info['edges'][0]

        _, X = BB_T3.get_X_points()
        BB_T3.rotate_around([0, 0, 1], T_site['angle'], True)

        R_matrix = rotation_matrix_from_vectors([0, 0, 1],
                                                T_site['align_v'])

        BB_T3.atom_pos = np.dot(BB_T3.atom_pos, R_matrix.T)

        BB_T3.replace_X('O')

        # Get the 3 atoms that are closer to T_site['position'])*a_conv
        X_pos_temp = sorted(X_pos, key=lambda x: np.linalg.norm(x - np.array(T_site['position'])*a_conv))

        X_center = np.array(X_pos_temp[:3]).mean(axis=0)

        BB_T3.shift(X_center)

        atom_types += BB_T3.atom_types
        atom_pos += BB_T3.atom_pos.tolist()
        atom_labels += BB_T3.atom_labels.tolist()

        T4 = BB_T3.copy()
        T4.rotate_around([0, 0, 1], 180, True)

        atom_types += T4.atom_types
        atom_pos += T4.atom_pos.tolist()
        atom_labels += T4.atom_labels.tolist()

        T2 = BB_T3.copy()
        T2.rotate_around([0, 0, 1], 90, True)
        T2.rotate_around([1, 0, 0], -90, True)

        atom_types += T2.atom_types
        atom_pos += T2.atom_pos.tolist()
        atom_labels += T2.atom_labels.tolist()

        T3 = BB_T3.copy()
        T3.rotate_around([0, 0, 1], -90, True)

        T3.atom_pos *= np.array([1, 1, -1])

        atom_types += T3.atom_types
        atom_pos += T3.atom_pos.tolist()
        atom_labels += T3.atom_labels.tolist()

        # Translate all atoms to inside the cell
        for i, pos in enumerate(atom_pos):
            for j, coord in enumerate(pos):
                if coord < 0:
                    atom_pos[i][j] += a_conv

        # Remove the X atoms from the list
        X_index = np.where(np.array(atom_types) == 'X')[0]

        self.atom_types = [atom_types[i] for i in range(len(atom_types)) if i not in X_index]
        self.atom_pos = [atom_pos[i] for i in range(len(atom_pos)) if i not in X_index]
        self.atom_labels = [atom_labels[i] for i in range(len(atom_labels)) if i not in X_index]

        atom_types, atom_labels, atom_pos = [], [], []
        for n_int in range(int(self.stacking)):
            int_direction = np.array([0, 1, 0]) * d_param_base * n_int

            atom_types += self.atom_types
            atom_pos += (np.array(self.atom_pos) + int_direction).tolist()
            atom_labels += self.atom_labels

        self.atom_types = atom_types
        self.atom_pos = atom_pos
        self.atom_labels = atom_labels

        StartingFramework = Structure(
            self.cellMatrix,
            self.atom_types,
            self.atom_pos,
            coords_are_cartesian=True,
            site_properties={'source': self.atom_labels}
        ).get_sorted_structure()

        StartingFramework.translate_sites(
            np.ones(len(self.atom_types)).astype(int).tolist(),
            [0, 0, 0],
            frac_coords=True,
            to_unit_cell=True
            )

        dict_structure = StartingFramework.as_dict()

        self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)
        self.cellParameters = np.array([dict_structure['lattice']['a'],
                                        dict_structure['lattice']['b'],
                                        dict_structure['lattice']['c'],
                                        dict_structure['lattice']['alpha'],
                                        dict_structure['lattice']['beta'],
                                        dict_structure['lattice']['gamma']]).astype(float)

        self.atom_types = [i['label'] for i in dict_structure['sites']]
        self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
        self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
        self.n_atoms = len(dict_structure['sites'])
        self.composition = StartingFramework.formula

        StartingFramework.to('TESTE_BOR.cif', fmt='cif')

        dist_matrix = StartingFramework.distance_matrix

        # Check if there are any atoms closer than 0.8 A
        for i in range(len(dist_matrix)):
            for j in range(i+1, len(dist_matrix)):
                if dist_matrix[i][j] < self.dist_threshold:
                    raise BondLenghError(i, j, dist_matrix[i][j], self.dist_threshold)

        # Get the simmetry information of the generated structure
        symm = SpacegroupAnalyzer(StartingFramework,
                                  symprec=self.symm_tol,
                                  angle_tolerance=self.angle_tol)

        try:
            self.prim_structure = symm.get_primitive_standard_structure(keep_site_properties=True)

            dict_structure = symm.get_refined_structure(keep_site_properties=True).as_dict()

            self.cellMatrix = np.array(dict_structure['lattice']['matrix']).astype(float)
            self.cellParameters = np.array([dict_structure['lattice']['a'],
                                            dict_structure['lattice']['b'],
                                            dict_structure['lattice']['c'],
                                            dict_structure['lattice']['alpha'],
                                            dict_structure['lattice']['beta'],
                                            dict_structure['lattice']['gamma']]).astype(float)

            self.atom_types = [i['label'] for i in dict_structure['sites']]
            self.atom_pos = [i['xyz'] for i in dict_structure['sites']]
            self.atom_labels = [i['properties']['source'] for i in dict_structure['sites']]
            self.n_atoms = len(dict_structure['sites'])
            self.composition = self.prim_structure.formula

            self.logger.debug(self.prim_structure)

            self.lattice_type = symm.get_lattice_type()
            self.space_group = symm.get_space_group_symbol()
            self.space_group_n = symm.get_space_group_number()

            symm_op = symm.get_point_group_operations()
            self.hall = symm.get_hall()

        except Exception as e:
            self.logger.exception(e)

            self.lattice_type = 'Triclinic'
            self.space_group = 'P1'
            self.space_group_n = '1'

            symm_op = [1]
            self.hall = 'P 1'

        symm_text = get_framework_symm_text(self.name,
                                            str(self.lattice_type),
                                            str(self.hall[0:2]),
                                            str(self.space_group),
                                            str(self.space_group_n),
                                            len(symm_op))

        self.logger.info(symm_text)

        return [self.name,
                str(self.lattice_type),
                str(self.hall[0:2]),
                str(self.space_group),
                str(self.space_group_n),
                len(symm_op)]
