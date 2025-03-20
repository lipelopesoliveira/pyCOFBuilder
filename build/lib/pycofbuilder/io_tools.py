# -*- coding: utf-8 -*-
# Created by Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

"""
This module contains tools for input and output file manipulation used by pyCOFBuilder.
"""

import os

import warnings
import gemmi
import json

from datetime import date
import numpy as np

from ase.io import read
from pymatgen.io.cif import CifParser


from pycofbuilder.tools import (elements_dict,
                                cell_to_cellpar,
                                cellpar_to_cell,
                                get_fractional_to_cartesian_matrix,
                                get_cartesian_to_fractional_matrix,
                                get_kgrid,
                                smiles_to_xsmiles,
                                cell_to_ibrav)

from pycofbuilder.cjson import ChemJSON


def read_xyz(path: str, file_name: str, extxyz=False) -> tuple:
    """
    Reads a file in format `.xyz` from the given `path` and returns
    a list containg the N atom labels and a Nx3 array contaning
    the atoms coordinates.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the `xyz` file. Does not neet to contain the `.xyz` extention.
    extxyz : bool
        If True, the function will consider the extended xyz file format and use ase library to read the file.

    Returns
    -------
    atonTypes : list
        List of strings containing containg the N atom labels.
    cartPos : numpy array
        Nx3 array contaning the atoms coordinates
    cellMatrix : numpy array
        3x3 array contaning the cell vectors.
    """

    # Remove the extention if exists
    file_name = file_name.split('.')[0]

    # Check if the file exists
    if not os.path.exists(os.path.join(path, file_name + '.xyz')):
        raise FileNotFoundError(f'File {file_name} not found!')

    if extxyz:
        atoms = read(os.path.join(path, file_name + '.xyz'))

        atomTypes = atoms.get_chemical_symbols()  # type: ignore
        cartPos = atoms.get_positions()  # type: ignore
        cellMatrix = atoms.get_cell()  # type: ignore

    else:
        temp_file = open(os.path.join(path, file_name + '.xyz'), 'r').readlines()

        atoms = [i.split() for i in temp_file[2:]]

        atomTypes = [i[0] for i in atoms if len(i) > 1]
        cartPos = np.array([[float(i[1]), float(i[2]), float(i[3])] for i in atoms if len(i) > 1])
        cellMatrix = np.zeros((3, 3))

    return atomTypes, cartPos, cellMatrix


def read_pdb(path, file_name):
    """
    Reads a file in format `.pdb` from the `path` given and returns
    a list containg the N atom labels and a Nx3 array contaning
    the atoms coordinates.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the `pdb` file. Does not neet to contain the `.pdb` extention.

    Returns
    -------
    atonTypes : list
        List of strings containing containg the N atom labels.
    cartPos : numpy array
        Nx3 array contaning the atoms coordinates
    cellMatrix : numpy array
        3x3 array contaning the cell vectors.
    """

    # Remove the extention if exists
    file_name = file_name.split('.')[0]

    if not os.path.exists(os.path.join(path, file_name + '.pdb')):
        raise FileNotFoundError(f'File {file_name} not found!')

    temp_file = open(os.path.join(path, file_name + '.pdb'), 'r').read().splitlines()

    has_cell = any(['CRYST1' in i for i in temp_file])

    if has_cell:
        cellParameters = np.array([i.split()[1:] for i in temp_file if 'CRYST1' in i][0]).astype(float)
        cellMatrix = cellpar_to_cell(cellParameters)
    else:
        cellParameters = np.zeros(6)
        cellMatrix = np.zeros((3, 3))

    if any(['ATOM' in i for i in temp_file]):
        atomTypes = [i.split()[2] for i in temp_file if 'ATOM' in i]
        cartPos = np.array([i.split()[5:8] for i in temp_file if 'ATOM' in i]).astype(float)

    else:
        atomTypes = [i.split()[-1] for i in temp_file if 'HETATM' in i]
        cartPos = np.array([i.split()[3:6] for i in temp_file if 'HETATM' in i]).astype(float)

    return atomTypes, cartPos, cellMatrix


def read_gjf(path, file_name) -> tuple:
    """
    Reads a file in format `.gjf` from the `path` given and returns
    a list containg the N atom labels and a Nx3 array contaning
    the atoms coordinates.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the `gjf` file. Does not neet to contain the `.gjf` extention.

    Returns
    -------
    atonTypes : list
        List of strings containing containg the N atom labels.
    cartPos : numpy array
        Nx3 array contaning the atoms coordinates
    cellMatrix : numpy array
        3x3 array contaning the cell vectors.
    """

    # Remove the extention if exists
    file_name = file_name.split('.')[0]

    if not os.path.exists(os.path.join(path, file_name + '.gjf')):
        raise FileNotFoundError(f'File {file_name} not found!')

    temp_file = open(os.path.join(path, file_name + '.gjf'), 'r').readlines()
    temp_file = [i.split() for i in temp_file if i != '\n']

    atoms = [i for i in temp_file if i[0] in elements_dict()]

    atonTypes = [i[0] for i in atoms]
    cartPos = np.array([[float(i[1]), float(i[2]), float(i[3])] for i in atoms])

    cellMatrix = [i for i in temp_file if 'Tv' in i]

    if cellMatrix:
        cellMatrix = np.array([i[1:] for i in cellMatrix]).astype(float)
    else:
        cellMatrix = np.zeros((3, 3))

    return atonTypes, cartPos, cellMatrix


def read_cif(path, file_name, useASE=False, usePymatgen=False):
    """
    Reads a file in format `.cif` from the `path` given and returns
    a list containg the N atom labels and a Nx3 array contaning
    the atoms coordinates.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the `cif` file. Does not neet to contain the `.cif` extention.
    useASE : bool
        If True, the function will use ASE library to read the file.
    usePymatgen : bool
        If True, the function will use Pymatgen library to read the file.

    Returns
    -------
    atomTypes : list
        List of strings containing containg the N atom labels.
    cartPos : numpy array
        Nx3 array contaning the atoms coordinates
    cellMatrix : numpy array
        3x3 array contaning the cell vectors.
    partialCharges : list
        List of floats containing the N atom partial charges.
    """

    # Remove the extention if exists
    file_name = file_name.split('.')[0]

    if not os.path.exists(os.path.join(path, file_name + '.cif')):
        raise FileNotFoundError(f'File {file_name} not found!')

    # Check if the cif files is in P1 symmetry
    if not useASE and not usePymatgen:
        cif = gemmi.cif.read_file(os.path.join(path, file_name + '.cif')).sole_block()

        # Check if the cif is not with P1 symmetry
        symm_elements = 0
        for loop_name in ['_symmetry_equiv_pos_as_xyz', '_space_group_symop_operation_xyz']:
            if len(cif.find_loop(loop_name)) != 0:
                symm_elements = len(cif.find_loop(loop_name))

        if symm_elements == 0 or symm_elements > 1:
            warnings.warn('The CIF file is not in P1 symmetry. The structure will be read using pyMatGen.')
            usePymatgen = True

    if useASE:
        atoms = read(os.path.join(path, file_name + '.cif'))
        atomTypes = atoms.get_chemical_symbols()  # type: ignore
        cartPos = atoms.get_positions()  # type: ignore
        cellMatrix = np.array(atoms.get_cell())  # type: ignore
        partialCharges = atoms.get_initial_charges()  # type: ignore

    elif usePymatgen:
        structure = CifParser(os.path.join(path, file_name + '.cif')).get_structures(primitive=False)[0]

        atomTypes = [str(i) for i in structure.species]
        cartPos = structure.cart_coords
        cellMatrix = structure.lattice.matrix
        partialCharges = [0 for i in range(len(atomTypes))]

    else:
        cif = gemmi.cif.read_file(os.path.join(path, file_name + '.cif')).sole_block()
        a = float(cif.find_value('_cell_length_a').split('(')[0])
        b = float(cif.find_value('_cell_length_b').split('(')[0])
        c = float(cif.find_value('_cell_length_c').split('(')[0])
        beta = float(cif.find_value('_cell_angle_beta').split('(')[0])
        gamma = float(cif.find_value('_cell_angle_gamma').split('(')[0])
        alpha = float(cif.find_value('_cell_angle_alpha').split('(')[0])

        cellMatrix = cellpar_to_cell([a, b, c, alpha, beta, gamma])

        atomTypes = list(cif.find_values('_atom_site_type_symbol'))

        if len(atomTypes) == 0:
            atomTypes = list(cif.find_values('_atom_site_label'))
            atomTypes = [i.rstrip('0123456789') for i in atomTypes]

        atom_site_fract_x = np.array(cif.find_values('_atom_site_fract_x')).astype(float)
        atom_site_fract_y = np.array(cif.find_values('_atom_site_fract_y')).astype(float)
        atom_site_fract_z = np.array(cif.find_values('_atom_site_fract_z')).astype(float)

        cartPos = np.array([np.dot(cellMatrix, [i, j, k]) for i, j, k in zip(atom_site_fract_x,
                                                                             atom_site_fract_y,
                                                                             atom_site_fract_z)])

        partialCharges = np.array(cif.find_values('_atom_site_charge')).astype(float)

    return atomTypes, cartPos, cellMatrix, partialCharges


def save_csv(path: str, file_name: str, data: list, delimiter: str = ',', head: list = []) -> None:
    """
    Saves a file in format `.csv`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the `csv` file. Does not neet to contain the `.csv` extention.
    data : list
        Data to be saved.
    delimiter: str
        Delimiter of the columns. `,` is the default.
    head : list
        Names of the columns.
    """

    # Remove the extention if exists
    file_name = file_name.split('.')[0]
    file_name = os.path.join(path, file_name + '.csv')

    content = []

    if len(head) > 0:
        content.append(delimiter.join(head))

    for i in data:
        content.append(delimiter.join([str(j) for j in i]))

    with open(file_name, 'w') as f:
        f.write('\n'.join(content))


def save_xsf(path: str,
             file_name: str,
             atomTypes: list,
             atomPos: list,
             cellMatrix: list,
             frac_coords=False,
             **kwargs):
    """
    Save a file in format `.xsf` on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    atomTypes : list
        List of strings containing containg the N atom types
    atomPos : list
        Nx3 array contaning the atoms coordinates.
    cellMatrix : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    frac_coords : bool
        If True, the coordinates are in fractional coordinates.
    """

    file_name = file_name.split('.')[0]

    if len(cellMatrix) == 6:
        cellMatrix = cellpar_to_cell(cellMatrix)  # type: ignore

    if frac_coords:
        # Convert to fractional coordinates
        frac_matrix = get_fractional_to_cartesian_matrix(*cell_to_cellpar(cellMatrix))
        atomPos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atomPos]

    xsf_file = []
    xsf_file.append('CRYSTAL')
    xsf_file.append('PRIMVEC')

    for i in range(len(cellMatrix)):
        xsf_file.append(f'  {cellMatrix[i][0]:>15.9f}    {cellMatrix[i][1]:>15.9f}    {cellMatrix[i][2]:>15.9f}')

    xsf_file.append('   PRIMCOORD')
    xsf_file.append(f'           {len(atomPos)}           1')

    for i in range(len(atomPos)):
        xsf_file.append('{:3s}        {:>15.9f}    {:>15.9f}    {:>15.9f}'.format(atomTypes[i],
                                                                                  atomPos[i][0],
                                                                                  atomPos[i][1],
                                                                                  atomPos[i][2]))

    with open(os.path.join(path, file_name + '.xsf'), 'w') as f:
        f.write('\n'.join(xsf_file))


def save_pqr(path: str,
             file_name: str,
             atomTypes: list,
             atomPos: list,
             cell: list = [],
             partialCharges: list = [],
             frac_coords=False,
             atomLabels: list = [],
             bonds: list = [],
             bond_orders: list = []) -> None:
    """
    Save a file in format `.pqr` on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    atomTypes : list
        List of strings containing containg the N atom types
    atomPos : list
        Nx3 array contaning the atoms coordinates.
    cell : list, optional
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    partialCharges : list, optional
        List of strings containing containg the N atom partial charges.
    frac_coords : bool, optional
        If True, the coordinates are in fractional coordinates and will be converted to cartesian. Default is False.
    atomLabels : list, optional
        List of strings containing containg the N atom labels.
    bonds : list, optional
        List of lists containing the index of the bonded atoms and the bond length.
    bond_orders : list, optional
        List of integers containing the bond orders.
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 3:
        cell = cell_to_cellpar(cell)  # type: ignore

    if frac_coords:
        # Convert to fractional coordinates
        frac_matrix = get_fractional_to_cartesian_matrix(*cell)
        atomPos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atomPos]

    pqr_file = []
    pqr_file.append(f'TITLE       {file_name}')
    pqr_file.append('REMARK   4')
    pqr_file.append(f'REMARK   Created by pyCOFBuilder on {date.today()}')

    if len(cell) > 0:
        pqr_file.append('CRYST1{:>9.3f}{:>9.3f}{:>9.3f}{:>7.2f}{:>7.2f}{:>7.2f} P1'.format(*cell))

    if len(partialCharges) == 0:
        partialCharges = [0 for i in range(len(atomTypes))]

    if len(atomLabels) == 0:
        atomLabels = atomTypes

    atom_line = 'ATOM   {:>4} {:>2}   MOL A   0    {:>8.3f}{:>8.3f}{:>8.3f}{:>8.5f}   {:>15}'
    for i in range(len(atomPos)):
        pqr_file.append(atom_line.format(i + 1,
                                         atomTypes[i],
                                         atomPos[i][0],
                                         atomPos[i][1],
                                         atomPos[i][2],
                                         partialCharges[i],
                                         atomLabels[i]))

    if bonds and not bond_orders:
        bond_orders = [1 for i in range(len(bonds))]

    if bonds:
        for i in range(len(bonds)):
            for j in range(bond_orders[i]):
                pqr_file.append(f'CONECT {bonds[i][0] + 1:4} {bonds[i][1] + 1:4}')

    with open(os.path.join(path, file_name + '.pqr'), 'w') as f:
        f.write('\n'.join(pqr_file))


def save_pdb(path: str,
             file_name: str,
             atomTypes: list,
             atomPos: list,
             cell: list = [],
             frac_coords=False,
             atomLabels: list = [],
             bonds: list = [],
             bond_orders: list = [],
             **kwargs) -> None:
    """
    Save a file in format `.pdb` on the `path`.

    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    atomTypes : list
        List of strings containing containg the N atom types
    atomPos : list
        Nx3 array contaning the atoms coordinates.
    cell : list, optional
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    frac_coords : bool, optional
        If True, the coordinates are in fractional coordinates and will be converted to cartesian. Default is False.
    atomLabels : list, optional
        List of strings containing containg the N atom labels.
    bonds : list, optional
        List of lists containing the index of the bonded atoms and the bond length.
    bond_orders : list, optional
        List of integers containing the bond orders.
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 3:
        cell = cell_to_cellpar(cell)  # type: ignore

    if frac_coords:
        # Convert to fractional coordinates
        frac_matrix = get_fractional_to_cartesian_matrix(*cell)
        atomPos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atomPos]

    pdb_file = []
    pdb_file.append(f'TITLE       {file_name}')
    pdb_file.append(f'REMARK   Created by pyCOFBuilder on {date.today()}')

    if len(cell) > 0:
        pdb_file.append('CRYST1{:>9.3f}{:>9.3f}{:>9.3f}{:>7.2f}{:>7.2f}{:>7.2f} P1'.format(*cell))

    atom_line = 'ATOM   {:>4} {:>2}   MOL     {:>13.3f}{:>8.3f}{:>8.3f}  1.00  0.00  {:>11}'

    for i in range(len(atomPos)):
        pdb_file.append(atom_line.format(i+1,
                                         atomTypes[i],
                                         atomPos[i][0],
                                         atomPos[i][1],
                                         atomPos[i][2],
                                         atomLabels[i]))

    if bonds and not bond_orders:
        bond_orders = [1 for i in range(len(bonds))]

    if bonds:
        for i in range(len(bonds)):
            for j in range(bond_orders[i]):
                pdb_file.append(f'CONECT {bonds[i][0] + 1:4} {bonds[i][1] + 1:4}')

    with open(os.path.join(path, file_name + '.pdb'), 'w') as f:
        f.write('\n'.join(pdb_file))


def save_gjf(path: str,
             file_name: str,
             atomTypes: list,
             atomPos: list,
             cell: list = [],
             frac_coords=False,
             **kwargs) -> None:
    """
    Save a file in format `.gjf` on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    atomTypes : list
        List of strings containing containg the N atom types
    atomPos : list
        Nx3 array contaning the atoms coordinates.
    cell : list, optional
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    frac_coords : bool, optional
        If True, the coordinates are in fractional coordinates and will be converted to cartesian. Default is False.
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 6:
        cell = cellpar_to_cell(cell)  # type: ignore

    if frac_coords:
        # Convert to fractional coordinates
        frac_matrix = get_fractional_to_cartesian_matrix(*cell_to_cellpar(cell))
        atomPos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atomPos]

    gjf_file = []
    gjf_file.append('opt pm6')
    gjf_file.append('')

    gjf_file.append(file_name)
    gjf_file.append('')
    gjf_file.append('0 1')

    for i in range(len(atomPos)):
        gjf_file.append('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atomTypes[i],
                                                                     atomPos[i][0],
                                                                     atomPos[i][1],
                                                                     atomPos[i][2]))
    if cell:
        for i in cell:
            gjf_file.append('Tv   {:>15.7f}{:>15.7f}{:>15.7f}\n'.format(*i))

    gjf_file.append('')
    gjf_file.append('')

    with open(os.path.join(path, file_name + '.gjf'), 'w') as f:
        f.write('\n'.join(gjf_file))


def save_xyz(path: str,
             file_name: str,
             atomTypes: list,
             atomPos: list,
             cell: list = [],
             partialCharges: list = [],
             frac_coords=False,
             **kwargs) -> None:
    """
    Save a file in format `.pqr` on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    atomTypes : list
        List of strings containing containg the N atom types
    atomPos : list
        Nx3 array contaning the atoms coordinates.
    cell : list, optional
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    partialCharges : list, optional
        List of strings containing containg the N atom partial charges.
    frac_coords : bool, optional
        If True, the coordinates are in fractional coordinates and will be converted to cartesian. Default is False.
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 3:
        cell = cell_to_cellpar(cell)  # type: ignore

    if frac_coords:
        # Convert to fractional coordinates
        frac_matrix = get_fractional_to_cartesian_matrix(*cell)
        atomPos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atomPos]

    if len(partialCharges) == 0:
        partialCharges = [0.0 for i in range(len(atomTypes))]

    xyz_file = []
    xyz_file.append(f'{len(atomTypes)}')
    header = ''

    if len(cell) == 6:
        header += 'Lattice="{:>15.7f} {:>15.7f} {:>15.7f} {:>15.7f} {:>15.7f} {:>15.7f}" pbc="T T T"'.format(*cell)

    xyz_file.append(header + ' species:S:1:pos:R:3:charge:R:1')

    for i in range(len(atomTypes)):
        xyz_file.append('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}{:>15.7f}'.format(atomTypes[i],
                                                                            atomPos[i][0],
                                                                            atomPos[i][1],
                                                                            atomPos[i][2],
                                                                            partialCharges[i]))

    with open(os.path.join(path, file_name + '.xyz'), 'w') as f:
        f.write('\n'.join(xyz_file))


def save_turbomole(path: str,
                   file_name: str,
                   atomTypes: list,
                   atomPos: list,
                   cell: list = [],
                   frac_coords=False,
                   **kwargs) -> None:
    """
    Save the structure in Turbomole .coord format on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    atomTypes : list
        List of strings containing containg the N atom types
    atomPos : list
        Nx3 array contaning the atoms coordinates.
    cell : list, optional
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    frac_coords : bool, optional
        If True, the coordinates are in fractional coordinates and will be converted to cartesian. Default is False.
    """

    file_name = file_name.split('.')[0]

    if np.array(cell).shape == (3, 3):
        cell = cell_to_cellpar(cell)  # type: ignore

    if frac_coords:
        # Convert to fractional coordinates
        frac_matrix = get_fractional_to_cartesian_matrix(*cell)
        atomPos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atomPos]

    temp_file = ['$coord angs']

    for i in range(len(atomTypes)):
        temp_file.append('{:>15.7f}{:>15.7f}{:>15.7f}   {:<5s}'.format(atomPos[i][0],
                                                                       atomPos[i][1],
                                                                       atomPos[i][2],
                                                                       atomTypes[i]))

    temp_file.append('$periodic 3')
    temp_file.append('$cell')
    temp_file.append('{}  {}  {}  {}  {}  {}'.format(*cell))
    temp_file.append('$opt')
    temp_file.append('   engine=inertial')
    temp_file.append('$end')

    with open(os.path.join(path, file_name + '.coord'), 'w') as f:
        f.write('\n'.join(temp_file))


def save_vasp(path: str,
              file_name: str,
              atomTypes: list,
              atomPos: list,
              cell: list = [],
              frac_coords=False,
              **kwargs) -> None:
    """
    Save the structure in VASP .vasp format on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    atomTypes : list
        List of strings containing containg the N atom types
    atomPos : list
        Nx3 array contaning the atoms coordinates.
    cell : list, optional
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    frac_coords : bool, optional
        If True, the coordinates are in fractional coordinates and will be converted to cartesian. Default is False.
    """

    file_name = file_name.split('.')[0]

    if np.array(cell).shape == 6:
        cell = cellpar_to_cell(cell)  # type: ignore

    unique_atoms = []
    for i in atomTypes:
        if i not in unique_atoms:
            unique_atoms.append(i)

    composition_dict = {i: atomTypes.count(i) for i in unique_atoms}

    temp_file = []

    temp_file.append(file_name)
    temp_file.append('1.0')

    for i in range(3):
        temp_file.append('{:>15.7f}{:>15.7f}{:>15.7f}'.format(*cell[i]))

    temp_file.append(' '.join(composition_dict.keys()))
    temp_file.append(' '.join([str(i) for i in composition_dict.values()]))

    if frac_coords:
        temp_file.append('Direct')
    else:
        temp_file.append('Cartesian')

    for i in range(len(atomTypes)):
        temp_file.append('{:>15.7f}{:>15.7f}{:>15.7f}   {:<5s}'.format(atomPos[i][0],
                                                                       atomPos[i][1],
                                                                       atomPos[i][2],
                                                                       atomPos[i]))

    with open(os.path.join(path, file_name + '.vasp'), 'w') as f:
        f.write('\n'.join(temp_file))


def save_qe(path: str,
            file_name: str,
            cell: list,
            atomTypes: list,
            atomLabels: list,
            atomPos: list,
            frac_coords=False,
            calc_type: str = 'scf',
            kspacing: float = 0.3,
            **kwargs) -> None:
    """
    Save the structure in Quantum Espresso .pwscf format.

    The `input_dict` can be used to specify the input parameters for the
    QuantumESPRESSO calculation.

    This dictionary must contain the keys: `control`, `system`, `electrons`, and `ions`.
    Each of these keys must contain a dictionary with the corresponding input parameters.
    This dictionary can contain the kpoints item, with the kpoints grid as a list of 3 integers.
    Additionally, it can contain the kspacing item, with the kpoints spacing as a float. In this
    case the kpoints grid will be calculated automatically. By default, the kspacing is set to 0.3.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    atomTypes : list
        List of strings containing containg the N atom types
    atomLabels : list
        List of strings containing containg the N atom labels
    atomPos : list
        Nx3 array contaning the atoms coordinates.
    frac_coords : bool
        If True, the coordinates are in fractional coordinates.
    calc_type : str
        Type of calculation. Can be 'scf', 'vc-relax', 'relax', 'md', 'vc-md', 'vc-tddft', 'tddft'.
    kspacing : float
        Kpoints spacing in 1/Angstrom.
    """

    if len(cell) == 6:
        cell_matrix = cellpar_to_cell(cell)
    else:
        cell_matrix = cell

    ibrav_dict = cell_to_ibrav(cell_matrix)

    input_dict = {}

    input_dict['control'] = {
        'prefix': f"'{file_name}'",
        'calculation': f"{calc_type}",
        'restart_mode': "'from_scratch'",
        'wf_collect': '.true.',
        'pseudo_dir': "'$PSEUDO_DIR'",
        'outdir': "'$SCRATCH_DIR'",
        'verbosity': "'high'",
        'tstress': '.true.',
        'tprnfor': '.true.',
        'etot_conv_thr': '1.0d-5',
        'forc_conv_thr': '1.0d-6',
        'nstep': 1000}

    input_dict['system'] = {
        'nat': len(atomTypes),
        'ntyp': len(set(atomTypes)),
        'ecutwfc': 40,
        'ecutrho': 360,
        'vdw_corr': "'grimme-d3'",
        'occupations': "'smearing'",
        **ibrav_dict}

    input_dict['electrons'] = {
        'conv_thr': 1.0e-9,
        'electron_maxstep': 100,
        'mixing_beta': 0.3}

    if calc_type == 'vc-relax':

        input_dict['ions'] = {
            'ion_dynamics': "'bfgs'"}

        input_dict['cell'] = {
            'cell_dynamics': "'bfgs'",
            'cell_dofree': "'all'"}

    # If the kpoints grid is not specified, calculate it automatically
    if 'k_points' not in input_dict.keys():
        if 'kspacing' not in input_dict.keys():
            input_dict['kspacing'] = kspacing
        input_dict['kpoints'] = get_kgrid(cell_matrix, input_dict['kspacing'])

    with open(os.path.join(path, file_name + '.pwscf'), 'w') as f:
        f.write('&CONTROL\n')
        for key in input_dict['control']:
            f.write(f"  {key} = {input_dict['control'][key]}\n")
        f.write('/\n\n')

        f.write('&SYSTEM\n')
        for key in input_dict['system']:
            f.write(f"  {key} = {input_dict['system'][key]}\n")
        f.write('/\n\n')

        f.write('&ELECTRONS\n')
        for key in input_dict['electrons']:
            f.write(f"  {key} = {input_dict['electrons'][key]}\n")
        f.write('/\n\n')

        if calc_type == 'vc-relax':
            f.write('&IONS\n')
            for key in input_dict['ions']:
                f.write(f"  {key} = {input_dict['ions'][key]}\n")
            f.write('/\n\n')

            f.write('&CELL\n')
            for key in input_dict['cell']:
                f.write(f"  {key} = {input_dict['cell'][key]}\n")
            f.write('/\n\n')

        f.write('ATOMIC_SPECIES\n')
        for atom in set(atomTypes):
            f.write(f" {atom}   {elements_dict()[atom]:>9.5f}  {atom}.PSEUDO.UPF\n")
        f.write('\n')

        # f.write('CELL_PARAMETERS (angstrom) \n')
        # for v in cell_matrix:
        #     f.write(f'{v[0]:>15.9f}      {v[1]:>15.9f}      {v[2]:>15.9f}\n')
        # f.write('\n')

        if frac_coords:
            coords_type = 'crystal'
        else:
            coords_type = 'angstrom'

        f.write(f'ATOMIC_POSITIONS ({coords_type})\n')

        for i, atom in enumerate(atomPos):
            f.write('{:<5s}{:>15.9f}{:>15.9f}{:>15.9f}   ! {:5}\n'.format(atomTypes[i],
                                                                          atom[0],
                                                                          atom[1],
                                                                          atom[2],
                                                                          atomLabels[i]))

        f.write('\n')
        f.write('K_POINTS automatic\n')
        f.write(' {} {} {} 1 1 1\n'.format(*input_dict['kpoints']))


def save_cif(path: str,
             file_name: str,
             atomTypes: list,
             atomPos: list,
             cell: list = [],
             atomLabels: list = [],
             partialCharges: list[float] = [],
             bonds: list = [],
             bondTypes: list = [],
             frac_coords=False,
             renormalize_charges=False,
             **kwargs) -> None:
    """
    Save a file in format `.cif` on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    atomTypes : list
        List of strings containing containg the N atom types
    atomPos : list
        Nx3 array contaning the atoms coordinates.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    atomLabels : list
        List of strings containing containg the N atom labels
    partialCharges : list
        List of strings containing containg the N atom partial charges.
    bonds : list
        List of lists containing the index of the bonded atoms and the bond length.
    bondTypes : list
        List of strings containing the bond types. Can be singe (S), double (D), triple (T), or aromatic (A).
    frac_coords : bool
        If True, the coordinates are in fractional coordinates.
    renormalize_charges : bool
        If True, the charges will be renormalized to the total charge of the system.
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 3:
        a, b, c, alpha, beta, gamma = cell_to_cellpar(cell)
    if len(cell) == 6:
        a, b, c, alpha, beta, gamma = cell

    if renormalize_charges and partialCharges:
        partialCharges = np.array(partialCharges) - np.mean(partialCharges)  # type: ignore

    if not partialCharges:
        partialCharges = [0 for i in range(len(atomTypes))]

    if not atomLabels:
        atomLabels = [f'{atomTypes[i]}{i+1}' for i in range(len(atomTypes))]

    cif_text = f"""\
data_{file_name}

_audit_creation_date     {date.today().strftime("%Y-%d-%m")}
_audit_creation_method   pyCOFBuilder

_chemical_name_common                  '{file_name}'
_cell_length_a                          {a:>10.6f}
_cell_length_b                          {b:>10.6f}
_cell_length_c                          {c:>10.6f}
_cell_angle_alpha                       {alpha:>6.2f}
_cell_angle_beta                        {beta:>6.2f}
_cell_angle_gamma                       {gamma:>6.2f}
_space_group_name_H-M_alt               'P 1'
_space_group_IT_number                  1

loop_
_symmetry_equiv_pos_as_xyz
   'x, y, z'

loop_
   _atom_site_label
   _atom_site_type_symbol
   _atom_site_fract_x
   _atom_site_fract_y
   _atom_site_fract_z
   _atom_site_charge
"""

    if not frac_coords:
        # Convert to fractional coordinates
        frac_matrix = get_cartesian_to_fractional_matrix(a, b, c, alpha, beta, gamma)
        atomPos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atomPos]

    for i, pos in enumerate(atomPos):
        u, v, w = pos[0], pos[1], pos[2]

        cif_text += '{:<15}    {} {:>15.9f} {:>15.9f} {:>15.9f} {:>10.5f}\n'.format(
                atomLabels[i],
                atomTypes[i],
                u,
                v,
                w,
                partialCharges[i])

    if bonds:
        cif_text += '\nloop_\n'
        cif_text += '_geom_bond_atom_site_label_1\n'
        cif_text += '_geom_bond_atom_site_label_2\n'
        cif_text += '_geom_bond_distance\n'
        cif_text += '_geom_bond_site_symmetry_2\n'
        cif_text += '_ccdc_geom_bond_type\n'

        if not bondTypes:
            bondTypes = ['S' for i in range(len(bonds))]

        for i, bond in enumerate(bonds):
            cif_text += f'{atomLabels[bond[0]]:10} {atomLabels[bond[1]]:10} {bond[2]:.5f}  .  {bondTypes[i]}\n'

    with open(os.path.join(path, file_name + '.cif'), 'w') as f:
        f.write(cif_text)


def save_chemjson(path: str,
                  file_name: str,
                  cell: list,
                  atomTypes: list,
                  atomLabels: list,
                  atomPos: list,
                  partialCharges: list = [],
                  bonds: list = [],
                  frac_coords=False):
    """
    Save a file in format `.json` on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    atomTypes : list
        List of strings containing containg the N atom types
    atomLabels : list
        List of strings containing containg the N atom labels
    atomPos : list
        Nx3 array contaning the atoms coordinates.
    partialCharges : list
        List of floats containing containg the N atom partial charges.
    bonds : list
        List of lists containing the index of the bonded atoms and the bond length.
    frac_coords : bool
        If True, the coordinates are in fractional coordinates.
    """

    file_name = file_name.split('.')[0]

    data = ChemJSON()

    if len(cell) == 6:
        data.set_cell_parameters(cell)
    elif len(cell) == 3:
        data.set_cell_matrix(cell)
    else:
        raise ValueError('Cell parameters not provided.')

    data.set_atomic_types(atomTypes)

    if frac_coords:
        data.set_fractional_positions(atomPos)
    else:
        data.set_cartesian_positions(atomPos)

    if partialCharges:
        pass  # TODO: Add partial charges
    if atomLabels:
        pass  # TODO: Add atom labels
    if bonds:
        pass  # TODO: Add bonds

    data.write_cjson(path, file_name + '.cjson')


def generate_mol_dict(path, file_name, name, code, smiles):

    xsmiles, xsmiles_label, composition = smiles_to_xsmiles(smiles)

    if file_name.endswith('gjf'):
        atom_types, atom_pos = read_gjf(path, file_name)
    elif file_name.endswith('xyz'):
        atom_types, atom_pos = read_xyz(path, file_name)

    mol_dict = {
        "name": name,
        "smiles": smiles,
        "code": code,
        "xsmiles": xsmiles,
        "xsmiles_label": xsmiles_label,
        "formula": composition,
        "atoms": {
            "elements": {"elementType": atom_types},
            "coords": {"3d": atom_pos.tolist()}
        }
    }

    print(mol_dict)

    # Write the dictionary to a json file
    with open(os.path.join(path, file_name.split('.')[0] + '.json'), 'w') as f:
        json.dump(mol_dict, f, indent=4)
