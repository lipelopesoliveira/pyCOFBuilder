# -*- coding: utf-8 -*-
# Created by Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

"""
This module contains tools for input and output file manipulation used by pyCOFBuilder.
"""

import os
from datetime import date
import numpy as np

from pymatgen.io.cif import CifParser

import simplejson
from pycofbuilder.tools import (elements_dict,
                                cell_to_cellpar,
                                cellpar_to_cell,
                                get_fractional_to_cartesian_matrix,
                                get_cartesian_to_fractional_matrix,
                                get_kgrid,
                                formula_from_atom_list,
                                smiles_to_xsmiles,
                                cell_to_ibrav)


def save_csv(path, file_name, data, delimiter=',', head=False):
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
    head : str
        Names of the columns.
    """

    # Remove the extention if exists
    file_name = file_name.split('.')[0]
    file_name = os.path.join(path, file_name + '.csv')

    file_temp = open(file_name, 'w')
    if head is not False:
        file_temp.write(head)

    for i in range(len(data)):
        file_temp.write(delimiter.join([str(j) for j in data[i]]) + '\n')

    file_temp.close()


def read_xyz(path, file_name):
    """
    Reads a file in format `.xyz` from the `path` given and returns
    a list containg the N atom labels and a Nx3 array contaning
    the atoms coordinates.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the `xyz` file. Does not neet to contain the `.xyz` extention.

    Returns
    -------
    atom_labels : list
        List of strings containing containg the N atom labels.
    atom_pos : numpy array
        Nx3 array contaning the atoms coordinates
    """

    # Remove the extention if exists
    file_name = file_name.split('.')[0]

    if os.path.exists(os.path.join(path, file_name + '.xyz')):
        temp_file = open(os.path.join(path, file_name + '.xyz'), 'r').readlines()

        atoms = [i.split() for i in temp_file[2:]]

        atom_labels = [i[0] for i in atoms if len(i) > 1]
        atom_pos = np.array([[float(i[1]), float(i[2]), float(i[3])] for i in atoms if len(i) > 1])

        return atom_labels, atom_pos
    else:
        print(f'File {file_name} not found!')
        return None


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
    atom_labels : list
        List of strings containing containg the N atom labels.
    atom_pos : numpy array
        Nx3 array contaning the atoms coordinates
    """

    # Remove the extention if exists
    file_name = file_name.split('.')[0]

    if not os.path.exists(os.path.join(path, file_name + '.pdb')):
        raise FileNotFoundError(f'File {file_name} not found!')

    temp_file = open(os.path.join(path, file_name + '.pdb'), 'r').read().splitlines()

    cellParameters = np.array([i.split()[1:] for i in temp_file if 'CRYST1' in i][0]).astype(float)

    AtomTypes = [i.split()[2] for i in temp_file if 'ATOM' in i]
    CartPos = np.array([i.split()[4:7] for i in temp_file if 'ATOM' in i]).astype(float)

    return cellParameters, AtomTypes, CartPos


def read_gjf(path, file_name):
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
    atom_labels : list
        List of strings containing containg the N atom labels.
    atom_pos : numpy array
        Nx3 array contaning the atoms coordinates
    """
    # Remove the extention if exists
    file_name = file_name.split('.')[0]

    if os.path.exists(os.path.join(path, file_name + '.gjf')):

        temp_file = open(os.path.join(path, file_name + '.gjf'), 'r').readlines()
        temp_file = [i.split() for i in temp_file if i != '\n']

        atoms = [i for i in temp_file if i[0] in elements_dict()]

        atom_labels = [i[0] for i in atoms]
        atom_pos = np.array([[float(i[1]), float(i[2]), float(i[3])] for i in atoms])

        return atom_labels, atom_pos
    else:
        print(f'File {file_name} not found!')
        return None


def read_cif(path, file_name):
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

    Returns
    -------
    cell : numpy array
        3x3 array contaning the cell vectors.
    atom_labels : list
        List of strings containing containg the N atom labels.
    atom_pos : numpy array
        Nx3 array contaning the atoms coordinates
    charges : list
        List of strings containing containg the N atom partial charges.
    """

    # Remove the extention if exists
    file_name = file_name.split('.')[0]

    if os.path.exists(os.path.join(path, file_name + '.cif')):

        temp_file = open(os.path.join(path, file_name + '.cif'), 'r').readlines()
        cell = []
        atom_label = []
        atom_pos = []
        charges = []
        has_charges = False
        for i in temp_file:
            if 'cell_length_a' in i:
                cell += [float(i.split()[-1])]
            if 'cell_length_b' in i:
                cell += [float(i.split()[-1])]
            if 'cell_length_c' in i:
                cell += [float(i.split()[-1])]
            if 'cell_angle_alpha' in i:
                cell += [float(i.split()[-1])]
            if '_cell_angle_beta' in i:
                cell += [float(i.split()[-1])]
            if '_cell_angle_gamma' in i:
                cell += [float(i.split()[-1])]
            if '_atom_site_charge' in i:
                has_charges = True

        for i in temp_file:
            line = i.split()
            if len(line) > 1 and line[0] in elements_dict().keys():
                atom_label += [line[0]]
                atom_pos += [[float(j) for j in line[2:5]]]
                if has_charges:
                    charges += [float(line[-1])]
        cell = cellpar_to_cell(cell)

        return cell, atom_label, atom_pos, charges
    else:
        print(f'File {file_name} not found!')
        return None


def save_xsf(path: str,
             file_name: str,
             cell: list,
             atom_types: list,
             atom_labels: list,
             atom_pos: list,
             atom_charges: list = None,
             bonds: list = None,
             bond_orders: list = None,
             frac_coords=False):
    """
    Save a file in format `.xsf` on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    atom_types : list
        List of strings containing containg the N atom types
    atom_label : list
        List of strings containing containg the N atom labels
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    atom_charges : list
        List of strings containing containg the N atom partial charges.
    bonds : list
        List of lists containing the index of the bonded atoms and the bond length.
    frac_coords : bool
        If True, the coordinates are in fractional coordinates.
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 6:
        cell = cellpar_to_cell(cell)

    if frac_coords:
        # Convert to fractional coordinates
        frac_matrix = get_fractional_to_cartesian_matrix(*cell_to_cellpar(cell))
        atom_pos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atom_pos]

    xsf_file = open(os.path.join(path, file_name + '.xsf'), 'w')
    xsf_file.write(' CRYSTAL\n')
    xsf_file.write('  PRIMVEC\n')

    for i in range(len(cell)):
        xsf_file.write(f'  {cell[i][0]:>15.9f}    {cell[i][1]:>15.9f}    {cell[i][2]:>15.9f}\n')

    xsf_file.write('   PRIMCOORD\n')
    xsf_file.write(f'           {len(atom_pos)}           1\n')

    for i in range(len(atom_pos)):
        xsf_file.write('{:3s}        {:>15.9f}    {:>15.9f}    {:>15.9f}\n'.format(atom_types[i],
                                                                                   atom_pos[i][0],
                                                                                   atom_pos[i][1],
                                                                                   atom_pos[i][2]))

    xsf_file.close()


def save_pqr(path: str,
             file_name: str,
             cell: list,
             atom_types: list,
             atom_labels: list,
             atom_pos: list,
             atom_charges: list = None,
             bonds: list = None,
             bond_orders: list = None,
             frac_coords=False):
    """
    Save a file in format `.pqr` on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    atom_types : list
        List of strings containing containg the N atom types
    atom_label : list
        List of strings containing containg the N atom labels
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    atom_charges : list
        List of strings containing containg the N atom partial charges.
    bonds : list
        List of lists containing the index of the bonded atoms and the bond length.
    frac_coords : bool
        If True, the coordinates are in fractional coordinates.
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 3:
        cell = cell_to_cellpar(cell)

    if frac_coords:
        # Convert to fractional coordinates
        frac_matrix = get_fractional_to_cartesian_matrix(*cell)
        atom_pos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atom_pos]

    pqr_file = open(os.path.join(path, file_name + '.pqr'), 'w')
    pqr_file.write(f'TITLE       {file_name}  \n')
    pqr_file.write('REMARK   4\n')
    pqr_file.write('CRYST1{:>9.3f}{:>9.3f}{:>9.3f}{:>7.2f}{:>7.2f}{:>7.2f} P1\n'.format(cell[0],
                                                                                        cell[1],
                                                                                        cell[2],
                                                                                        cell[3],
                                                                                        cell[4],
                                                                                        cell[5]))

    if atom_charges is None:
        atom_line = 'ATOM   {:>4} {:>2}   MOL A   0    {:>8.3f}{:>8.3f}{:>8.3f}   {:>15}\n'
        for i in range(len(atom_pos)):
            pqr_file.write(atom_line.format(i + 1,
                                            atom_types[i],
                                            atom_pos[i][0],
                                            atom_pos[i][1],
                                            atom_pos[i][2],
                                            atom_types[i]))
    else:
        atom_line = 'ATOM   {:>4} {:>2}   MOL A   0    {:>8.3f}{:>8.3f}{:>8.3f}{:>8.5f}   {:>15}\n'
        for i in range(len(atom_pos)):
            pqr_file.write(atom_line.format(i + 1,
                                            atom_types[i],
                                            atom_pos[i][0],
                                            atom_pos[i][1],
                                            atom_pos[i][2],
                                            atom_charges[i],
                                            atom_types[i]))

    if bonds and not bond_orders:
        bond_orders = [1 for i in range(len(bonds))]

    if bonds:
        for i in range(len(bonds)):
            for j in range(bond_orders[i]):
                pqr_file.write(f'CONECT {bonds[i][0] + 1:4} {bonds[i][1] + 1:4}\n')

    pqr_file.close()


def save_pdb(path: str,
             file_name: str,
             cell: list,
             atom_types: list,
             atom_labels: list,
             atom_pos: list,
             atom_charges: list = None,
             bonds: list = None,
             bond_orders: list = None,
             frac_coords=False):
    """
    Save a file in format `.pdb` on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    atom_types : list
        List of strings containing containg the N atom types
    atom_label : list
        List of strings containing containg the N atom labels
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    atom_charges : list
        List of strings containing containg the N atom partial charges.
    bonds : list
        List of lists containing the index of the bonded atoms and the bond length.
    frac_coords : bool
        If True, the coordinates are in fractional coordinates.
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 3:
        cell = cell_to_cellpar(cell)

    if frac_coords:
        # Convert to fractional coordinates
        frac_matrix = get_fractional_to_cartesian_matrix(*cell)
        atom_pos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atom_pos]

    pdb_file = open(os.path.join(path, file_name + '.pdb'), 'w')
    pdb_file.write(f'TITLE       {file_name}  \n')
    pdb_file.write('REMARK   pyCOFBuilder\n')
    pdb_file.write('CRYST1{:>9.3f}{:>9.3f}{:>9.3f}{:>7.2f}{:>7.2f}{:>7.2f} P1\n'.format(cell[0],
                                                                                        cell[1],
                                                                                        cell[2],
                                                                                        cell[3],
                                                                                        cell[4],
                                                                                        cell[5]))

    atom_line = 'ATOM   {:>4} {:>2}   MOL     {:>13.3f}{:>8.3f}{:>8.3f}  1.00  0.00  {:>11}\n'

    for i in range(len(atom_pos)):
        pdb_file.write(atom_line.format(i+1,
                                        atom_types[i],
                                        atom_pos[i][0],
                                        atom_pos[i][1],
                                        atom_pos[i][2],
                                        atom_types[i]))

    if bonds and not bond_orders:
        bond_orders = [1 for i in range(len(bonds))]

    if bonds:
        for i in range(len(bonds)):
            for j in range(bond_orders[i]):
                pdb_file.write(f'CONECT {bonds[i][0] + 1:4} {bonds[i][1] + 1:4}\n')

    pdb_file.close()


def save_gjf(path: str,
             file_name: str,
             cell: list,
             atom_types: list,
             atom_labels: list,
             atom_pos: list,
             atom_charges: list = None,
             bonds: list = None,
             bond_orders: list = None,
             frac_coords=False,
             header: str = 'opt pm6'):
    """
    Save a file in format `.pqr` on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    atom_types : list
        List of strings containing containg the N atom types
    atom_label : list
        List of strings containing containg the N atom labels
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    atom_charges : list
        List of strings containing containg the N atom partial charges.
    bonds : list
        List of lists containing the index of the bonded atoms and the bond length.
    frac_coords : bool
        If True, the coordinates are in fractional coordinates.
    header : str
        Parameters for Gaussian calculations.
    """
    if len(cell) == 6:
        cell = cellpar_to_cell(cell)

    if frac_coords:
        # Convert to fractional coordinates
        frac_matrix = get_fractional_to_cartesian_matrix(*cell_to_cellpar(cell))
        atom_pos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atom_pos]

    file_name = file_name.split('.')[0]

    temp_file = open(os.path.join(path, file_name + '.gjf'), 'w')
    temp_file.write(f'%chk={file_name}.chk \n')
    temp_file.write(f'# {header}\n')
    temp_file.write('\n')
    temp_file.write(f'{file_name}\n')
    temp_file.write('\n')
    temp_file.write('0 1 \n')

    for i in range(len(atom_types)):
        temp_file.write('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atom_types[i],
                                                                     atom_pos[i][0],
                                                                     atom_pos[i][1],
                                                                     atom_pos[i][2]))
    if cell is not None:
        for i in range(len(cell)):
            temp_file.write('Tv   {:>15.7f}{:>15.7f}{:>15.7f}\n'.format(*cell[i]))

    temp_file.write('\n\n')
    temp_file.close()


def save_xyz(path: str,
             file_name: str,
             atom_types: list,
             atom_pos: list,
             atom_labels: list = None,
             cell: list = None,
             atom_charges: list = None,
             bonds: list = None,
             bond_orders: list = None,
             frac_coords=False):
    """
    Save a file in format `.xyz` on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    atom_types : list
        List of strings containing containg the N atom types
    atom_label : list
        List of strings containing containg the N atom labels
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    atom_charges : list
        List of strings containing containg the N atom partial charges.
    bonds : list
        List of lists containing the index of the bonded atoms and the bond length.
    frac_coords : bool
        If True, the coordinates are in fractional coordinates.
    """

    if cell:
        cell = cell_to_cellpar(cell) if len(cell) == 3 else cell

    if frac_coords:
        # Convert to fractional coordinates
        frac_matrix = get_fractional_to_cartesian_matrix(cell)
        atom_pos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atom_pos]

    file_name = file_name.split('.')[0]

    temp_file = open(os.path.join(path, file_name + '.xyz'), 'w')
    temp_file.write(f'{len(atom_types)}\n')

    if cell is None:
        temp_file.write(f'{file_name}\n')
    else:
        temp_file.write(f'{cell[0]}  {cell[1]}  {cell[2]}  {cell[3]}  {cell[4]}  {cell[5]}\n')

    for i in range(len(atom_types)):
        temp_file.write('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atom_types[i],
                                                                     atom_pos[i][0],
                                                                     atom_pos[i][1],
                                                                     atom_pos[i][2]))

    temp_file.close()


def save_turbomole(path: str,
                   file_name: str,
                   cell: list,
                   atom_types: list,
                   atom_labels: list,
                   atom_pos: list,
                   atom_charges: list = None,
                   bonds: list = None,
                   bond_orders: list = None,
                   frac_coords=False):
    """
    Save the structure in Turbomole .coord format on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    atom_types : list
        List of strings containing containg the N atom types
    atom_label : list
        List of strings containing containg the N atom labels
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    atom_charges : list
        List of strings containing containg the N atom partial charges.
    bonds : list
        List of lists containing the index of the bonded atoms and the bond length.
    frac_coords : bool
        If True, the coordinates are in fractional coordinates.
    """

    if np.array(cell).shape == (3, 3):
        cell = cell_to_cellpar(cell)

    if frac_coords:
        # Convert to fractional coordinates
        frac_matrix = get_fractional_to_cartesian_matrix(*cell)
        atom_pos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atom_pos]

    with open(os.path.join(path, file_name + '.coord'), 'w') as temp_file:
        temp_file.write('$coord angs\n')

        for i in range(len(atom_types)):
            temp_file.write('{:>15.7f}{:>15.7f}{:>15.7f}   {:<5s}\n'.format(atom_pos[i][0],
                                                                            atom_pos[i][1],
                                                                            atom_pos[i][2],
                                                                            atom_types[i]))

        temp_file.write('$periodic 3\n')
        temp_file.write('$cell\n')
        temp_file.write('{}  {}  {}  {}  {}  {}\n'.format(*cell))
        temp_file.write('$opt\n')
        temp_file.write('   engine=inertial\n')
        temp_file.write('$end\n')


def save_vasp(path: str,
              file_name: str,
              cell: list,
              atom_types: list,
              atom_labels: list,
              atom_pos: list,
              atom_charges: list = None,
              bonds: list = None,
              bond_orders: list = None,
              frac_coords=False):
    """
    Save the structure in VASP .vasp format on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    atom_types : list
        List of strings containing containg the N atom types
    atom_label : list
        List of strings containing containg the N atom labels
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    atom_charges : list
        List of strings containing containg the N atom partial charges.
    bonds : list
        List of lists containing the index of the bonded atoms and the bond length.
    frac_coords : bool
        If True, the coordinates are in fractional coordinates.
    """

    if np.array(cell).shape == 6:
        cell = cellpar_to_cell(cell)

    unique_atoms = []
    for i in atom_types:
        if i not in unique_atoms:
            unique_atoms.append(i)

    composition_dict = {i: atom_types.count(i) for i in unique_atoms}

    with open(os.path.join(path, file_name + '.vasp'), 'w') as temp_file:
        temp_file.write(f'{file_name}\n')
        temp_file.write('1.0\n')

        for i in range(3):
            temp_file.write('{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(cell[i][0],
                                                                   cell[i][1],
                                                                   cell[i][2]))

        temp_file.write(' '.join(composition_dict.keys()) + '\n')
        temp_file.write(' '.join([str(i) for i in composition_dict.values()]) + '\n')

        if frac_coords:
            temp_file.write('Direct\n')
        else:
            temp_file.write('Cartesian\n')

        for i in range(len(atom_types)):
            temp_file.write('{:>15.7f}{:>15.7f}{:>15.7f}   {:<5s}\n'.format(atom_pos[i][0],
                                                                            atom_pos[i][1],
                                                                            atom_pos[i][2],
                                                                            atom_types[i]))


def save_qe(path: str,
            file_name: str,
            cell: list,
            atom_types: list,
            atom_labels: list,
            atom_pos: list,
            atom_charges: list = None,
            bonds: list = None,
            bond_orders: list = None,
            frac_coords=False,
            calc_type: str = 'scf',
            kspacing: float = 0.3):
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
    atom_types : list
        List of strings containing containg the N atom types
    atom_label : list
        List of strings containing containg the N atom labels
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    atom_charges : list
        List of strings containing containg the N atom partial charges.
    bonds : list
        List of lists containing the index of the bonded atoms and the bond length.
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
        'nat': len(atom_types),
        'ntyp': len(set(atom_types)),
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
        for atom in set(atom_types):
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

        for i, atom in enumerate(atom_pos):
            f.write('{:<5s}{:>15.9f}{:>15.9f}{:>15.9f}   ! {:5}\n'.format(atom_types[i],
                                                                          atom[0],
                                                                          atom[1],
                                                                          atom[2],
                                                                          atom_labels[i]))

        f.write('\n')
        f.write('K_POINTS automatic\n')
        f.write(' {} {} {} 1 1 1\n'.format(*input_dict['kpoints']))


def convert_cif_2_qe(out_path, file_name):
    """
    Convert a cif file to a Quantum Espresso input file

    Parameters
    ----------
    out_path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.cif` extention.
    """

    cell, atom_labels, atom_pos, _ = read_cif(out_path, file_name, has_charges=False)

    print(cell, atom_labels, atom_pos)

    save_qe(out_path,
            file_name,
            cell,
            atom_labels,
            atom_pos,
            coords_are_cartesian=True,
            supercell=False,
            angs=False,
            ecut=40,
            erho=360,
            k_dist=0.3)


def save_chemjson(path: str,
                  file_name: str,
                  cell: list,
                  atom_types: list,
                  atom_labels: list,
                  atom_pos: list,
                  atom_charges: list = None,
                  bonds: list = None,
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
    atom_types : list
        List of strings containing containg the N atom types
    atom_label : list
        List of strings containing containg the N atom labels
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    atom_charges : list
        List of strings containing containg the N atom partial charges.
    bonds : list
        List of lists containing the index of the bonded atoms and the bond length.
    frac_coords : bool
        If True, the coordinates are in fractional coordinates.
    """

    file_name = file_name.split('.')[0]
    if len(cell) == 6:
        CellParameters = cell
        CellMatrix = None
    if len(cell) == 3:
        CellParameters = None
        CellMatrix = cell

    chemJSON = create_structure_CJSON(StructureName=file_name.split('.')[0],
                                      CellParameters=CellParameters,
                                      CellMatrix=CellMatrix,
                                      AtomTypes=atom_types,
                                      AtomPositions=atom_pos,
                                      AtomLabels=atom_labels,
                                      CartesianPositions=not frac_coords,
                                      BondIndexes=bonds)

    write_json(path, file_name, chemJSON)


def save_cif(path: str,
             file_name: str,
             cell: list,
             atom_types: list,
             atom_labels: list,
             atom_pos: list,
             atom_charges: list = None,
             bonds: list = None,
             frac_coords=False):
    """
    Save a file in format `.cif` on the `path`.

    Parameters
    ----------
    path : str
        Path to the save the file.
    file_name : str
        Name of the file. Does not neet to contain the extention.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    atom_types : list
        List of strings containing containg the N atom types
    atom_label : list
        List of strings containing containg the N atom labels
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    atom_charges : list
        List of strings containing containg the N atom partial charges.
    bonds : list
        List of lists containing the index of the bonded atoms and the bond length.
    frac_coords : bool
        If True, the coordinates are in fractional coordinates.
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 3:
        a, b, c, alpha, beta, gamma = cell_to_cellpar(cell)
    if len(cell) == 6:
        a, b, c, alpha, beta, gamma = cell

    if atom_labels is None:
        atom_labels = [''] * len(atom_types)

    cif_text = f"""\
data_{file_name}

_audit_creation_date     {date.today().strftime("%Y-%d-%m")}
_audit_creation_method   pyCOFBuilder
_audit_author_name       '{os.getlogin()}'

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
"""

    if atom_charges:
        cif_text += '   _atom_site_charge\n'

    if frac_coords is False:
        # Convert to fractional coordinates
        frac_matrix = get_cartesian_to_fractional_matrix(a, b, c, alpha, beta, gamma)
        atom_pos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atom_pos]

    for i in range(len(atom_pos)):
        u, v, w = atom_pos[i][0], atom_pos[i][1], atom_pos[i][2]
        if atom_charges:
            atom_labels[i] = f"{atom_types[i]}{str(i + 1)}_{atom_labels[i]}"
            cif_text += '{:<15}    {} {:>15.9f} {:>15.9f} {:>15.9f} {:>10.5f}\n'.format(
                f"{atom_types[i]}{str(i + 1)}_{atom_labels[i]}",
                atom_types[i],
                u,
                v,
                w,
                atom_charges[i])
        else:
            atom_labels[i] = f"{atom_types[i]}{str(i + 1)}_{atom_labels[i]}"
            cif_text += '{:<15}    {} {:>15.9f} {:>15.9f} {:>15.9f}\n'.format(
                f"{atom_types[i]}{str(i + 1)}_{atom_labels[i]}",
                atom_types[i],
                u,
                v,
                w)

    if bonds:
        cif_text += '\nloop_\n'
        cif_text += '_geom_bond_atom_site_label_1\n'
        cif_text += '_geom_bond_atom_site_label_2\n'
        cif_text += '_geom_bond_distance\n'

        for bond in bonds:
            cif_text += f'{atom_labels[bond[0]]:10} {atom_labels[bond[1]]:10} {bond[2]:.5f}\n'

    # Write cif_text to file
    cif_file = open(os.path.join(path, file_name + '.cif'), 'w')
    cif_file.write(cif_text)
    cif_file.close()


def convert_json_2_cif(origin_path, file_name, destiny_path, charge_type='None'):
    """
    Convert a file in format `.json` to `.cif`.

    Parameters
    ----------
    origin_path : str
        Path to the '.json' file.
    file_name : str
        Name of the file. Does not neet to contain the `.json` extention.
    destiny_path : str
        path where the `.cif` file will be saved.
    """

    framework_JSON = read_json(origin_path, file_name)

    cell = framework_JSON['geometry']['cell_matrix']
    atom_labels = framework_JSON['geometry']['atom_labels']
    atom_pos = framework_JSON['geometry']['atom_pos']

    if charge_type + '_charges' in list(framework_JSON['system'].keys()):
        partial_charges = framework_JSON['geometry'][charge_type + '_charges']
    else:
        partial_charges = False

    save_cif(destiny_path,
             file_name,
             cell,
             atom_labels,
             atom_pos,
             partial_charges,
             frac_coords=False)


def convert_gjf_2_xyz(path, file_name):

    file_name = file_name.split('.')[0]

    atom_labels, atom_pos = read_gjf(path, file_name + '.gjf')

    save_xyz(path, file_name + '.xyz', atom_labels, atom_pos)


def convert_xyz_2_gjf(path, file_name):

    file_name = file_name.split('.')[0]

    atom_labels, atom_pos = read_xyz(path, file_name + '.xyz')

    save_gjf(path=path,
             file_name=file_name + '.gjf',
             atom_types=atom_labels,
             atom_pos=atom_pos,
             cell=[10, 10, 10, 90, 90, 90])


def convert_cif_2_xyz(path, file_name, supercell=[1, 1, 1]):

    file_name = file_name.split('.')[0]

    structure = CifParser(os.path.join(path, file_name + '.cif')).get_structures(primitive=True)[0]

    structure.make_supercell([[supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, supercell[2]]])

    dict_sctructure = structure.as_dict()

    a, b, c = dict_sctructure['lattice']['a']
    b = dict_sctructure['lattice']['b']
    c = dict_sctructure['lattice']['c']

    alpha = round(dict_sctructure['lattice']['alpha'])
    beta = round(dict_sctructure['lattice']['beta'])
    gamma = round(dict_sctructure['lattice']['gamma'])

    atom_labels = [i['label'] for i in dict_sctructure['sites']]

    atom_pos = [i['xyz'] for i in dict_sctructure['sites']]

    temp_file = open(os.path.join(path, file_name + '.xyz'), 'w')
    temp_file.write(f'{len(atom_labels)} \n')

    temp_file.write(f'{a}  {b}  {c}  {alpha}  {beta}  {gamma}\n')

    for i in range(len(atom_labels)):
        temp_file.write('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atom_labels[i],
                                                                     atom_pos[i][0],
                                                                     atom_pos[i][1],
                                                                     atom_pos[i][2]))

    temp_file.close()


def write_json(path, name, COF_json):

    name = name.split('.')[0]

    if os.path.exists(path) is not True:
        os.mkdir(path)

    save_path = os.path.join(path, name + '.cjson')

    with open(save_path, 'w', encoding='utf-8') as f:
        simplejson.dump(COF_json,
                        f,
                        ensure_ascii=False,
                        separators=(',', ':'),
                        indent=2,
                        ignore_nan=True)


def read_json(path, name):

    cof_path = os.path.join(path, name + '.json')

    with open(cof_path, 'r') as r:
        json_object = simplejson.loads(r.read())

    return json_object


def create_COF_json(name) -> dict:
    """
    Create a empty dictionary with the COF information.
    """

    system_info = 'Informations about the system.'
    geometry_info = 'Informations about the geometry.'
    optimization_info = 'Information about the optimization process.'
    adsorption_info = 'Information about the adsorption simulation experiments on RASPA2'
    textural_info = 'Information about the textural properties'
    spectrum_info = 'Information about spectra simulation.'
    experimental_info = 'Experimental data DRX, FTIR, ssNMR, UV-VIS...'

    COF_json = {'system': {'description': system_info,
                           'name': name,
                           'geo_opt': True,
                           'execution_times_seconds': {}},
                'geometry': {'description': geometry_info},
                'optimization': {'description': optimization_info},
                'adsorption': {'description': adsorption_info},
                'textural': {'description': textural_info},
                'spectrum': {'description': spectrum_info},
                'experimental': {'description': experimental_info}
                }

    return COF_json


def create_empty_CJSON() -> dict:
    """
    Create a dictionary with the structure information to be saved using the
    chemical JSON format.
    """

    chemJSON = {
        "chemicalJson": 1,
        "name": "",
        "formula": "",
        "unitCell": {
            "a": 0.0,
            "b": 0.0,
            "c": 0.0,
            "alpha": 0.0,
            "beta":  0.0,
            "gamma": 0.0,
            "cellVectors": []
        },
        "atoms": {
            "elements": {
                "number": [],
                "type": []
                },
            "coords": {
                "3d": [],
                "3dFractional": []
            },
            "formalCharges": [],
            "labels": []
        },
        "bonds": {
                "connections": {
                    "index": []
                },
                "order": []
            },
        "PartialCharges": {},
        "properties": {
            "molecularMass": 0,
            "totalCharge": 0,
            "spinMultiplicity": 1,
            "totalEnergy": 0,
            "bandGap": 0,
        },
        "spectra": {},
        "vibrations": {},
        "metadata": {},
    }

    return chemJSON


def create_structure_CJSON(StructureName: str,
                           CellParameters: list = None,
                           CellMatrix: list = None,
                           AtomTypes: list = None,
                           AtomLabels: list = [],
                           AtomPositions: list = None,
                           CartesianPositions: bool = False,
                           BondIndexes: list = [],
                           BondOrders: list = [],
                           PartialCharges: dict = {},
                           ) -> dict:
    """
    Creates a dictionary with the structure information to be saved using the
    chemical JSON format.

    Parameters
    ----------
    StructureName : str
        Name of the structure.
    CellParameters : list
        List with the cell parameters.
    CellMatrix : list
        List with the cell matrix. Optional
    AtomTypes : list
        List with the atom types.
    AtomLabels : list
        List with the atom labels.
    AtomPositions : list
        List with the atom positions.
    CartesianPositions : bool
        If True, the coordinates are in cartesian coordinates.
    BondIndexes : list
        List with the bonds indexes and bond length.
    BondOrders : list
        List with the bond orders.
    PartialCharges : dict
        Dictionary with the partial charges.
    """

    chemJSON = create_empty_CJSON()

    chemJSON['name'] = StructureName
    chemJSON['formula'] = formula_from_atom_list(AtomTypes)

    if CellParameters is not None:
        CellMatrix = cellpar_to_cell(CellParameters)
    else:
        CellParameters = cell_to_cellpar(CellMatrix)
        CellMatrix = np.array(CellMatrix)

    chemJSON['unitCell']['a'] = CellParameters[0]
    chemJSON['unitCell']['b'] = CellParameters[1]
    chemJSON['unitCell']['c'] = CellParameters[2]
    chemJSON['unitCell']['alpha'] = CellParameters[3]
    chemJSON['unitCell']['beta'] = CellParameters[4]
    chemJSON['unitCell']['gamma'] = CellParameters[5]
    chemJSON['unitCell']['cellVectors'] = CellMatrix.flatten().tolist()

    AtomNumbers = [elements_dict(property="atomic_number")[i] for i in AtomTypes]
    chemJSON['atoms']['elements']['number'] = AtomNumbers
    chemJSON['atoms']['elements']['type'] = AtomTypes

    chemJSON['atoms']['elements']['labels'] = AtomLabels

    if CartesianPositions:
        chemJSON['atoms']['coords']['3d'] = np.array(AtomPositions).flatten().tolist()
        V_frac = get_cartesian_to_fractional_matrix(*CellParameters)
        FracPosition = np.array([np.dot(V_frac, atom) for atom in AtomPositions]).flatten().tolist()
        chemJSON['atoms']['coords']['3dFractional'] = FracPosition

    else:
        chemJSON['atoms']['coords']['3dFractional'] = np.array(AtomPositions).flatten().tolist()
        V_cart = get_fractional_to_cartesian_matrix(*CellParameters)
        CartPosition = np.array([np.dot(V_cart, atom) for atom in AtomPositions]).flatten().tolist()
        chemJSON['atoms']['coords']['3d'] = CartPosition

    if PartialCharges != {}:
        chemJSON['atoms']['PartialCharges'] = PartialCharges

    bond_indexes = [[i[0], i[1]] for i in BondIndexes]
    bond_indexes = [item for row in bond_indexes for item in row]

    bond_orders = BondOrders if BondOrders != [] else [1] * len(bond_indexes)

    chemJSON['bonds']['connections']['index'] = bond_indexes
    chemJSON['bonds']['order'] = bond_orders

    return chemJSON


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

    write_json(path, file_name.split('.')[0], mol_dict)
