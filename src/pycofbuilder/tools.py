# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020

@author: Felipe Lopes de Oliveira
"""

import os
import simplejson
import numpy as np
from scipy.spatial import distance


def elements_dict(property='atomic_mass'):
    '''Returns a dictionary containing the elements symbol and its selected property.

    Parameters
    ----------
    prop : string
        The desired property can be:
            - "full_name"
            - "atomic_number"
            - "atomic_mass"
            - "polarizability"
            - "pauling_electronegativity"
            - "thermo_electronegativity"
            - "mulliken_electronegativity"
            - "sanderson_electronegativity"
            - "allen_electronegativity"
            - "ghosh_electronegativity"
            - "martynov_batsanov_electronegativity"
            - "atomic_radius"
            - "covalent_radius"
            - "vdw_radius"

    Returns
    -------
    prop_dic : dictionary
        A dictionary containing the elements symbol and its respective property.
    '''

    file_name = os.path.join(os.path.dirname(__file__), 'data', 'periodic_table.json')

    with open(file_name, 'r') as f:
        periodic_table = simplejson.load(f)

    prop_list = periodic_table['H'].keys()

    # Check if the property is valid
    if property not in prop_list:
        raise ValueError('Invalid property. Valid properties are: ' + ', '.join(prop_list))

    prop_dic = {}
    for element in periodic_table:
        prop_dic[element] = periodic_table[element][property]

    return prop_dic


def unit_vector(vector):
    """Return a unit vector in the same direction as x."""
    y = np.array(vector, dtype='float')
    return y / np.linalg.norm(y)


def angle(v1, v2, unit='degree'):
    """
    Calculates the angle between two vectors v1 and v2.
    ----------
    v1 : array
        (N,1) matrix with N dimensions
    v2 : array
        (N,1) matrix with N dimensions
    unit : str
        Unit of the output. Could be 'degree', 'radians' or 'cos'.
    Returns
    -------
    angle : float
        Angle in the selected unit.
    """
    unit_vector1 = unit_vector(v1)
    unit_vector2 = unit_vector(v2)

    dot_product = np.dot(unit_vector1, unit_vector2)

    if unit == 'degree':
        angle = np.arccos(dot_product) * 180. / np.pi
    if unit == 'radians':
        angle = np.arccos(dot_product)
    if unit == 'cos':
        angle = dot_product
    return angle


def rotation_matrix_from_vectors(vec1, vec2):
    '''
    Find the rotation matrix that aligns vec1 to vec2
    ----------
    vec1 : array
        (3,3) array
    vec2 : array
        (3,3) array
    Returns
    -------
    rotation_matrix : array
        A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    '''
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    if s != 0:
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
        return rotation_matrix
    else:
        return np.identity(3)


def translate_inside(matrix):
    """
    Translate the matrix to the center of the cell
    """
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] >= 1:
                matrix[i][j] -= 1
    return matrix


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    Parameters
    ----------
    V : array
        (N,D) matrix, where N is points and D is dimension.
    W : array
        (N,D) matrix, where N is points and D is dimension.
    Returns
    -------
    rmsd : float
        Root-mean-square deviation between the two vectors
    """
    diff = np.array(V) - np.array(W)
    N = len(V)
    return np.sqrt((diff * diff).sum() / N)


def cell_to_cellpar(cell, radians=False):
    """Returns the cell parameters [a, b, c, alpha, beta, gamma]
    given a 3x3 cell matrix [v1, v2, v3]

    Angles are in degrees unless radian=True is used.
    """
    lengths = [np.linalg.norm(v) for v in cell]
    angles = []
    for i in range(3):
        j = i - 1
        k = i - 2
        ll = lengths[j] * lengths[k]
        if ll > 1e-16:
            x = np.dot(cell[j], cell[k]) / ll
            angle = 180.0 / np.pi * np.arccos(x)
        else:
            angle = 90.0
        angles.append(angle)
    if radians:
        angles = [angle * np.pi / 180 for angle in angles]
    return np.array(lengths + angles)


def cellpar_to_cell(cellpar, ab_normal=(0, 0, 1), a_direction=None):
    """Return a 3x3 cell matrix from cellpar=[a,b,c,alpha,beta,gamma].

    Angles must be in degrees.

    The returned cell is orientated such that a and b
    are normal to `ab_normal` and a is parallel to the projection of
    `a_direction` in the a-b plane.

    Default `a_direction` is (1,0,0), unless this is parallel to
    `ab_normal`, in which case default `a_direction` is (0,0,1).

    The returned cell has the vectors va, vb and vc along the rows. The
    cell will be oriented such that va and vb are normal to `ab_normal`
    and va will be along the projection of `a_direction` onto the a-b
    plane.
    Example:
    >>> cell = cellpar_to_cell([1, 2, 4, 10, 20, 30], (0, 1, 1), (1, 2, 3))
    >>> np.round(cell, 3)
    array([[ 0.816, -0.408,  0.408],
            [ 1.992, -0.13 ,  0.13 ],
            [ 3.859, -0.745,  0.745]])
    """
    if a_direction is None:
        if np.linalg.norm(np.cross(ab_normal, (1, 0, 0))) < 1e-5:
            a_direction = (0, 0, 1)
        else:
            a_direction = (1, 0, 0)

    # Define rotated X,Y,Z-system, with Z along ab_normal and X along
    # the projection of a_direction onto the normal plane of Z.
    ad = np.array(a_direction)
    Z = unit_vector(ab_normal)
    X = unit_vector(ad - np.dot(ad, Z) * Z)
    Y = np.cross(Z, X)

    # Express va, vb and vc in the X,Y,Z-system
    alpha, beta, gamma = 90., 90., 90.
    if isinstance(cellpar, (int, float)):
        a = b = c = cellpar
    elif len(cellpar) == 1:
        a = b = c = cellpar[0]
    elif len(cellpar) == 3:
        a, b, c = cellpar
    else:
        a, b, c, alpha, beta, gamma = cellpar

    # Handle orthorhombic cells separately to avoid rounding errors
    eps = 2 * np.spacing(90.0, dtype=np.float64)  # around 1.4e-14
    # alpha
    if abs(abs(alpha) - 90) < eps:
        cos_alpha = 0.0
    else:
        cos_alpha = np.cos(alpha * np.pi / 180.0)
    # beta
    if abs(abs(beta) - 90) < eps:
        cos_beta = 0.0
    else:
        cos_beta = np.cos(beta * np.pi / 180.0)
    # gamma
    if abs(gamma - 90) < eps:
        cos_gamma = 0.0
        sin_gamma = 1.0
    elif abs(gamma + 90) < eps:
        cos_gamma = 0.0
        sin_gamma = -1.0
    else:
        cos_gamma = np.cos(gamma * np.pi / 180.0)
        sin_gamma = np.sin(gamma * np.pi / 180.0)

    # Build the cell vectors
    va = a * np.array([1, 0, 0])
    vb = b * np.array([cos_gamma, sin_gamma, 0])
    cx = cos_beta
    cy = (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    cz_sqr = 1. - cx * cx - cy * cy
    assert cz_sqr >= 0
    cz = np.sqrt(cz_sqr)
    vc = c * np.array([cx, cy, cz])

    # Convert to the Cartesian x,y,z-system
    abc = np.vstack((va, vb, vc))
    T = np.vstack((X, Y, Z))
    cell = np.dot(abc, T)

    return cell


def get_fractional_to_cartesian_matrix(cell_a: float,
                                       cell_b: float,
                                       cell_c: float,
                                       alpha: float,
                                       beta: float,
                                       gamma: float,
                                       angle_in_degrees: bool = True):
    """
    Return the transformation matrix that converts fractional coordinates to
    cartesian coordinates.

    Parameters
    ----------
    a, b, c : float
        The lengths of the edges.
    alpha, gamma, beta : float
        The angles between the sides.
    angle_in_degrees : bool
        True if alpha, beta and gamma are expressed in degrees.
    Returns
    -------
    T_matrix : ndarray
        The 3x3 rotation matrix. ``V_cart = np.dot(T_matrix, V_frac)``.
    """
    if angle_in_degrees:
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)

    cosa = np.cos(alpha)
    cosb = np.cos(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)

    volume = np.sqrt(1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg)

    T_matrix = np.zeros((3, 3))

    T_matrix[0, 0] = cell_a
    T_matrix[0, 1] = cell_b * cosg
    T_matrix[0, 2] = cell_c * cosb
    T_matrix[1, 1] = cell_b * sing
    T_matrix[1, 2] = cell_c * (cosa - cosb * cosg) / sing
    T_matrix[2, 2] = cell_c * volume / sing

    return T_matrix


def get_cartesian_to_fractional_matrix(a: float,
                                       b: float,
                                       c: float,
                                       alpha: float,
                                       beta: float,
                                       gamma: float,
                                       angle_in_degrees: bool = True):
    """
    Return the transformation matrix that converts cartesian coordinates to
    fractional coordinates.
    Parameters
    ----------
    a, b, c : float
        The lengths of the edges.
    alpha, gamma, beta : float
        The angles between the sides.
    angle_in_degrees : bool
        True if alpha, beta and gamma are expressed in degrees.
    Returns
    -------
    T_matrix : array_like
        The 3x3 rotation matrix. ``R_frac = np.dot(T_matrix, R_cart)``.
    """
    if angle_in_degrees:
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)

    cosa = np.cos(alpha)
    cosb = np.cos(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)

    volume = np.sqrt(1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg)

    T_matrix = np.zeros((3, 3))

    T_matrix[0, 0] = 1.0 / a
    T_matrix[0, 1] = -cosg / (a * sing)
    T_matrix[0, 2] = (cosa * cosg - cosb) / (a * volume * sing)
    T_matrix[1, 1] = 1.0 / (b * sing)
    T_matrix[1, 2] = (cosb * cosg - cosa) / (b * volume * sing)
    T_matrix[2, 2] = sing / (c * volume)

    return T_matrix


def get_reciprocal_vectors(cell):
    '''
    Get the reciprocal vectors of a cell given in cell parameters of cell vectors
    ----------
    cell : array
        (3,1) array for cell vectors or (6,1) array for cell parameters
    Returns
    -------
    b1 : array
        (3,1) array containing b_1 vector in the reciprocal space
    b2 : array
        (3,1) array containing b_2 vector in the reciprocal space
    b3 : array
        (3,1) array containing b_3 vector in the reciprocal space
    '''
    if len(cell) == 3:
        v1, v2, v3 = cell
    if len(cell) == 6:
        v1, v2, v3 = cellpar_to_cell(cell)

    vol = np.dot(v1, np.cross(v2, v3))

    b1 = 2 * np.pi * np.cross(v2, v3) / vol
    b2 = 2 * np.pi * np.cross(v3, v1) / vol
    b3 = 2 * np.pi * np.cross(v1, v2) / vol

    return b1, b2, b3


def get_kgrid(cell, distance=0.3):
    '''Get the k-points grid in the reciprocal space with a given distance for a
    cell given in cell parameters of cell vectors.
    ----------
    cell : array
        (3,1) array for cell vectors or (6,1) array for cell parameters
    distance : float
        distance between the points in the reciprocal space
    Returns
    -------
    kx : int
        Number of points in the x direction on reciprocal space
    ky : int
        Number of points in the y direction on reciprocal space
    kz : int
        Number of points in the z direction on reciprocal space
    '''

    b1, b2, b3 = get_reciprocal_vectors(cell)

    b = np.array([np.linalg.norm(b1), np.linalg.norm(b2), np.linalg.norm(b3)])

    kx = np.ceil(b[0]/distance).astype(int)
    ky = np.ceil(b[1]/distance).astype(int)
    kz = np.ceil(b[2]/distance).astype(int)

    return kx, ky, kz


def create_CellBox(A, B, C, alpha, beta, gamma):
    """Creates the CellBox using the same expression as RASPA."""

    tempd = (np.cos(alpha) - np.cos(gamma) * np.cos(beta)) / np.sin(gamma)

    ax = A
    ay = 0
    az = 0
    bx = B * np.cos(gamma)
    by = B * np.sin(gamma)
    bz = 0
    cx = C * np.cos(beta)
    cy = C * tempd
    cz = C * np.sqrt(1 - np.cos(beta) ** 2 - tempd ** 2)

    CellBox = np.array([[ax, ay, az],
                        [bx, by, bz],
                        [cx, cy, cz]])

    return CellBox


def calculate_UnitCells(cell, cutoff):
    '''
    Calculate the number of unit cell repetitions so that all supercell lengths are larger than
    twice the interaction potential cut-off radius.

    RASPA considers the perpendicular directions the directions perpendicular to the `ab`, `bc`,
    and `ca` planes. Thus, the directions depend on who the crystallographic vectors `a`, `b`,
    and `c` are and the length in the perpendicular directions would be the projections
    of the crystallographic vectors on the vectors `a x b`, `b x c`, and `c x a`.
    (here `x` means cross product)
    ----------
    cell_matrix : array
        (3,3) cell vectors or (6,1)
    Returns
    -------
    superCell
        (3,1) list containg the number of repiting units in `x`, `y`, `z` directions.
    '''

    # Make sure that the cell is in the format of cell matrix
    if len(cell) == 6:
        cell_box = cellpar_to_cell(cell)
    if len(cell) == 3:
        cell_box = cell

    # Pre-calculate the cross products
    axb = np.cross(cell_box[0], cell_box[1])
    bxc = np.cross(cell_box[1], cell_box[2])
    cxa = np.cross(cell_box[2], cell_box[0])

    # Calculates the cell volume
    V = np.dot(np.cross(cell_box[0], cell_box[1]), cell_box[2])

    # Calculate perpendicular widths
    cx = V / np.linalg.norm(bxc)
    cy = V / np.linalg.norm(cxa)
    cz = V / np.linalg.norm(axb)

    # Calculate UnitCells array
    supercell = np.ceil(2.0 * cutoff / np.array([cx, cy, cz])).astype(int)

    return supercell


def cellpar_to_lammpsbox(a: float,
                         b: float,
                         c: float,
                         alpha: float,
                         beta: float,
                         gamma: float,
                         angle_in_degrees: bool = True):
    """
    Return the box parameters lx, ly, lz, xy, xz, yz for LAMMPS data input.
    Parameters
    ----------
    a, b, c : float
        The lengths of the edges.
    alpha, gamma, beta : float
        The angles between the sides.
    angle_in_degrees : bool
        True if alpha, beta and gamma are expressed in degrees.
    Returns
    -------
    r : array_like
        The 1x6 array with the box parameters 'lx', 'ly', 'lz', 'xy', 'xz', 'yz'.
    """
    if angle_in_degrees:
        alpha = alpha*(np.pi/180)
        beta = beta*(np.pi/180)
        gamma = gamma*(np.pi/180)

    lx = a
    xy = b * np.cos(gamma)
    xz = c * np.cos(beta)
    ly = np.sqrt(b ** 2 - xy ** 2)
    yz = (b * c * np.cos(alpha) - xy * xz) / ly
    lz = np.sqrt(c ** 2 - xz ** 2 - yz ** 2)

    return np.array([lx, ly, lz, xy, xz, yz])


def find_index(element, e_list):
    '''Finds the index of a given element in a list
    ----------
    element : string
        String containing the label of the element in e_list
    e_list : list
        List with the atom labels
    Returns
    ----------
    i : int
        The index of element in the e_list
    '''
    for i in range(len(e_list)):
        if np.array_equal(e_list[i], element):
            return i


def change_X_atoms(atom_labels, atom_pos, bond_atom):
    ''' Changes the X atom for the desired bond_atom or remove it if bond_atom == 'R'.
    ----------
    atom_labels : list
        List containing the atom labels
    atom_pos : list
        List containing the atom position
    Returns
    ----------
    labels : list
        List containing the processed atom labels
    pos : list
        List containing the processed atom position
    '''
    label, pos = [], []

    for i in range(len(atom_labels)):
        if atom_labels[i] == 'X' and bond_atom != 'R':
            label += [bond_atom]
            pos += [atom_pos[i]]
        if atom_labels[i] != 'X':
            label += [atom_labels[i]]
            pos += [atom_pos[i]]

    return label, pos


def find_bond_atom(cof_name):
    '''Finds the type of atom that the program heve
    to substitute X based on the building blocks'''

    bb1, bb2, net, stacking = cof_name.split('-')
    conect_1 = bb1.split('_')[2]
    conect_2 = bb2.split('_')[2]

    bond_dict = {'NH2': 'N',
                 'CONHNH2': 'N',
                 'BOH2': 'B',
                 'Cl': 'R',
                 'Br': 'R'}

    for group in list(bond_dict.keys()):
        if group in [conect_1, conect_2]:
            return bond_dict[group]


def closest_atom(label_1, pos_1, labels, pos):
    '''Finds the closest atom to a given atom'''

    list_labels = []
    list_pos = []

    for i in range(len(labels)):
        if labels[i] != label_1:
            list_labels += [labels[i]]
            list_pos += [pos[i]]

    closest_index = distance.cdist([pos_1], list_pos).argmin()

    closest_label = list_labels[closest_index]
    closest_position = list_pos[closest_index]
    euclidian_distance = np.linalg.norm(pos_1 - list_pos[closest_index])

    return closest_label, closest_position, euclidian_distance


def closest_atom_struc(label_1, pos_1, labels, pos):
    '''Finds the closest atom on the structure to a given atom'''

    list_labels = []
    list_pos = []
    for i in range(len(labels)):
        if labels[i] != label_1:
            if 'C' in labels[i]:
                list_labels += [labels[i]]
                list_pos += [pos[i]]

    closest_index = distance.cdist([pos_1], list_pos).argmin()

    closet_label = list_labels[closest_index]
    closet_position = list_pos[closest_index]
    euclidian_distance = np.linalg.norm(pos_1 - list_pos[closest_index])

    return closet_label, closet_position, euclidian_distance


def get_bond_atom(connector_1: str, connector_2: str) -> str:
    '''
    Get the atom that will be used to bond two building blocks.
    '''

    bond_dict = {'NH2': 'N',
                 'CONHNH2': 'N',
                 'CHNNH2': 'N',
                 'BOH2': 'B',
                 'Cl': 'C',
                 'Br': 'C',
                 'CHCN': 'C'}

    bond_atom = None
    for group in list(bond_dict.keys()):
        if group in [connector_1, connector_2]:
            bond_atom = bond_dict[group]

    return bond_atom


def print_framework_name(name, lattice, hall, space_group, space_number, symm_op):
    '''Print the results of the created structure'''
    print('{:<60s} {:^12s} {:<4s} {:^4s} #{:^5s}   {:^2} sym. op.'.format(name,
                                                                          lattice,
                                                                          hall.lstrip('-'),
                                                                          space_group,
                                                                          space_number,
                                                                          symm_op))


def print_command(text, verbose, match):
    if verbose in match:
        print(text)


def formula_from_atom_list(AtomLabels: list) -> str:
    """
    Create a string with the formula of the structure from the list of atoms.

    Parameters
    ----------
    AtomLabels : list
        List of strings containing the atom labels.

    Returns
    -------
    formula : str
        String with the formula of the structure.
    """

    formula = ''
    for i in set(AtomLabels):
        formula += i + str(AtomLabels.count(i))

    return formula
