# -*- coding: utf-8 -*-
# Created by Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

"""
This module contains the tools used by pyCOFBuilder.
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
    norm = y / np.linalg.norm(y)
    return norm


def angle(v1, v2, unit='degree'):
    """
    Calculates the angle between two vectors v1 and v2.

    Parameters
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

    Parameters
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
        kmat = np.array([[0, -v[2], v[1]],
                         [v[2], 0, -v[0]],
                         [-v[1], v[0], 0]])

        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
        return rotation_matrix
    else:
        return np.identity(3)


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
    given a 3x3 cell matrix.

    Angles are in degrees unless radian=True is used.

    Parameters
    ----------
    cell : array
        (3,3) matrix of cell vectors v1, v2, and v3
    radians : bool
        Return the cell angles in radians

    Returns
    -------
    cellpar : array
        (6,1) vector with the cell parameters
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

    # Corvet to radians if radians is True
    if radians:
        angles = [angle * np.pi / 180 for angle in angles]

    return np.array(lengths + angles)


def cellpar_to_cell(cellpar, ab_normal=(0, 0, 1), a_direction=None):
    """Return a 3x3 cell matrix from cell parameters (a,b,c,alpha,beta, and gamma).

    Angles must be in degrees.

    The returned cell is orientated such that a and b are normal to `ab_normal` and a is
    parallel to the projection of `a_direction` in the a-b plane.

    Default `a_direction` is (1,0,0), unless this is parallel to
    `ab_normal`, in which case default `a_direction` is (0,0,1).

    The returned cell has the vectors va, vb and vc along the rows. The
    cell will be oriented such that va and vb are normal to `ab_normal`
    and va will be along the projection of `a_direction` onto the a-b
    plane.

    Parameters
    ----------
    cellpar : array
        (6,1) vector with the cell parameters
    ab_normal : array
        Normal vector between a and b cell vectors. Default: (0, 0, 1)
    a_direction : array
        Specific direction for the a vector. Default: None
    Returns
    -------
    cell : array
        (3,3) matrix of cell vectors v1, v2, and v3
    """
    if a_direction is None:
        if np.linalg.norm(np.cross(ab_normal, (1, 0, 0))) < 1e-5:
            a_direction = (0, 0, 1)
        else:
            a_direction = (1, 0, 0)

    # Define rotated X,Y,Z-system, with Z along ab_normal and X along the
    # projection of a_direction onto the normal plane of Z.
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
                                       angle_in_degrees: bool = True) -> np.array(float):
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
                                       angle_in_degrees: bool = True) -> np.array(float):
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
    T_matrix : np.array
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


def get_reciprocal_vectors(cell) -> tuple:
    """
    Get the reciprocal vectors of a cell given in cell parameters of cell vectors.

    Parameters
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
    """

    if len(cell) == 3:
        v1, v2, v3 = cell
    if len(cell) == 6:
        v1, v2, v3 = cellpar_to_cell(cell)

    vol = np.dot(v1, np.cross(v2, v3))

    b1 = 2 * np.pi * np.cross(v2, v3) / vol
    b2 = 2 * np.pi * np.cross(v3, v1) / vol
    b3 = 2 * np.pi * np.cross(v1, v2) / vol

    return b1, b2, b3


def get_kgrid(cell, distance=0.3) -> tuple:
    """
    Get the k-points grid in the reciprocal space with a given distance for a
    cell given in cell parameters of cell vectors.

    Parameters
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
    """

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

    Parameters
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
    """
    Finds the index of a given element in a list

    Parameters
    ----------
    element : string
        String containing the label of the element in e_list
    e_list : list
        List with the atom labels
    Returns
    ----------
    i : int
        The index of element in the e_list
    """

    index = None
    for i in range(len(e_list)):
        if np.array_equal(e_list[i], element):
            index = i
            break
    return index


def change_X_atoms(atom_labels, atom_pos, bond_atom) -> tuple:
    '''
    Changes the X atom for the desired bond_atom or remove it if bond_atom == 'R'.

    Parameters
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


def closest_atom(label_1: str, pos_1: list, labels: list, pos: list):
    '''
    Find the closest atom to a given atom

    Parameters
    ----------
    label_1 : string
        String containing the label of the atom
    pos_1 : list
        Array containing the position of the atom
    labels : list
        List containing the all the atom labels on the structure
    pos : list
        List containing the all the atom positions on the structure

    Returns
    ----------
    closest_label : string
        String containing the label of the closest atom
    closest_position : array
        Array containing the position of the closest atom
    euclidian_distance : float
        Euclidian distance between the two atoms
    '''

    list_labels = []
    list_pos = []

    for i in range(len(labels)):
        if labels[i] != label_1:
            list_labels += [labels[i]]
            list_pos += [pos[i]]

    if len(list_pos) == 0:
        return None, np.array([0, 0, 0]), None

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

    bond_dict = {
        'NH2': 'N',
        'NHOH': 'N',
        'COCHCHOH': 'N',
        'CONHNH2': 'N',
        'CHNNH2': 'N',
        'COOH': 'N',
        'BOH2': 'B',
        'OH2': 'B',
        'Cl': 'X',
        'Br': 'X',
        'CH2CN': 'C',
        'CH3': 'C'
    }

    bond_atom = None
    for group in list(bond_dict.keys()):
        if group in [connector_1, connector_2]:
            bond_atom = bond_dict[group]

    return bond_atom


def get_framework_symm_text(name, lattice, hall, space_group, space_number, symm_op):
    '''Get the text for the framework symmop'''
    text = '{:<60s} {:^12s} {:<4s} {:^4s} #{:^5s}   {:^2} sym. op.'.format(name,
                                                                           lattice,
                                                                           hall.lstrip('-'),
                                                                           space_group,
                                                                           space_number,
                                                                           symm_op)
    return text


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


def smiles_to_xsmiles(smiles_string: str) -> str:
    '''
    Converts a SMILES string to an extended SMILES string with labels

    Parameters
    ----------
    smiles_string : str
        SMILES string to be converted

    Returns
    -------
    xsmiles : str
        Extended SMILES string with special labels
    xsmiles_label : str
        xsmiles labels for images with the special labels
    composition : str
        String containing the composition
    '''
    SPECIAL_ATOMS = ['Q', 'R', 'X']
    REGULAR_ATOMS = ['C', 'N', 'H', 'O', 'S', 'B', 'F']

    xsmiles = ''
    labels = []
    atom_list = []

    for i, letter in enumerate(smiles_string):

        if letter in SPECIAL_ATOMS:
            xsmiles += '*'
            labels.append(letter)
            if letter == 'R':
                atom_list.append(smiles_string[i:i+2])
            else:
                atom_list.append(letter)

        elif letter.isnumeric():
            if smiles_string[i-1] == 'R':
                labels[-1] = labels[-1] + letter
            else:
                xsmiles += letter

        elif letter in REGULAR_ATOMS:
            xsmiles += letter
            labels += ['']
            atom_list.append(letter)

        else:
            xsmiles += letter

    # Generate the xsmiles label
    xsmiles_label = '|$' + ';'.join(labels) + '$|'

    # Generate the composition
    composition = formula_from_atom_list(atom_list)

    return xsmiles, xsmiles_label, composition


def ibrav_to_cell(ibrav, celldm1, celldm2, celldm3, celldm4, celldm5, celldm6):
    """
    Convert a value of ibrav to a cell.

    Parameters
    ----------
    ibrav : int
    celldmx: float

    Returns
    -------
    cell : matrix
        The cell as a 3x3 numpy array
    """

    alat = celldm1 * 0.5291772105638411

    if ibrav == 1:
        cell = np.identity(3) * alat
    elif ibrav == 2:
        cell = np.array([[-1.0, 0.0, 1.0],
                         [0.0, 1.0, 1.0],
                         [-1.0, 1.0, 0.0]]) * (alat / 2)
    elif ibrav == 3:
        cell = np.array([[1.0, 1.0, 1.0],
                         [-1.0, 1.0, 1.0],
                         [-1.0, -1.0, 1.0]]) * (alat / 2)
    elif ibrav == -3:
        cell = np.array([[-1.0, 1.0, 1.0],
                         [1.0, -1.0, 1.0],
                         [1.0, 1.0, -1.0]]) * (alat / 2)
    elif ibrav == 4:
        cell = np.array([[1.0, 0.0, 0.0],
                         [-0.5, 0.5 * 3**0.5, 0.0],
                         [0.0, 0.0, celldm3]]) * alat
    elif ibrav == 5:
        tx = ((1.0 - celldm4) / 2.0)**0.5
        ty = ((1.0 - celldm4) / 6.0)**0.5
        tz = ((1 + 2 * celldm4) / 3.0)**0.5
        cell = np.array([[tx, -ty, tz],
                         [0, 2 * ty, tz],
                         [-tx, -ty, tz]]) * alat
    elif ibrav == -5:
        ty = ((1.0 - celldm4) / 6.0)**0.5
        tz = ((1 + 2 * celldm4) / 3.0)**0.5
        a_prime = alat / 3**0.5
        u = tz - 2 * 2**0.5 * ty
        v = tz + 2**0.5 * ty
        cell = np.array([[u, v, v],
                         [v, u, v],
                         [v, v, u]]) * a_prime
    elif ibrav == 6:
        cell = np.array([[1.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.0, 0.0, celldm3]]) * alat
    elif ibrav == 7:
        cell = np.array([[1.0, -1.0, celldm3],
                         [1.0, 1.0, celldm3],
                         [-1.0, -1.0, celldm3]]) * (alat / 2)
    elif ibrav == 8:
        cell = np.array([[1.0, 0.0, 0.0],
                         [0.0, celldm2, 0.0],
                         [0.0, 0.0, celldm3]]) * alat
    elif ibrav == 9:
        cell = np.array([[1.0 / 2.0, celldm2 / 2.0, 0.0],
                         [-1.0 / 2.0, celldm2 / 2.0, 0.0],
                         [0.0, 0.0, celldm3]]) * alat
    elif ibrav == -9:
        cell = np.array([[1.0 / 2.0, -celldm2 / 2.0, 0.0],
                         [1.0 / 2.0, celldm2 / 2.0, 0.0],
                         [0.0, 0.0, celldm3]]) * alat
    elif ibrav == 10:
        cell = np.array([[1.0 / 2.0, 0.0, celldm3 / 2.0],
                         [1.0 / 2.0, celldm2 / 2.0, 0.0],
                         [0.0, celldm2 / 2.0, celldm3 / 2.0]]) * alat
    elif ibrav == 11:
        cell = np.array([[1.0 / 2.0, celldm2 / 2.0, celldm3 / 2.0],
                         [-1.0 / 2.0, celldm2 / 2.0, celldm3 / 2.0],
                         [-1.0 / 2.0, -celldm2 / 2.0, celldm3 / 2.0]]) * alat
    elif ibrav == 12:
        sinab = (1.0 - celldm4**2)**0.5
        cell = np.array([[1.0, 0.0, 0.0],
                         [celldm2 * celldm4, celldm2 * sinab, 0.0],
                         [0.0, 0.0, celldm3]]) * alat
    elif ibrav == -12:
        sinac = (1.0 - celldm5**2)**0.5
        cell = np.array([[1.0, 0.0, 0.0],
                         [0.0, celldm2, 0.0],
                         [celldm3 * celldm5, 0.0, celldm3 * sinac]]) * alat
    elif ibrav == 13:
        sinab = (1.0 - celldm4**2)**0.5
        cell = np.array([[1.0 / 2.0, 0.0, -celldm3 / 2.0],
                         [celldm2 * celldm4, celldm2 * sinab, 0.0],
                         [1.0 / 2.0, 0.0, celldm3 / 2.0]]) * alat
    elif ibrav == 14:
        sinab = (1.0 - celldm4**2)**0.5
        v3 = [celldm3 * celldm5,
              celldm3 * (celldm6 - celldm5 * celldm4) / sinab,
              celldm3 * ((1 + 2 * celldm6 * celldm5 * celldm4
                          - celldm6**2 - celldm5**2 - celldm4**2)**0.5) / sinab]
        cell = np.array([[1.0, 0.0, 0.0],
                         [celldm2 * celldm4, celldm2 * sinab, 0.0],
                         v3]) * alat
    else:
        raise NotImplementedError('ibrav = {0} is not implemented'.format(ibrav))

    return cell


def equal_value(val1, val2, threshold=1e-3) -> bool:
    '''
    Determine if two values are equal based on a given threshold.
    '''
    return abs(val1 - val2) <= threshold


def classify_unit_cell(cell, thr=1e-3) -> str:
    '''
    Determine the bravais lattice based on the cell lattice.
    The cell lattice can be the cell parameters as (6,1) array or
    the cell vectors as (3x3) array.

    Bravais lattice can be cubic, tetragonal, orthorhombic, hexagonal,
    monoclinic, or triclinic.

    Parameters
    ----------
    cell : array
        Array with the cell vectors or parameters
    threshold: float
        Numeric threshold for the analysis. Default: 1e-3

    Returns
    -------
    cell_type : string
        Bravais lattice.
    '''

    if len(cell) == 3:
        a, b, c, alpha, beta, gamma = cell_to_cellpar(cell)

    cell_type = None

    if equal_value(alpha, 90, thr) and equal_value(beta, 90, thr) and equal_value(gamma, 90, thr):
        if equal_value(a, b, thr) and equal_value(b, c, thr):
            cell_type = "cubic"
        if equal_value(a, b, thr) and not equal_value(a, c, thr):
            cell_type = "tetragonal"
        else:
            cell_type = "orthorhombic"
    elif equal_value(alpha, 90, thr) and equal_value(beta, 90, thr) and equal_value(gamma, 120, thr):
        if equal_value(a, b, thr):
            cell_type = "hexagonal"
    elif equal_value(alpha, 90, thr) or equal_value(beta, 90, thr) or equal_value(gamma, 90, thr):
        if not equal_value(a, b, thr) and not equal_value(b, c, thr) and not equal_value(a, c, thr):
            cell_type = "monoclinic"
    else:
        cell_type = "triclinic"

    return cell_type


def cell_to_ibrav(cell):
    '''
    Return the ibrav number for a given cell.
    '''

    if len(cell) == 3:
        a, b, c, alpha, beta, gamma = cell_to_cellpar(cell)
    else:
        a, b, c, alpha, beta, gamma = cell

    cell_type = classify_unit_cell(cell)

    if cell_type == 'cubic':
        celldm = {'ibrav': 1,
                  'celldm(1)': a / 0.5291772105638411}
    elif cell_type == 'hexagonal':
        celldm = {'ibrav': 4,
                  'celldm(1)': a / 0.5291772105638411,
                  'celldm(3)': c / a}
    elif cell_type == 'tetragonal':
        celldm = {'ibrav': 6,
                  'celldm(1)': a / 0.5291772105638411,
                  'celldm(3)': c / a}
    elif cell_type == 'orthorhombic':
        celldm = {'ibrav': 8,
                  'celldm(1)': a / 0.5291772105638411,
                  'celldm(2)': b / a,
                  'celldm(3)': c / a}
    elif cell_type == 'monoclinic':
        celldm = {'ibrav': 12,
                  'celldm(1)': a / 0.5291772105638411,
                  'celldm(2)': b / a,
                  'celldm(3)': c / a,
                  'celldm(4)': np.cos(np.deg2rad(beta))}
    else:
        celldm = {'ibrav': 14,
                  'celldm(1)': a / 0.5291772105638411,
                  'celldm(2)': b / a,
                  'celldm(3)': c / a,
                  'celldm(4)': np.cos(np.deg2rad(alpha)),
                  'celldm(5)': np.cos(np.deg2rad(beta)),
                  'celldm(6)': np.cos(np.deg2rad(gamma))}

    return celldm


def is_bonded(atom1: str, atom2: str, dist: float, cutoff: float = 1.3):
    """
    Determine if two atoms are bonded based on the distance between them.
    Two atoms are considered bonded if the distance between them is less than
    the sum of their covalent radii multiplied by a cutoff factor.

    Parameters
    ----------
    atom1 : str
        Label of the first atom
    atom2 : str
        Label of the second atom
    dist : float
        Distance between the two atoms
    cutoff : float
        Cutoff factor for the covalent radii. Default: 1.3
    """

    periodic_table = elements_dict(property='covalent_radius')

    # Get the covalent radii of the two atoms
    cr_1 = periodic_table[atom1]
    cr_2 = periodic_table[atom2]

    # Calculate max bond distance
    max_bond_distance = (cr_1 + cr_2) * cutoff

    if dist < 0.6:
        print('Distance between atoms is less than 0.6 Å. Check if the structure is correct.')

    # Check if the distance is less than the cutoff
    if 0.6 < dist <= max_bond_distance:
        return True
    else:
        return False


def get_bonds(structure, cutoff=1.3):
    """
    Get the bonded atoms in a structure based on the distance between them.

    Parameters
    ----------
    structure : pymatgen.Structure
        Structure object of pymatgen
    cutoff : float
        Cutoff factor for the covalent radii. Default: 1.3 Å
    """

    atom_types = [i.as_dict()['species'][0]['element'] for i in structure.sites]

    # Get bonded atoms
    center_indices, points_indies, _, bond_distances = structure.get_neighbor_list(5)
    bonded_atoms = np.array([center_indices, points_indies, bond_distances]).T

    bonded_atoms = [i for i in bonded_atoms if is_bonded(atom_types[int(i[0])], atom_types[int(i[1])], i[2], cutoff)]

    bonded_atoms = [[int(i[0]), int(i[1]), i[2]] for i in bonded_atoms]

    return bonded_atoms
