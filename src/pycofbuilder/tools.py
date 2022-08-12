# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020

@author: Felipe Lopes de Oliveira
"""

import os
import numpy as np
from scipy.spatial import distance
try:
    from pymatgen.io.cif import CifParser
except Exception:
    print('WARNING: Could no import CifParser from pymatgen.',
          'The conversion from cif to xyz and COF generation may not work properlly')
    CIF_PARSER_IMPORTED = False
import simplejson


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

    with open(os.path.join(os.path.dirname(__file__), 'data', 'periodic_table.json'), 'r') as f:
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
    r : ndarray
        The 3x3 rotation matrix. ``V_cart = np.dot(r, V_frac)``.
    """
    if angle_in_degrees:
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)
    cosa = np.cos(alpha)
    cosb = np.cos(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)
    volume = 1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
    volume = np.sqrt(volume)
    r = np.zeros((3, 3))
    r[0, 0] = cell_a
    r[0, 1] = cell_b * cosg
    r[0, 2] = cell_c * cosb
    r[1, 1] = cell_b * sing
    r[1, 2] = cell_c * (cosa - cosb * cosg) / sing
    r[2, 2] = cell_c * volume / sing
    return r


def get_cartesian_to_fractional_matrix(a, b, c, alpha, beta, gamma, angle_in_degrees=True):
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
    r : array_like
        The 3x3 rotation matrix. ``V_frac = np.dot(r, V_cart)``.
    """
    if angle_in_degrees:
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)
    cosa = np.cos(alpha)
    cosb = np.cos(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)
    volume = 1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
    volume = np.sqrt(volume)
    r = np.zeros((3, 3))
    r[0, 0] = 1.0 / a
    r[0, 1] = -cosg / (a * sing)
    r[0, 2] = (cosa * cosg - cosb) / (a * volume * sing)
    r[1, 1] = 1.0 / (b * sing)
    r[1, 2] = (cosb * cosg - cosa) / (b * volume * sing)
    r[2, 2] = sing / (c * volume)
    return r


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
    b1 = 2*np.pi*np.cross(v2, v3)/vol
    b2 = 2*np.pi*np.cross(v3, v1)/vol
    b3 = 2*np.pi*np.cross(v1, v2)/vol

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
    kx = np.ceil(b[0]/distance)
    ky = np.ceil(b[1]/distance)
    kz = np.ceil(b[2]/distance)

    return kx, ky, kz


def create_CellBox(A, B, C, alpha, beta, gamma):
    """Creates the CellBox using the same expression as RASPA."""
    tempd = np.cos(alpha) - np.cos(gamma) * np.cos(beta) / np.sin(gamma)
    ax = A
    ay = 0
    az = 0
    bx = B * np.cos(gamma)
    by = B * np.sin(gamma)
    bz = 0 
    cx = C * np.cos(beta)
    cy = C * tempd
    cz = C * np.sqrt(1 - np.cos(beta) ** 2 - tempd ** 2 )
    
    CellBox = np.array([[ax, ay, az], [bx, by, bz], [cx, cy, cz]])
    
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
    SuperCell
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


def cellpar_to_lammpsbox(a, b, c, alpha, beta, gamma, angle_in_degrees=True):
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
    ly = np.sqrt( b**2 - xy **2)
    yz = (b * c * np.cos(alpha) - xy * xz) / ly
    lz = np.sqrt(c**2 - xz**2 - yz**2)

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
    '''Finds the type of atom that the program heve to substitute X based on the building blocks'''
    bb1, bb2, net, stacking = cof_name.split('-')
    conect_1 = bb1.split('_')[2]
    conect_2 = bb2.split('_')[2]

    bond_dict = {'NH2':'N',
                 'CONHNH2':'N',
                 'BOH2': 'B',
                 'Cl': 'R',
                 'Br':'R'}

    for group in list(bond_dict.keys()):
        if group in [conect_1, conect_2]:
            return bond_dict[group]


def closest_atom(label_1, pos_1, labels, pos):

    list_labels = []
    list_pos = []

    for i in range(len(labels)):
        if labels[i] != label_1:
            list_labels += [labels[i]]
            list_pos += [pos[i]]

    closest_index = distance.cdist([pos_1], list_pos).argmin()

    return list_labels[closest_index], list_pos[closest_index], np.linalg.norm(pos_1 - list_pos[closest_index])


def closest_atom_struc(label_1, pos_1, labels, pos):
    list_labels = []
    list_pos = []
    for i in range(len(labels)):
        if labels[i] != label_1:
            if 'C' in labels[i]:
                list_labels += [labels[i]]
                list_pos += [pos[i]]
    closest_index = distance.cdist([pos_1], list_pos).argmin()

    return list_labels[closest_index], list_pos[closest_index], np.linalg.norm(pos_1-list_pos[closest_index])


def print_result(name, lattice, hall, space_group, space_number, symm_op):
    '''Print the results of the created structure'''
    print('{:<60s} {:^12s} {:<4s} {:^4s} #{:^5s}   {:^2} sym. op.'.format(name, lattice, hall.lstrip('-'), space_group, space_number, symm_op))


def print_comand(text, verbose, match):
    if verbose in match:
        print(text)

############# Reads and save files #####################

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
    file_name = file_name.split('.')[0] # Remove the extention if exists
    file_name = os.path.join(path, file_name + '.csv')

    file_temp = open(file_name, 'w')
    if head is not False:
        file_temp.write(head)
    for i in range(len(data)):
        file_temp.write(delimiter.join([str(j) for j in data[i]]) + '\n')

    file_temp.close()


def read_xyz_file(path, file_name):
    """
    Reads a file in format `.xyz` from the `path` given and returns a list containg the N atom labels and 
    a Nx3 array contaning the atoms coordinates.

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
    
    file_name = file_name.split('.')[0] # Remove the extention if exists

    if os.path.exists(os.path.join(path, file_name + '.xyz')):
        temp_file = open(os.path.join(path, file_name + '.xyz'), 'r').readlines()
        
        atoms = [i.split() for i in temp_file[2:]]

        atom_labels = [i[0] for i in atoms if len(i) > 1]
        atom_pos = np.array([[float(i[1]), float(i[2]), float(i[3])] for i in atoms if len(i) > 1])

        return atom_labels, atom_pos
    else:
        print(f'File {file_name} not found!')
        return None


def read_gjf_file(path, file_name):
    """
    Reads a file in format `.gjf` from the `path` given and returns a list containg the N atom labels and 
    a Nx3 array contaning the atoms coordinates.

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

    file_name = file_name.split('.')[0] # Remove the extention if exists

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
    Reads a file in format `.cif` from the `path` given and returns a list containg the N atom labels and 
    a Nx3 array contaning the atoms coordinates.

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

    file_name = file_name.split('.')[0] # Remove the extention if exists

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


def save_xsf(path, file_name, cell, atom_label, atom_pos):
    """
    Save a file in format `.xsf` on the `path`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.xsf` extention. 
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters. 
    atom_label : list
        List of strings containing containg the N atom label. 
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 6:
        cell = cellpar_to_cell(cell)

    xsf_file = open(os.path.join(path, file_name + '.xsf'), 'w')
    xsf_file.write(' CRYSTAL\n')
    xsf_file.write('  PRIMVEC\n')

    for i in range(len(cell)):
        xsf_file.write(f'  {cell[i][0]:>5.9f}    {cell[i][1]:>5.9f}    {cell[i][2]:>5.9f}\n')

    xsf_file.write('   PRIMCOORD\n')
    xsf_file.write(f'           {len(atom_pos)}           1\n')

    for i in range(len(atom_pos)):
        xsf_file.write(f'{atom_label[i]}        {atom_pos[i][0]:>5.9f}    {atom_pos[i][1]:>5.9f}    {atom_pos[i][2]:>5.9f}\n')

    xsf_file.close()


def save_pqr(path, file_name, cell, atom_label, atom_pos, partial_charges=False):
    """
    Save a file in format `.pqr` on the `path`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.pqr` extention. 
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters. 
    atom_label : list
        List of strings containing containg the N atom partial charges. 
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    partial_charges: list
        List of strings containing containg the N atom partial charges. 
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 3:
        cell = cell_to_cellpar(cell)

    pqr_file = open(os.path.join(path, file_name + '.pqr'), 'w')
    pqr_file.write(f'TITLE       {file_name}  \n')
    pqr_file.write('REMARK   4\n')
    pqr_file.write(f'CRYST1{cell[0]:>9.3f}{cell[1]:>9.3f}{cell[2]:>9.3f}{cell[3]:>7.2f}{cell[4]:>7.2f}{cell[5]:>7.2f} P1\n')

    if partial_charges is not False:
        for i in range(len(atom_pos)):
            pqr_file.write(f'ATOM   {i+1:>4} {atom_label[i]:>2}   MOL A   0    {atom_pos[i][0]:>8.3f}{atom_pos[i][1]:>8.3f}{atom_pos[i][2]:>8.3f}{partial_charges:>8.5f}                {atom_label[i]}\n')
    if partial_charges is False:
        for i in range(len(atom_pos)):
            pqr_file.write(f'ATOM   {i+1:>4} {atom_label[i]:>2}   MOL A   0    {atom_pos[i][0]:>8.3f}{atom_pos[i][1]:>8.3f}{atom_pos[i][2]:>8.3f}                {atom_label[i]}\n')

    pqr_file.close()


def save_pdb(path, file_name, cell, atom_label, atom_pos):
    """
    Save a file in format `.pdb` on the `path`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.pdb` extention. 
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters. 
    atom_label : list
        List of strings containing containg the N atom partial charges. 
    atom_pos : list
        Nx3 array contaning the atoms coordinates in cartesian form.
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 3:
        cell = cell_to_cellpar(cell)

    pdb_file = open(os.path.join(path, file_name + '.pdb'), 'w')
    pdb_file.write(f'TITLE       {file_name}  \n')
    pdb_file.write('REMARK   pyCOFBuilder\n')
    pdb_file.write(f'CRYST1{cell[0]:>9.3f}{cell[1]:>9.3f}{cell[2]:>9.3f}{cell[3]:>7.2f}{cell[4]:>7.2f}{cell[5]:>7.2f} P1\n')

    for i in range(len(atom_pos)):
        pdb_file.write(f'ATOM   {i+1:>4} {atom_label[i]:>2}   MOL          {atom_pos[i][0]:>8.3f}{atom_pos[i][1]:>8.3f}{atom_pos[i][2]:>8.3f}  1.00  0.00           {atom_label[i]}\n')

    pdb_file.close()


def save_gjf(path, file_name, atom_labels, atom_pos, text='opt pm6'):
    """
    Save a file in format `.gjf` on the `path`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.gjf` extention. 
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters. 
    atom_label : list
        List of strings containing containg the N atom partial charges. 
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    text : str
        Parameters for Gaussian calculations.
    """

    file_name = file_name.split('.')[0]

    temp_file = open(os.path.join(path, file_name), 'w')
    temp_file.write(f'%chk={file_name}.chk \n')
    temp_file.write(f'# {text}\n')
    temp_file.write('\n')
    temp_file.write(f'{file_name}\n')
    temp_file.write('\n')
    temp_file.write('0 1 \n')

    for i in range(len(atom_labels)):
        temp_file.write('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atom_labels[i],
                                                                     atom_pos[i][0],
                                                                     atom_pos[i][1],
                                                                     atom_pos[i][2]))

    temp_file.write('\n')
    temp_file.write('\n')
    temp_file.close()


def save_xyz(path, file_name, atom_labels, atom_pos, cell=None):
    """
    Save a file in format `.xyz` on the `path`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.xyz` extention. 
    atom_label : list
        List of strings containing containg the N atom partial charges. 
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters. 
    """

    file_name = file_name.split('.')[0]

    temp_file = open(os.path.join(path, file_name + '.xyz'), 'w')
    temp_file.write(f'{len(atom_labels)}\n')

    if cell is None:
        temp_file.write(f'{file_name}\n')
    else:
        if len(cell) == 3:
            cell = cell_to_cellpar(cell)
        temp_file.write(f'{cell[0]}  {cell[1]}  {cell[2]}  {cell[3]}  {cell[4]}  {cell[5]}\n')

    for i in range(len(atom_labels)):
        temp_file.write('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atom_labels[i],
                                                                     atom_pos[i][0],
                                                                     atom_pos[i][1],
                                                                     atom_pos[i][2]))

    temp_file.close()


def save_json(path, file_name, cell, atom_labels, atom_pos):
    """
    Save a file in format `.json` on the `path`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.cif` extention. 
    atom_label : list
        List of strings containing containg the N atom partial charges. 
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters. 
    """

    file_name = file_name.split('.')[0]

    cof_json = create_COF_json(file_name)

    if len(cell) == 3:
        cell_par = cell_to_cellpar(np.array(cell)).tolist()
        cell_par =  [round(i, 10) for i in cell_par]

    if len(cell) == 6:
        cell_par = cell
        cell = cellpar_to_cell(cell_par).tolist()

    cof_json['system']['geo_opt'] = False

    cof_json['geometry']['cell_matrix'] = cell
    cof_json['geometry']['cell_parameters'] = cell_par
    cof_json['geometry']['atom_labels'] = atom_labels
    cof_json['geometry']['atom_pos'] = atom_pos

    write_json(path, file_name, cof_json)


def save_cif(path,
             file_name,
             cell,
             atom_labels,
             atom_pos,
             partial_charges=False,
             frac_coords=True ):
    """
    Save a file in format `.cif` on the `path`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.cif` extention.
    atom_label : list
        List of strings containing containg the N atom partial charges.
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters.
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 3:
        a, b, c, alpha, beta, gamma = cell_to_cellpar(cell)
    if len(cell) == 6:
        a, b, c, alpha, beta, gamma = cell

    cif_text = f"""\
data_{file_name}

_audit_creation_date
_audit_creation_method pyCOFBuilder
_audit_author_name  'Felipe Lopes de Oliveira'

_chemical_name_common                  '{file_name}'
_cell_length_a                          {a:>10.6f}
_cell_length_b                          {b:>10.6f}
_cell_length_c                          {c:>10.6f}
_cell_angle_alpha                       {alpha:>6.2f}
_cell_angle_beta                        {beta:>6.2f}
_cell_angle_gamma                       {gamma:>6.2f}
_space_group_name_H-M_alt               P 1
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

    if partial_charges is not False:
        cif_text += '   _atom_site_charge\n'

    if frac_coords is False:
        # Convert to fractional coordinates
        frac_matrix = get_cartesian_to_fractional_matrix(a, b, c, alpha, beta, gamma)
        atom_pos = [np.dot(frac_matrix, [i[0], i[1], i[2]]) for i in atom_pos]

    for i in range(len(atom_pos)):
        u, v, w = atom_pos[i][0], atom_pos[i][1], atom_pos[i][2]
        if partial_charges is not False:
            cif_text += f'{atom_labels[i]}    {atom_labels[i]} {u:>15.9f} {v:>15.9f} {w:>15.9f} {partial_charges[i]:>10.5f}\n'
        else:
            cif_text += f'{atom_labels[i]}    {atom_labels[i]} {u:>15.9f} {v:>15.9f} {w:>15.9f}\n'

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

    atom_labels, atom_pos = read_gjf_file(path, file_name + '.gjf')

    save_xyz(path, file_name + '.xyz', atom_labels, atom_pos)


def convert_xyz_2_gjf(path, file_name):

    file_name = file_name.split('.')[0]

    atom_labels, atom_pos = read_xyz_file(path, file_name + '.xyz')

    save_xyz(path, file_name + '.gjf', atom_labels, atom_pos)


def convert_cif_2_xyz(path, file_name, supercell=[1, 1, 1]):

    file_name = file_name.split('.')[0]

    if CIF_PARSER_IMPORTED is not False:

        structure = CifParser(os.path.join(path, file_name + '.cif')).get_structures(primitive=True)[0]

        structure.make_supercell([[supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, supercell[2]]])

        dict_sctructure = structure.as_dict()
        a, b, c = dict_sctructure['lattice']['a'], dict_sctructure['lattice']['b'], dict_sctructure['lattice']['c']
        alpha = round(dict_sctructure['lattice']['alpha'])
        beta = round(dict_sctructure['lattice']['beta'])
        gamma = round(dict_sctructure['lattice']['gamma'])

        atom_labels = [i['label'] for i in dict_sctructure['sites']]

        atom_pos = [i['xyz'] for i in dict_sctructure['sites']]

    if CIF_PARSER_IMPORTED is False:
        cell, atom_labels, atom_pos, charges = read_cif(path, file_name)
        a, b, c, alpha, beta, gamma = cell

    temp_file = open(os.path.join(path, file_name + '.xyz'), 'w')
    temp_file.write(f'{len(atom_labels)} \n')

    temp_file.write(f'{a}  {b}  {c}  {alpha}  {beta}  {gamma}\n')

    for i in range(len(atom_labels)):
        temp_file.write('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atom_labels[i], atom_pos[i][0], atom_pos[i][1], atom_pos[i][2]))

    temp_file.close()

########################### JSON related ##########################  

def write_json(path, name, COF_json):

    name = name.split('.')[0]

    if os.path.exists(path) is not True:
        os.mkdir(path)

    save_path = os.path.join(path, name + '.json')

    with open(save_path, 'w', encoding='utf-8') as f:
        simplejson.dump(COF_json, f, ensure_ascii=False, separators=(',', ':'), indent=2, ignore_nan=True)


def read_json(path, cof_name):

    cof_path = os.path.join(path, cof_name + '.json')

    with open(cof_path, 'r') as r:
        json_object = simplejson.loads(r.read())

    return json_object


def create_COF_json(name):

    system_info = 'Informations about the system such as name, if it is optimized and other relevant information.'
    geometry_info = 'Informations about the geometry: cell parameters, cell matrix, atomic positions, partial charges, bond orders, simmetry information'
    optimization_info = 'Information about the optimization process such as level of calculations, optimization schema and optimization steps.'
    adsorption_info = 'Information about the adsorption simulation experiments on RASPA2'
    textural_info = 'Information about the textural calculations of the structure such as specific area, pore volume, void fraction.'
    spectrum_info = 'Information about spectra simulation like DRX, FTIR, ssNMR, UV-VIS, Band dispersion, Phonon dispersion...'
    experimental_info = 'Experimental data DRX, FTIR, ssNMR, UV-VIS...'

    COF_json = {'system':{'description': system_info,
                          'name': name,
                          'geo_opt': False,
                          'execution_times_seconds': {}},
                'geometry':{'description':geometry_info},
                'optimization':{'description':optimization_info},
                'adsorption':{'description':adsorption_info},
                'textural':{'description':textural_info},
                'spectrum':{'description':spectrum_info},
                'experimental':{'description':experimental_info}
                }

    return COF_json
