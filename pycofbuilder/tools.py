# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020

@author: lipel
"""

import os
import numpy as np
import math
try:
    from pymatgen.io.cif import CifParser
except Exception:
    print('WARNING: Could no import CifParser from pymatgen the conversion from cif to xyz and COF generation may not work properlly')
    cif_parser_imported = False

def elements_dict():

    element_symbols = [
    'H', 'He',  # Period 1
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',  # Period 2
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',  # Period 3
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',  # Period 4
    'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',  # Period 4
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',  # Period 5
    'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',  # Period 5
    'Cs', 'Ba',  # Period 6
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',  # Lanthanides
    'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',  # Lanthanides
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',  # Period 6
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',  # Period 6
    'Fr', 'Ra',  # Period 7
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',  # Actinides
    'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',  # Actinides
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',  # Period 7
    'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo', # Period 7
    'X', 'Q', 'R', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6']  #Specific labels for the programm

    """Returns a dictionary containing the elements and its respective atomic mass in g/mol"""
    return {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.0107, 'N': 14.0067,
              'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
              'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
              'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045, 'Fe': 55.845, 
              'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64,
              'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585,
              'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
              'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6,
              'I': 126.90447, 'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116,
              'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
              'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04,
              'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
              'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804,
              'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0254, 'Ac': 227.0278,
              'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 'Am': 243.0614,
              'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0951,
              'No': 259.1009, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271, 'Bh': 270, 'Hs': 269, 'Mt': 278,
              'Ds': 281, 'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294,
              'X': 0.0, 'Q': 0.0, 'R': 0.0, 'R1': 0.0, 'R2': 0.0, 'R3': 0.0, 'R4': 0.0, 'R5': 0.0, 'R6': 0.0}

def unit_vector(x):
    """Return a unit vector in the same direction as x."""
    y = np.array(x, dtype='float')
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
    """Returns the cell parameters [a, b, c, alpha, beta, gamma].

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

def get_fractional_to_cartesian_matrix(a, b, c, alpha, beta, gamma, angle_in_degrees=True):
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
    r : array_like
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
    r[0, 0] = a
    r[0, 1] = b * cosg
    r[0, 2] = c * cosb
    r[1, 1] = b * sing
    r[1, 2] = c * (cosa - cosb * cosg) / sing
    r[2, 2] = c * volume / sing
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
    kx = math.ceil(b[0]/distance)
    ky = math.ceil(b[1]/distance)
    kz = math.ceil(b[2]/distance)

    return kx, ky, kz

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
    label = []
    pos = []

    for i in range(len(atom_labels)):
        if atom_labels[i] == 'X' and bond_atom != 'R':
            label += [bond_atom]
            pos += [atom_pos[i]]
        if atom_labels[i] != 'X':
            label += [atom_labels[i]]
            pos += [atom_pos[i]]

    return label, pos 

def find_bond_atom(cof_name):
    bb1, bb2, net, stacking = cof_name.split('-')
    conect_1 = bb1.split('_')[2]
    conect_2 = bb2.split('_')[2]

    if 'NH2' in [conect_1, conect_2]:
        return 'N'
    if 'B(OH)2' in [conect_1, conect_2]:
        return 'B'
    if 'Cl' in [conect_1, conect_2]:
        return 'R'
    if 'Br' in [conect_1, conect_2]:
        return 'R'

def print_result(name, lattice, hall, space_group, space_number, symm_op):
    '''Print the results of the created structure'''
    print('{:<60s} {:^12s} {:<4s} {:^4s} #{:^5s}   {:^2} sym. op.'.format(name, lattice, hall.lstrip('-'), space_group, space_number, symm_op))

############# Reads and save files #####################

def read_xyz_file(path, file_name):
    '''Lê um arquivo em formato .xyz e retorna uma lista com os átomos e um array Nx3 contendo as coordenadas dos N átomos'''
    
    file_name = file_name.split('.')[0]

    if os.path.exists(os.path.join(path, file_name + '.xyz')):
        temp_file = open(os.path.join(path, file_name + '.xyz'), 'r').readlines()

        n_atoms = int(temp_file[0].rstrip('\n'))
        atoms = [i.split() for i in temp_file[2:]]

        atom_labels = [i[0] for i in atoms if len(i) > 1]
        atom_pos = np.array([[float(i[1]), float(i[2]), float(i[3])] for i in atoms if len(i) > 1])

        connectivity = len([i for i in atom_labels if 'X' in i])

        if connectivity == 0:
            print('Non X point could be found!')

        return atom_labels, atom_pos, n_atoms, connectivity
    else:
        print(f'File {file_name} not found!')

def read_gjf_file(path, file_name):
    '''Lê um arquivo em formato .gjf e retorna uma lista com os átomos e um array Nx3 contendo as coordenadas dos N átomos'''

    temp_file = open(os.path.join(path, file_name + '.gjf'), 'r').readlines()
    temp_file = [i.split() for i in temp_file if i != '\n']

    atoms = [i for i in temp_file if i[0] in elements_dict()]

    atom_labels = [i[0] for i in atoms]
    atom_pos = np.array([[float(i[1]), float(i[2]), float(i[3])] for i in atoms])

    n_atoms = len(atom_labels)
    connectivity = len([i for i in atom_labels if 'X' in i])

    return atom_labels, atom_pos

def convert_gjf_2_xyz(path, file_name):

    file_name = file_name.split('.')[0]

    atom_labels, atom_pos = read_gjf_file(path, file_name + '.gjf')

    save_xyz(path, file_name + '.xyz', atom_labels, atom_pos)

def convert_xyz_2_gjf(path, file_name):

    file_name = file_name.split('.')[0]

    atom_labels, atom_pos, n_atoms, connectivity = read_xyz_file(path, file_name + '.xyz')

    save_xyz(path, file_name + '.gjf', atom_labels, atom_pos)

def convert_cif_2_xyz(path, file_name, supercell=[1, 1, 1]):

    file_name = file_name.split('.')[0]

    if cif_parser_imported is not False:

        structure = CifParser(os.path.join(path, file_name + '.cif')).get_structures(primitive=True)[0]

        structure.make_supercell([[supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, supercell[2]]])

        dict_sctructure = structure.as_dict()
        a, b, c = dict_sctructure['lattice']['a'], dict_sctructure['lattice']['b'], dict_sctructure['lattice']['c']
        alpha = round(dict_sctructure['lattice']['alpha'])
        beta = round(dict_sctructure['lattice']['beta'])
        gamma = round(dict_sctructure['lattice']['gamma'])

        atom_labels = [i['label'] for i in dict_sctructure['sites']]

        atom_pos = [i['xyz'] for i in dict_sctructure['sites']]

    if cif_parser_imported is False:
        cell, atom_labels, atom_pos, charges = read_cif(path, file_name)
        a, b, c, alpha, beta, gamma = cell

    temp_file = open(os.path.join(path, file_name + '.xyz'), 'w')
    temp_file.write(f'{len(atom_labels)} \n')

    temp_file.write(f'{a}  {b}  {c}  {alpha}  {beta}  {gamma}\n')

    for i in range(len(atom_labels)):
        temp_file.write('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atom_labels[i], atom_pos[i][0], atom_pos[i][1], atom_pos[i][2]))

    temp_file.close()


def save_csv(file_name, data, delimiter=';', head=False):

    file_temp = open(file_name, 'w')
    if head is not False:
        file_temp.write(head)
    for i in range(len(data)):
        file_temp.write(delimiter.join([str(j) for j in data[i]]) + '\n')

    file_temp.close()

def save_xsf(file_path, file_name, cell, atom_pos):

    file_name = file_name.split('.')[0]

    xsf_file = open(os.path.join(file_path, file_name + '.xsf'), 'w')
    xsf_file.write(' CRYSTAL\n')
    xsf_file.write('  PRIMVEC\n')

    for i in range(len(cell)):
        xsf_file.write(f'  {cell[i][0]:<.9f}    {cell[i][1]:<.9f}    {cell[i][2]:<.9f}\n')

    xsf_file.write('   PRIMCOORD\n')
    xsf_file.write(f'           {len(atom_pos)}           1\n')

    for i in range(len(atom_pos)):
        xsf_file.write(f'{atom_pos[i][0]}        {atom_pos[i][1]:<.9f}    {atom_pos[i][2]:<.9f}    {atom_pos[i][3]:<.9f}\n')

    xsf_file.close()

def save_gjf(file_path, file_name, atom_labels, atom_pos, text='opt pm6'):

    file_name = file_name.split('.')[0]

    temp_file = open(os.path.join(file_path, file_name), 'w')
    temp_file.write(f'%chk={file_name[:-4]}.chk \n')
    temp_file.write(f'# {text}\n')
    temp_file.write('\n')
    temp_file.write(f'{file_name}\n')
    temp_file.write('\n')
    temp_file.write('0 1 \n')

    for i in range(len(atom_labels)):
        temp_file.write('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atom_labels[i], atom_pos[i][0], atom_pos[i][1], atom_pos[i][2]))

    temp_file.write('\n')
    temp_file.write('\n')
    temp_file.close()

def save_xyz(file_path, file_name, atom_labels, atom_pos):

    file_name = file_name.split('.')[0]

    temp_file = open(os.path.join(file_path, file_name + '.xyz'), 'w')
    temp_file.write(f'{len(atom_labels)} \n')
    temp_file.write(f'{file_name[:-4]} rotated \n')

    for i in range(len(atom_labels)):
        temp_file.write('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atom_labels[i], atom_pos[i][0], atom_pos[i][1], atom_pos[i][2]))

    temp_file.close()
    
def read_cif(file_path, file_name):
    tmp = open(os.path.join(file_path, file_name), 'r').readlines()
    cell = []
    atom_label = []
    atom_pos = []
    charges = []
    has_charges = False
    for i in tmp:
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

    for i in tmp:
        line = i.split()
        if len(line) > 1 and line[0] in elements_dict().keys():
            atom_label += [line[0]]
            atom_pos += [[float(j) for j in line[2:-1]]]
            charges += [float(line[-1])]
    cell = cellpar_to_cell(cell)

    return cell, atom_label, atom_pos, charges
        

def save_cif(file_path, file_name, cell, atom_labels, atom_pos, partial_charges=False, frac_coords=True ):

    file_name = file_name.split('.')[0]

    if len(cell) == 3:
        a, b, c, alpha, beta, gamma = cell_to_cellpar(cell)
    if len(cell) == 6:
        a, b, c, alpha, beta, gamma = cell

    r = get_cartesian_to_fractional_matrix(a, b, c, alpha, beta, gamma)

    cif_file = open(os.path.join(file_path, file_name + '.cif'), 'w')

    cif_file.write(f'data_{file_name}\n')
    cif_file.write(f'_chemical_name_common                  \'{file_name}\'\n')
    cif_file.write(f'_cell_length_a                         {a:.6f}\n')
    cif_file.write(f'_cell_length_b                         {b:.6f}\n')
    cif_file.write(f'_cell_length_c                         {c:.6f}\n')
    cif_file.write(f'_cell_angle_alpha                      {alpha:.2f}\n')
    cif_file.write(f'_cell_angle_beta                       {beta:.2f}\n')
    cif_file.write(f'_cell_angle_gamma                      {gamma:.2f}\n')
    cif_file.write('_space_group_name_H-M_alt              \'P 1\'\n')
    cif_file.write('_space_group_IT_number                 1\n')
    cif_file.write('\n')
    cif_file.write('loop_\n')
    cif_file.write('_symmetry_equiv_pos_as_xyz\n')
    cif_file.write('   \'x, y, z\'\n')
    cif_file.write('\n')
    cif_file.write('loop_\n')
    cif_file.write('   _atom_site_label\n')
    cif_file.write('   _atom_site_type_symbol\n')
    cif_file.write('   _atom_site_fract_x\n')
    cif_file.write('   _atom_site_fract_y\n')
    cif_file.write('   _atom_site_fract_z\n')
    if partial_charges is not False:
        cif_file.write('   _atom_site_charge\n')

    if frac_coords == False:
        atom_pos = [np.dot(r, [i[0], i[1], i[2]]) for i in atom_pos]

    for i in range(len(atom_pos)):
        u, v, w = atom_pos[i][0], atom_pos[i][1], atom_pos[i][2]
        if partial_charges is not False:
            cif_file.write(f'{atom_labels[i]}    {atom_labels[i]}    {u:<.9f}    {v:<.9f}    {w:<.9f}   {partial_charges[i]:.5f}\n')
        else:
            cif_file.write(f'{atom_labels[i]}    {atom_labels[i]}    {u:<.9f}    {v:<.9f}    {w:<.9f}\n')

    cif_file.close()
    

if __name__ == '__main__':
    main()