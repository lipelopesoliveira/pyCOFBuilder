from typing import Sequence

import numpy as np
from numpy.typing import NDArray


def read_mol_file(
    file_path,
) -> tuple[list, list[list[float]], list[float], list[tuple[int, int]], list[int]]:
    """
    Reads a .mol file and extracts atom types, Cartesian positions, partial charges, bonds, and bond types.

    Parameters
    ----------
    file_path : str
        Path to the .mol file.
    Returns
    -------
    atomTypes : list
        List of atom types.
    cartPos : list of list of float
        List of Cartesian positions for each atom.
    partialCharges : list of float
        List of partial charges for each atom.
    bonds : list of tuple of int
        List of bonds represented as tuples of atom indices.
    bondTypes : list of int
        List of bond types.
    """

    with open(file_path, "r") as f:
        mol_data = f.read().splitlines()

    n_atoms, n_bonds = map(int, mol_data[3].split()[:2])
    atoms_non_processed = mol_data[4 : 4 + n_atoms]

    atomTypes, cartPos, partialCharges = [], [], []

    for atom_line in atoms_non_processed:
        atom_line = atom_line.split()
        atomTypes.append(atom_line[3])
        cartPos.append(np.array(atom_line[:3], dtype=float))
        partialCharges.append(float(atom_line[4]))

    # Replace "F" by "R" on atomTypes
    bonds_non_processed = mol_data[4 + n_atoms : 4 + n_atoms + n_bonds]

    bonds: list[tuple[int, int]] = [
        (int(bond_line.split()[0]), int(bond_line.split()[1]))
        for bond_line in bonds_non_processed
    ]
    bondTypes = [int(bond_line.split()[2]) for bond_line in bonds_non_processed]

    return atomTypes, np.array(cartPos).tolist(), partialCharges, bonds, bondTypes


def rotation_matrix_from_vectors(vec1: NDArray, vec2: NDArray) -> NDArray:
    """
    Calculates the rotation matrix that aligns vec1 to vec2.

    Parameters
    ----------
    vec1 : NDArray
        Source vector (shape: 3,).
    vec2 : NDArray
        Destination vector (shape: 3,).

    Returns
    -------
    NDArray
        A 3x3 rotation matrix.
    """
    # Normalize vectors
    norm_v1 = np.linalg.norm(vec1)
    norm_v2 = np.linalg.norm(vec2)

    if np.isclose(norm_v1, 0) or np.isclose(norm_v2, 0):
        raise ValueError("Input vectors must have non-zero magnitude.")

    a = (vec1 / norm_v1).reshape(3)
    b = (vec2 / norm_v2).reshape(3)

    if np.allclose(a, b):
        return np.eye(3)
    if np.allclose(a, -b):
        # 180-degree rotation case. Returns -I (inversion) as a simplification.
        return -np.eye(3)

    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)

    # Skew-symmetric cross-product matrix
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

    # Rodriguez rotation formula
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s**2))

    return rotation_matrix
