import os
import numpy as np
from numpy.testing import assert_allclose
import pytest

from pycofbuilder.io_tools import (read_xyz, read_pdb, read_gjf, read_cif)

_ROOTDIR = os.path.abspath(os.path.dirname(__file__))


def test_read_xyz():
    """Test the read_xyz function."""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    atomTypes, cartPos, cellMatrix = read_xyz(file_path, 'CH4.xyz')

    # Check the number of atoms
    assert len(atomTypes) == 5

    # Check the number of coordinates
    assert cartPos.shape == (5, 3)

    # Check the first atom
    assert atomTypes[0] == 'C'

    # Check the first coordinate
    assert_allclose(cartPos[0], [0.000000, 0.000000, 0.000000])

    # Check the last atom
    assert atomTypes[-1] == 'H'

    # Check the last coordinate
    assert_allclose(cartPos[-1], [-0.62473386, 0.62473386, -0.62473386])

    # Check cellMatrix shape
    assert cellMatrix.shape == (3, 3)

    # Check cellMatrix values
    assert_allclose(cellMatrix, np.zeros((3, 3)), atol=1e-5)


def test_read_xyz_ase():
    """Test the read_xyz function with ase."""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    atomTypes, cartPos, cellMatrix = read_xyz(file_path, 'CH4.xyz', extxyz=True)

    # Check the number of atoms
    assert len(atomTypes) == 5

    # Check the number of coordinates
    assert cartPos.shape == (5, 3)

    # Check the first atom
    assert atomTypes[0] == 'C'

    # Check the first coordinate
    assert_allclose(cartPos[0], [0.000000, 0.000000, 0.000000])

    # Check the last atom
    assert atomTypes[-1] == 'H'

    # Check the last coordinate
    assert_allclose(cartPos[-1], [-0.62473386, 0.62473386, -0.62473386])

    # Check cellMatrix shape
    assert cellMatrix.shape == (3, 3)

    # Check cellMatrix values
    assert_allclose(cellMatrix, np.zeros((3, 3)), atol=1e-5)


def test_read_extxyz():
    """Test the read_xyz function with ase."""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    atomTypes, cartPos, cellMatrix = read_xyz(file_path, 'CH4_extxyz.xyz', extxyz=True)

    # Check the number of atoms
    assert len(atomTypes) == 5

    # Check the number of coordinates
    assert cartPos.shape == (5, 3)

    # Check the first atom
    assert atomTypes[0] == 'C'

    # Check the first coordinate
    assert_allclose(cartPos[0], [0.000000, 0.000000, 0.000000])

    # Check the last atom
    assert atomTypes[-1] == 'H'

    # Check the last coordinate
    assert_allclose(cartPos[-1], [-0.62473386, 0.62473386, -0.62473386])

    # Check cellMatrix shape
    assert cellMatrix.shape == (3, 3)

    # Check cellMatrix values
    assert_allclose(cellMatrix, np.zeros((3, 3)), atol=1e-5)


def test_read_extxyz_false():
    """Test the read_xyz function with ase."""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    atomTypes, cartPos, cellMatrix = read_xyz(file_path, 'CH4_extxyz.xyz', extxyz=False)

    # Check the number of atoms
    assert len(atomTypes) == 5

    # Check the number of coordinates
    assert cartPos.shape == (5, 3)

    # Check the first atom
    assert atomTypes[0] == 'C'

    # Check the first coordinate
    assert_allclose(cartPos[0], [0.000000, 0.000000, 0.000000])

    # Check the last atom
    assert atomTypes[-1] == 'H'

    # Check the last coordinate
    assert_allclose(cartPos[-1], [-0.62473386, 0.62473386, -0.62473386])

    # Check cellMatrix shape
    assert cellMatrix.shape == (3, 3)

    # Check cellMatrix values
    assert_allclose(cellMatrix, np.zeros((3, 3)), atol=1e-5)


def test_read_xyz_exception():
    """Test the FileNotFoundError of read_xyz function"""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    with pytest.raises(FileNotFoundError):
        read_xyz(file_path, 'non-existent.xyz')


def test_read_pdb_default():
    """Test the read_pdb function."""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    atomTypes, cartPos, cellMatrix = read_pdb(file_path, 'CH4.pdb')

    # Check the number of atoms
    assert len(atomTypes) == 5

    # Check the number of coordinates
    assert cartPos.shape == (5, 3)

    # Check the first atom
    assert atomTypes[0] == 'C'

    # Check the first coordinate
    assert_allclose(cartPos[0], [0.000000, 0.000000, 0.000000], atol=1e-3)

    # Check the last atom
    assert atomTypes[-1] == 'H'

    # Check the last coordinate
    assert_allclose(cartPos[-1], [-0.62473386, 0.62473386, -0.62473386], atol=1e-3)

    # Check cellMatrix shape
    assert cellMatrix.shape == (3, 3)

    # Check cellMatrix values
    assert_allclose(cellMatrix, np.zeros((3, 3)), atol=1e-5)


def test_read_pdb_vesta():
    """Test the read_pdb function."""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    atomTypes, cartPos, cellMatrix = read_pdb(file_path, 'CH4_VESTA.pdb')

    # Check the number of atoms
    assert len(atomTypes) == 5

    # Check the number of coordinates
    assert cartPos.shape == (5, 3)

    # Check the first atom
    assert atomTypes[0] == 'C'

    # Check the first coordinate
    assert_allclose(cartPos[0], [0.000000, 0.000000, 0.000000], atol=1e-3)

    # Check the last atom
    assert atomTypes[-1] == 'H'

    # Check the last coordinate
    assert_allclose(cartPos[-1], [-0.62473386, 0.62473386, -0.62473386], atol=1e-3)

    # Check cellMatrix shape
    assert cellMatrix.shape == (3, 3)

    # Check cellMatrix values
    assert_allclose(cellMatrix, np.zeros((3, 3)), atol=1e-5)


def test_read_pdb_exception():
    """Test the FileNotFoundError of read_pdb function"""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    with pytest.raises(FileNotFoundError):
        read_pdb(file_path, 'non-existent.pdb')


def test_read_gjf_mol():
    """Test the read_gjf function for molecules."""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    atomTypes, cartPos, cellMatrix = read_gjf(file_path, 'CH4.gjf')

    # Check the number of atoms
    assert len(atomTypes) == 5

    # Check the number of coordinates
    assert cartPos.shape == (5, 3)

    # Check the first atom
    assert atomTypes[0] == 'C'

    # Check the first coordinate
    assert_allclose(cartPos[0], [0.000000, 0.000000, 0.000000], atol=1e-3)

    # Check the last atom
    assert atomTypes[-1] == 'H'

    # Check the last coordinate
    assert_allclose(cartPos[-1], [-0.62473386, 0.62473386, -0.62473386], atol=1e-3)

    # Check cellMatrix shape
    assert cellMatrix.shape == (3, 3)

    # Check cellMatrix values
    assert_allclose(cellMatrix, np.zeros((3, 3)), atol=1e-5)


def test_read_gjf_crystal():
    """Test the read_gjf function for crystals."""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    atomTypes, cartPos, cellMatrix = read_gjf(file_path, 'RIO_13.gjf')

    # Check the number of atoms
    assert len(atomTypes) == 42

    # Check the number of coordinates
    assert cartPos.shape == (42, 3)

    # Check the first atom
    assert atomTypes[0] == 'H'

    # Check the first coordinate
    assert_allclose(cartPos[0], [1.76105047, 5.64392724, 1.80730000], atol=1e-3)

    # Check the last atom
    assert atomTypes[-1] == 'O'

    # Check the last coordinate
    assert_allclose(cartPos[-1], [4.90080521, 3.55499367, 1.80730000], atol=1e-3)

    # Check cellMatrix shape
    assert cellMatrix.shape == (3, 3)

    RIO_13_cell = np.array(
        [[15.05300000, 0.00000000, 0.00000000],
         [-7.52650000, 13.03628040, 0.00000000],
         [0.00000000, 0.00000000, 3.61460000]]
    )

    # Check cellMatrix values
    assert_allclose(cellMatrix, RIO_13_cell, atol=1e-5)


def test_read_gjf_exception():
    """Test the FileNotFoundError of read_gjf function"""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    with pytest.raises(FileNotFoundError):
        read_gjf(file_path, 'non-existent.gjf')


def test_read_cif_symmetry_default():
    """Test the read_cif function for files with symmetry."""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    atomTypes, cartPos, cellMatrix, partialCharges = read_cif(file_path, 'RIO_13_symmetry.cif', usePymatgen=True)

    # Check the number of atoms
    assert len(atomTypes) == 42

    # Check the number of coordinates
    assert cartPos.shape == (42, 3)

    # Check cellMatrix shape
    assert cellMatrix.shape == (3, 3)

    RIO_13_cell = np.array(
        [[15.05300000, 0.00000000, 0.00000000],
         [-7.52650000, 13.03628040, 0.00000000],
         [0.00000000, 0.00000000, 3.61460000]]
    )

    # Check cellMatrix values
    assert_allclose(cellMatrix, RIO_13_cell, atol=1e-5)


def test_read_cif_symmetry_ase():
    """Test the read_cif function for files with symmetry."""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    atomTypes, cartPos, cellMatrix, partialCharges = read_cif(file_path, 'RIO_13_symmetry.cif', useASE=True)

    # Check the number of atoms
    assert len(atomTypes) == 42

    # Check the number of coordinates
    assert cartPos.shape == (42, 3)

    # Check the first atom
    assert atomTypes[0] == 'H'

    # Check the first coordinate
    assert_allclose(cartPos[0], [1.76105047, 5.64392724, 1.80730000], atol=1e-3)

    # Check the last atom
    assert atomTypes[-1] == 'O'

    # Check the last coordinate
    assert_allclose(cartPos[-1], [4.90080521, 3.55499367, 1.80730000], atol=1e-3)

    # Check cellMatrix shape
    assert cellMatrix.shape == (3, 3)

    RIO_13_cell = np.array(
        [[15.05300000, 0.00000000, 0.00000000],
         [-7.52650000, 13.03628040, 0.00000000],
         [0.00000000, 0.00000000, 3.61460000]]
    )

    # Check cellMatrix values
    assert_allclose(cellMatrix, RIO_13_cell, atol=1e-5)


def test_read_cif_symmetry_pymatgen():
    """Test the read_cif function for files with symmetry."""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    atomTypes, cartPos, cellMatrix, partialCharges = read_cif(file_path, 'RIO_13_symmetry.cif', usePymatgen=True)

    # Check the number of atoms
    assert len(atomTypes) == 42

    # Check the number of coordinates
    assert cartPos.shape == (42, 3)

    # Check cellMatrix shape
    assert cellMatrix.shape == (3, 3)

    RIO_13_cell = np.array(
        [[15.05300000, 0.00000000, 0.00000000],
         [-7.52650000, 13.03628040, 0.00000000],
         [0.00000000, 0.00000000, 3.61460000]]
    )

    # Check cellMatrix values
    assert_allclose(cellMatrix, RIO_13_cell, atol=1e-5)


def test_read_cif_charges_default():
    """Test the read_cif function for files with symmetry."""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    atomTypes, cartPos, cellMatrix, partialCharges = read_cif(file_path, 'RIO_13_P1_charge.cif')

    # Check the number of atoms
    assert len(atomTypes) == 42

    # Check the number of coordinates
    assert cartPos.shape == (42, 3)

    # Check the first atom
    assert atomTypes[0] == 'H'

    # Check the first coordinate
    # assert_allclose(cartPos[0], [1.75914376, 5.64392724, 1.80730000], atol=1e-3)

    # Check the first partial charge
    assert_allclose(partialCharges[0], 0.03650, atol=1e-5)

    # Check the last atom
    assert atomTypes[-1] == 'O'

    # Check the last coordinate
    # assert_allclose(cartPos[-1], [9.52388257, 2.46672498, 1.80730000], atol=1e-3)

    # Check the first partial charge
    assert_allclose(partialCharges[-1], -0.38067, atol=1e-5)

    # Check cellMatrix shape
    assert cellMatrix.shape == (3, 3)

    RIO_13_cell = np.array(
        [[15.05300000, 0.00000000, 0.00000000],
         [-7.52650000, 13.03628040, 0.00000000],
         [0.00000000, 0.00000000, 3.61460000]]
    )

    # Check cellMatrix values
    assert_allclose(cellMatrix, RIO_13_cell, atol=1e-5)


def test_read_cif_charges_pymatgen():
    """Test the read_cif function for files with symmetry."""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    atomTypes, cartPos, cellMatrix, partialCharges = read_cif(file_path, 'RIO_13_P1_charge.cif', usePymatgen=True)

    # Check the number of atoms
    assert len(atomTypes) == 42

    # Check the number of coordinates
    assert cartPos.shape == (42, 3)

    # Check cellMatrix shape
    assert cellMatrix.shape == (3, 3)

    RIO_13_cell = np.array(
        [[15.05300000, 0.00000000, 0.00000000],
         [-7.52650000, 13.03628040, 0.00000000],
         [0.00000000, 0.00000000, 3.61460000]]
    )

    # Check cellMatrix values
    assert_allclose(cellMatrix, RIO_13_cell, atol=1e-5)


def test_read_cif_charges_ase():
    """Test the read_cif function for files with symmetry."""

    # Read the file
    file_path = os.path.join(_ROOTDIR, 'data')

    atomTypes, cartPos, cellMatrix, partialCharges = read_cif(file_path, 'RIO_13_P1_charge.cif', useASE=True)

    # Check the number of atoms
    assert len(atomTypes) == 42

    # Check the number of coordinates
    assert cartPos.shape == (42, 3)

    # Check the first atom
    assert atomTypes[0] == 'H'

    # Check the first coordinate
    assert_allclose(cartPos[0], [1.75914376, 5.64392724, 1.80730000], atol=1e-3)

    # Check the first partial charge
    assert_allclose(partialCharges[0], 0.0, atol=1e-5)

    # Check the last atom
    assert atomTypes[-1] == 'O'

    # Check the last coordinate
    assert_allclose(cartPos[-1], [9.52388257, 2.46672498, 1.80730000], atol=1e-3)

    # Check the first partial charge
    assert_allclose(partialCharges[-1], 0.0, atol=1e-5)

    # Check cellMatrix shape
    assert cellMatrix.shape == (3, 3)

    RIO_13_cell = np.array(
        [[15.05300000, 0.00000000, 0.00000000],
         [-7.52650000, 13.03628040, 0.00000000],
         [0.00000000, 0.00000000, 3.61460000]]
    )

    # Check cellMatrix values
    assert_allclose(cellMatrix, RIO_13_cell, atol=1e-5)
