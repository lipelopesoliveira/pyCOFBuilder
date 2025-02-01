import os
import numpy as np
from numpy.testing import assert_allclose
# import pytest

from pycofbuilder.io_tools import (read_xyz)

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
