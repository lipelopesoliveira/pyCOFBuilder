import numpy as np
from numpy.testing import assert_allclose

import sys
sys.path.insert(1, '/home/felipe/PRs/pyCOFBuilder/src')

from pycofbuilder.tools import (elements_dict,
                                unit_vector,
                                angle,
                                rotation_matrix_from_vectors,
                                rmsd)


def test_elements_dict() -> None:
    assert type(elements_dict()) == dict


def test_unit_vector_type() -> None:
    """Test that the unit_vector function returns a numpy array."""
    assert type(unit_vector([1, 0, 0])) == np.ndarray


def test_unit_vector_normalization():
    """Test that the unit_vector function returns a vector of length 1."""
    vec = np.array([3, 4, 0])
    uv = unit_vector(vec)
    # The norm of the returned vector should be 1.
    norm = np.linalg.norm(uv)
    assert_allclose(norm, 1.0, atol=1e-5)
    # Also check that the unit vector points in the same direction as the original.
    expected = vec / np.linalg.norm(vec)
    assert_allclose(uv, expected, atol=1e-5)


def test_angle_type():
    """Test that the angle function returns a float."""
    assert type(angle([1, 0, 0], [1, 0, 0])) == np.float64


def test_angle_same_vector():
    """Test that the angle between a vector and itself is zero."""
    v = np.array([1, 2, 3])

    # Angle in degrees should be 0
    assert np.isclose(angle(v, v, unit='degree'), 0.0, atol=1e-5)

    # Angle in radians should be 0
    assert np.isclose(angle(v, v, unit='radians'), 0.0, atol=1e-5)

    # The cosine of the angle between identical vectors should be 1
    assert np.isclose(angle(v, v, unit='cos'), 1.0, atol=1e-5)


def test_angle_orthogonal_vectors():
    """Test that the angle between two orthogonal vectors is 90 degrees (or π/2 radians)."""
    v1 = np.array([1, 0, 0])
    v2 = np.array([0, 1, 0])

    # 90 degrees
    assert np.isclose(angle(v1, v2, unit='degree'), 90.0, atol=1e-5)

    # π/2 radians
    assert np.isclose(angle(v1, v2, unit='radians'), np.pi/2, atol=1e-5)


def test_rotation_matrix_identity():
    """Test that the rotation matrix between a vector and itself is the identity matrix."""
    v = np.array([1, 2, 3])
    R = rotation_matrix_from_vectors(v, v)
    expected_identity = np.eye(3)
    assert_allclose(R, expected_identity, atol=1e-5)


def test_rotation_matrix_alignment():
    """Test that the rotation matrix correctly aligns one vector to another."""
    v1 = np.array([1, 0, 0])
    v2 = np.array([0, 1, 0])
    R = rotation_matrix_from_vectors(v1, v2)
    # Rotate v1 using the rotation matrix.
    rotated_v1 = np.dot(R, v1)
    # The expected result should be in the direction of v2 (normalized).
    expected_v2 = v2 / np.linalg.norm(v2)
    # Check that the rotated vector is close to expected.
    assert_allclose(rotated_v1, expected_v2, atol=1e-5)


def test_rmsd_identical():
    """Test that the RMSD between two identical sets of vectors is zero."""
    V = np.array([[0, 0, 0], [1, 1, 1]])
    W = np.array([[0, 0, 0], [1, 1, 1]])
    assert np.isclose(rmsd(V, W), 0.0, atol=1e-5)


def test_rmsd_nonzero():
    """Test the RMSD function with two different sets of vectors."""
    V = np.array([[0, 0, 0], [1, 1, 1]])
    W = np.array([[1, 1, 1], [2, 2, 2]])

    # Compute the expected RMSD manually.
    diff = V - W
    expected = np.sqrt((diff * diff).sum() / 2)
    assert np.isclose(rmsd(V, W), expected, atol=1e-5)
