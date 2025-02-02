# -*- coding: utf-8 -*-
# Created by Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

"""
This module contains the geometry tools used by pyCOFBuilder.
"""

import numpy as np


def distance(p1: list, p2: list, ndecimals=1):
    """Calculate the distance between two points in a N-dimentional space rounded to a given number of decimal places.

    The rounding is done to avoid that small distortions generate errors when comparing distances.

    Parameters
    ----------

    p1 : list
        The first point in the N-dimentional space.
    p2: list
        The second point in the N-dimentional space.
    ndecimals : int, optional
        The number of decimal places to round the distance. Default is 1.

    Returns
    -------
    float
        The distance between the two points.
    """
    return np.round(np.linalg.norm(np.array(p1) - np.array(p2)), ndecimals)


def angle_between(v1, v2):
    """Calculates the angle in radians between two vectors v1 and v2."""
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    return np.arccos(dot_product / (norm_v1 * norm_v2))


def is_equilateral(p1, p2, p3, tol=0.1):
    """Check if the three points form an equilateral triangle."""
    # Convert points to NumPy arrays
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)

    # Calculate the distances using NumPy
    d1 = np.linalg.norm(p1 - p2)
    d2 = np.linalg.norm(p2 - p3)
    d3 = np.linalg.norm(p3 - p1)

    # Check if all distances are equal (with a tolerance for floating point comparison)
    return np.isclose(d1, d2, rtol=tol) and np.isclose(d2, d3, rtol=tol)


def tetrahedron_volume(p1, p2, p3, p4):
    """Calculate the volume of the tetrahedron formed by four points."""
    # Create a matrix with the points and an additional row of ones
    matrix = np.array([
        [p1[0], p1[1], p1[2], 1],
        [p2[0], p2[1], p2[2], 1],
        [p3[0], p3[1], p3[2], 1],
        [p4[0], p4[1], p4[2], 1]
    ])

    # Calculate the determinant of the matrix
    det = np.linalg.det(matrix)

    # Volume is absolute value of the determinant divided by 6
    volume = abs(det) / 6
    return volume


def is_perfect_tetrahedron(p1, p2, p3, p4, tolerance=0.05):
    """Check if the tetrahedron formed by four points is perfect."""
    # Calculate all edge lengths
    edges = [
        distance(p1, p2),
        distance(p1, p3),
        distance(p1, p4),
        distance(p2, p3),
        distance(p2, p4),
        distance(p3, p4),
    ]

    # Find unique edge lengths
    unique_edges = set(edges)

    # A perfect tetrahedron should have exactly one unique edge length
    if len(unique_edges) != 1:
        return False

    edge_length = next(iter(unique_edges))  # Get the unique edge length

    # Check if all edges are within the tolerance of this length
    return all(np.isclose(edge, edge_length, rtol=tolerance) for edge in edges)


def is_distorted_tetrahedron(p1, p2, p3, p4, tolerance=0.05):
    """Check if four points form a tetrahedron with a given relative tolerance."""
    volume = tetrahedron_volume(p1, p2, p3, p4)
    return volume > tolerance  # Points must not be coplanar


def is_square(p1, p2, p3, p4, tolerance=0.05):
    """Check if four points form a square with a given relative tolerance."""
    dists = [
        distance(p1, p2, 0),
        distance(p1, p3, 0),
        distance(p1, p4, 0),
        distance(p2, p3, 0),
        distance(p2, p4, 0),
        distance(p3, p4, 0),
    ]

    unique_dists = set(dists)

    # A square has 2 unique distances: side length (4 times) and diagonal length (2 times)
    if len(unique_dists) != 2:
        return False

    side_length = min(unique_dists)
    diagonal_length = max(unique_dists)

    # Check if the counts are correct with tolerance
    return (dists.count(side_length) == 4 and
            dists.count(diagonal_length) == 2 and
            np.isclose(diagonal_length, 2 * side_length, rtol=tolerance))


def is_rectangle(p1, p2, p3, p4, tolerance=0.1):
    """Check if four points form a rectangle with a given relative tolerance."""
    dists = [
        distance(p1, p2),
        distance(p1, p3),
        distance(p1, p4),
        distance(p2, p3),
        distance(p2, p4),
        distance(p3, p4),
    ]

    unique_dists = set(dists)

    # A rectangle has 3 unique distances: length (2 times) and width (2 times), and diagonals (2 times)
    if len(unique_dists) != 3:
        return False

    side_length = sorted(list(unique_dists))[:2]
    diagonal = sorted(list(unique_dists))[-1]

    # Check if the counts are correct with tolerance
    return (dists.count(side_length[0]) == 2 and dists.count(side_length[1]) and dists.count(diagonal) == 2)


def is_rhombus(p1, p2, p3, p4, tolerance=0.1):
    """Check if four points form a rhombus with a given relative tolerance."""
    dists = [
        distance(p1, p2),
        distance(p1, p3),
        distance(p1, p4),
        distance(p2, p3),
        distance(p2, p4),
        distance(p3, p4),
    ]

    unique_dists = set(dists)

    # A rhombus at least 2 unique distances: (1 or 2) different sides and (1 or 2) diagonal lengths
    if len(unique_dists) < 2 or len(unique_dists) > 4:
        return False

    volume = tetrahedron_volume(p1, p2, p3, p4)

    return 2 <= len(unique_dists) <= 4 and volume < tolerance


def is_pyramid(p1, p2, p3, p4, p5, tolerance=0.05):
    """Check if five points form a pyramid with a square base."""
    # Check if the base is a square
    if not is_square(p1, p2, p3, p4, tolerance):
        return False

    # Distances from the apex (p5) to the base
    apex_edges = [
        distance(p1, p5),
        distance(p2, p5),
        distance(p3, p5),
        distance(p4, p5)
    ]

    # Check if all edges from the apex are the same length
    if len(set(apex_edges)) != 1:
        return False

    return True


def is_hexagon(p1, p2, p3, p4, p5, p6, tolerance=0.1):
    """Check if six points form a regular hexagon."""
    points = [p1, p2, p3, p4, p5, p6]

    # Check if all distances between consecutive points are equal
    distances = [distance(points[i], points[(i+1) % 6]) for i in range(6)]

    unique_distances = set(distances)

    if len(unique_distances) != 1:
        return False  # If the distances are not all equal, it's not a regular hexagon

    # Check if the angle between adjacent vectors is 120 degrees (pi/3 radians)
    for i in range(6):
        # Create vectors between consecutive points
        v1 = np.array([points[(i+1) % 6][0] - points[i][0], points[(i+1) % 6][1] - points[i][1]])
        v2 = np.array([points[(i+2) % 6][0] - points[(i+1) % 6][0], points[(i+2) % 6][1] - points[(i+1) % 6][1]])

        # Calculate angle between vectors
        angle = angle_between(v1, v2)

        # Check if the angle is approximately 120 degrees
        if not np.isclose(angle, np.pi / 3, atol=tolerance):
            return False

    return True


def is_octahedron(p1, p2, p3, p4, p5, p6, tolerance=0.1):
    """Check if six points form an octahedron."""
    # Distances between the points of an octahedron
    edges = [
        distance(p1, p2),
        distance(p1, p3),
        distance(p1, p4),
        distance(p1, p5),
        distance(p1, p6),
        distance(p2, p3),
        distance(p2, p4),
        distance(p2, p5),
        distance(p2, p6),
        distance(p3, p4),
        distance(p3, p5),
        distance(p3, p6),
        distance(p4, p5),
        distance(p4, p6),
        distance(p5, p6)
    ]

    # For an octahedron, we expect all edges to be of equal length
    if len(set(edges)) != 1:
        return False

    return True


def is_cube(p1, p2, p3, p4, p5, p6, p7, p8, tolerance=0.05):
    """Check if eight points form a cube."""
    edges = [
        distance(p1, p2), distance(p1, p3), distance(p1, p4), distance(p1, p5),
        distance(p2, p3), distance(p2, p6), distance(p3, p4), distance(p3, p7),
        distance(p4, p5), distance(p5, p6), distance(p6, p7), distance(p7, p8),
        distance(p8, p5)
    ]

    # Check if there are exactly 3 unique distances for the sides, diagonal, and face diagonal
    unique_edges = set(edges)

    if len(unique_edges) != 3:
        return False

    edge_lengths = sorted(list(unique_edges))

    # A cube has 12 equal edges, and face diagonals equal length
    return (edges.count(edge_lengths[0]) == 12)


def is_parallelepiped(p1, p2, p3, p4, p5, p6, p7, p8, tolerance=0.1):
    """Check if eight points form a parallelepiped."""
    # Calculate all pairwise distances between the points
    points = [p1, p2, p3, p4, p5, p6, p7, p8]
    distances = [distance(p1, p2) for i, p1 in enumerate(points) for j, p2 in enumerate(points) if i < j]

    # Get unique distances
    unique_distances = sorted(set(distances))

    # A parallelepiped has 3 unique edge lengths (for width, height, and depth)
    if len(unique_distances) != 3:
        return False

    # Check if the distances match the expected edge lengths (with tolerance)
    edge_lengths = sorted(list(unique_distances))

    # In a parallelepiped, we expect 4 edges of each length
    for edge_length in edge_lengths:
        if distances.count(edge_length) != 4:
            return False

    return True


def is_dodecahedron(*points, tolerance=0.1):
    """Check if twelve points form a dodecahedron."""
    if len(points) != 12:
        return False

    # Calculate all pairwise distances
    distances = [distance(p1, p2) for i, p1 in enumerate(points) for j, p2 in enumerate(points) if i < j]

    unique_distances = set(distances)

    # For a dodecahedron, we expect several unique distances (vertices and edge distances)
    # In general, dodecahedron has 30 edges and 20 faces
    return len(unique_distances) > 5
