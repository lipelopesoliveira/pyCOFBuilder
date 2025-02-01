# -*- coding: utf-8 -*-

from pycofbuilder.tools import (elements_dict,
                                unit_vector,
                                angle)

import unittest
import numpy as np


class BasicToolsTest(unittest.TestCase):
    """Basic test of Tools module."""

    def test_elements_dict(self) -> None:
        assert type(elements_dict()) == dict

    def test_unit_vector(self) -> None:
        assert type(unit_vector([1, 0, 0])) == np.ndarray
        assert np.array_equal(unit_vector([1, 0, 0]), np.array([1, 0, 0]))
        assert np.array_equal(unit_vector([0, 1, 0]), np.array([0, 1, 0]))
        assert np.array_equal(unit_vector([0, 0, 1]), np.array([0, 0, 1]))
        assert np.array_equal(unit_vector([1, 1, 1]), np.array([1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)]))

    def test_angle(self):
        assert type(angle([1, 0, 0], [1, 0, 0])) == np.float64
        assert angle([1, 0, 0], [0, 1, 0]) == 90.0
        assert angle([1, 0, 0], [0, 0, 1]) == 90.0
        assert angle([0, 1, 0], [0, 0, 1]) == 90.0
        assert angle([1, 0, 0], [1, 0, 0]) == 0.0
        assert angle([0, 1, 0], [0, 1, 0]) == 0.0

        # Test for 45 degrees
        assert angle([1, 0, 0], [1, 1, 0]) - 45.0 < 1e-6
        assert angle([1, 0, 0], [1, 0, 1]) - 45.0 < 1e-6
        assert angle([0, 1, 0], [1, 1, 0]) - 45.0 < 1e-6
        assert angle([0, 1, 0], [0, 1, 1]) - 45.0 < 1e-6
        assert angle([0, 0, 1], [1, 0, 1]) - 45.0 < 1e-6
        assert angle([0, 0, 1], [0, 1, 1]) - 45.0 < 1e-6

    def test_calculate_sides(self):
        pass

    def test_cell_to_cellpar(self):
        pass

class BasicIOToolsTest(unittest.TestCase):
    """Basic test of IO_Tools module."""

    def test_read_cif(self):
        pass


if __name__ == '__main__':
    unittest.main()
