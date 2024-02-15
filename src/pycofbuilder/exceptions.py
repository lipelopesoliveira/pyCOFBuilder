# -*- coding: utf-8 -*-
# Created by Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

"""
This file implements custom exceptions for the pycofbuilder package.
"""


class BondLenghError(Exception):
    """Exception raised when the bond length is shorter than 0.8 angstrom."""
    def __init__(self, atom_1=None, atom_2=None, dist=0, threshold=0.8):
        self.message = 'WARNING: Atoms {} and {} are closer than {} A, {}'.format(atom_1,
                                                                                  atom_2,
                                                                                  dist,
                                                                                  threshold)

    def __str__(self):
        return str(self.message)


class BBConnectivityError(Exception):
    """Exception raised when the building block connectivity is not valid."""
    def __init__(self, connectivity=None, found_connectivity=None):
        self.message = 'ERROR: The building block connectivity should be {} buy is {}'.format(connectivity,
                                                                                              found_connectivity)

    def __str__(self):
        return str(self.message)
