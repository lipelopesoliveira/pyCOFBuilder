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


class ConnectionGroupError(Exception):
    """Exception raised when the connection group is not valid."""
    def __init__(self, conn_1=None, conn_2=None):
        self.message = 'ERROR: The connection group {} not compatible with {}'.format(conn_1,
                                                                                      conn_2)

    def __str__(self):
        return str(self.message)


class MissingXError(Exception):
    """Exception raised when the building block has missing X atoms."""
    def __init__(self):
        self.message = 'No X points found in the structure!'

    def __str__(self):
        return str(self.message)
