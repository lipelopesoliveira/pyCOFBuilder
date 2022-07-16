__author__ = "Felipe Lopes de Oliveira"
__license__ = "MIT"

try:
    from ._version import version
    __version__ = version
except ImportError:
    __version__ = '0.0.2'

import sys

# This code is written for Python 3.
if sys.version_info[0] != 3:
    raise Exception("pyCOFBuilder requires Python 3.")
    sys.exit(1)

# Import tools class
import pycofbuilder.tools as Tools

# Import BuildingBlocks class
from pycofbuilder.building_block import Building_Block

# Import BuildingBlocks class
from pycofbuilder.reticulum import Reticulum

# Import Core class
from pycofbuilder.core import *


__all__ = ['Tools',
           'Building_Block',
           'Reticulum'
           ]

import os
_ROOT = os.path.abspath(os.path.dirname(__file__))
