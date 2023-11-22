# -*- coding: utf-8 -*-
# Created by Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

import os
import sys

# This code is written for Python 3
if sys.version_info[0] != 3:
    raise Exception("Sorry but pyCOFBuilder requires Python 3.")
    sys.exit(1)

# Import BuildingBlocks class
from pycofbuilder.building_block import BuildingBlock

# Import Framework class
from pycofbuilder.framework import Framework

# Import Tools
import pycofbuilder.tools as Tools
import pycofbuilder.io_tools as IO_Tools

__all__ = [
    'BuildingBlock',
    'Framework',
    'Tools',
    'IO_Tools',
    'ChemJSON'
    ]

_ROOT = os.path.abspath(os.path.dirname(__file__))

__author__ = "Felipe Lopes de Oliveira"
__license__ = "MIT"
__version__ = '0.0.5'
__email__ = "felipe.lopes@nano.ufrj.br"
__status__ = "Development"
