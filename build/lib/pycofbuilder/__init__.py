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

# Import ChemJSON
import pycofbuilder.cjson as ChemJSON

# Import Exceptions
import pycofbuilder.exceptions as Exceptions

# Import Logger
import pycofbuilder.logger as Logger

__all__ = [
    'BuildingBlock',
    'Framework',
    'Tools',
    'IO_Tools',
    'ChemJSON',
    'Exceptions',
    'Logger'
    ]

_ROOT = os.path.abspath(os.path.dirname(__file__))

__author__ = "Felipe Lopes de Oliveira"
__license__ = "MIT"
__version__ = '0.0.8.1'
__email__ = "felipe.lopes@nano.ufrj.br"
__status__ = "Development"
