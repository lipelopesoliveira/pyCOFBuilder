# -*- coding: utf-8 -*-

import os
import sys

# This code is written for Python 3
if sys.version_info[0] != 3:
    raise Exception("Sorry but pyCOFBuilder requires Python 3.")
    sys.exit(1)

# Import tools class
import pycofbuilder.tools as Tools

# Import BuildingBlocks class
from pycofbuilder.building_block import Building_Block

# Import BuildingBlocks class
from pycofbuilder.framework import Framework

# Import Core class
#from pycofbuilder.core import *

__all__ = ['Tools',
           'Building_Block',
           'Framework'
           ]

_ROOT = os.path.abspath(os.path.dirname(__file__))

__author__  = "Felipe Lopes de Oliveira"
__license__ = "MIT"
__version__ = '0.0.4'
__email__   = "felipe.lopes@nano.ufrj.br"
__status__  = "Development"
