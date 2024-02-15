# -*- coding: utf-8 -*-
# Created by Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

"""
The dictionary containing the definitions nets for a Framework buiding
"""

import numpy as np

TOPOLOGY_DICT = {
    'HCB': {
        'a': 2*np.cos(np.radians(30)),
        'b': 2*np.cos(np.radians(30)),
        'c': 3.6,
        'alpha': 90,
        'beta': 90,
        'gamma': 120,
        'vertice_connectivity': 3,
        'edge_connectivity': 0,
        'vertices': [
            {'position': [0, 0, 0], 'angle': 0},
            {'position': [0, np.sqrt(3)/3, 0], 'angle': 180}
            ],
        'edges': []
        },
    'HCB_A': {
        'a': 2*np.cos(np.radians(30))*2,
        'b': 2*np.cos(np.radians(30))*2,
        'c': 3.6,
        'alpha': 90,
        'beta': 90,
        'gamma': 120,
        'vertice_connectivity': 3,
        'edge_connectivity': 3,
        'vertices': [
            {'position': [0, 0, 0], 'angle': 0},
            {'position': [0, np.sqrt(3)/3, 0], 'angle': 180}
            ],
        'edges': [
            {'position': [0, np.sqrt(3)/6, 0], 'angle': 0},
            {'position': [-1/4, 5*np.sqrt(3)/12, 0], 'angle': 120},
            {'position': [1/4, 5*np.sqrt(3)/12, 0], 'angle': 240}
            ]
        },
    'SQL': {
        'a': 2/np.sqrt(2),
        'b': 2/np.sqrt(2),
        'c': 3.6,
        'alpha': 90,
        'beta': 90,
        'gamma': 90,
        'vertice_connectivity': 4,
        'edge_connectivity': 0,
        'vertices': [
            {'position': [0, 0, 0], 'angle': 0},
            {'position': [1/2, 1/2, 0], 'angle': 0}
            ],
        'edges': []
        },
    'SQL_A': {
        'a': 4/np.sqrt(2),
        'b': 4/np.sqrt(2),
        'c': 3.6,
        'alpha': 90,
        'beta': 90,
        'gamma': 90,
        'vertice_connectivity': 4,
        'edge_connectivity': 2,
        'vertices': [
            {'position': [0, 0, 0], 'angle': 0},
            {'position': [1/2, 1/2, 0], 'angle': 0}
            ],
        'edges': [
            {'position': [1/4, 1/4, 0], 'angle': 45},
            {'position': [3/4, 1/4, 0], 'angle': 135},
            {'position': [3/4, 3/4, 0], 'angle': 225},
            {'position': [1/4, 3/4, 0], 'angle': 315},
            ]
        },
    'KGD': {
        'a': 2*np.cos(np.radians(30)),
        'b': 2*np.cos(np.radians(30)),
        'c': 3.6,
        'alpha': 90,
        'beta': 90,
        'gamma': 120,
        'vertice_connectivity': 6,
        'edge_connectivity': 3,
        'vertices': [
            {'position': [0, 0, 0], 'angle': 0},
            ],
        'edges': [
            {'position': [0, np.sqrt(3)/3, 0], 'angle': -180},
            {'position': [0.5, np.sqrt(3)/6, 0], 'angle': 0}
            ]
        },
    'HXL_A': {
        'a': 2,
        'b': 2,
        'c': 3.6,
        'alpha': 90,
        'beta': 90,
        'gamma': 120,
        'vertice_connectivity': 6,
        'edge_connectivity': 0,
        'vertices': [
            {'position': [0, 0, 0], 'angle': 30},
            ],
        'edges': [
            {'position': [1/4, np.sqrt(3)/4, 0], 'angle': 30},
            {'position': [0.5, 0, 0], 'angle': 90},
            {'position': [-1/4, np.sqrt(3)/4, 0], 'angle': -30}
            ]
        },
    'KGM': {
        'a': 4,
        'b': 4,
        'c': 3.6,
        'alpha': 90,
        'beta': 90,
        'gamma': 120,
        'vertice_connectivity': 4,
        'edge_connectivity': 0,
        'vertices': [
            {'position': [1/4, np.sqrt(3)/4, 0], 'angle': 30},
            {'position': [1/2, 0, 0], 'angle': -90},
            {'position': [-1/4, np.sqrt(3)/4, 0], 'angle': -30}
            ],
        'edges': []
        },
    'KGM_A': {
        'a': 4,
        'b': 4,
        'c': 3.6,
        'alpha': 90,
        'beta': 90,
        'gamma': 120,
        'vertice_connectivity': 4,
        'edge_connectivity': 2,
        'vertices': [
            {'position': [1/4, np.sqrt(3)/4, 0], 'angle': 30},
            {'position': [1/2, 0, 0], 'angle': -90},
            {'position': [-1/4, np.sqrt(3)/4, 0], 'angle': -30}
            ],
        'edges': [
            {'position': [3/8, np.sqrt(3)/8, 0], 'angle': -30},
            {'position': [1/8, 3*np.sqrt(3)/8, 0], 'angle': -30},
            {'position': [5/8, np.sqrt(3)/8, 0], 'angle': 30},
            {'position': [-1/8, np.sqrt(3)/8, 0], 'angle': 30},
            {'position': [4/8, np.sqrt(3)/4, 0], 'angle': 90},
            {'position': [0, np.sqrt(3)/4, 0], 'angle': 90},
            ]
        },
    'FXT': {
        'a': 1,
        'b': 1,
        'c': 3.6,
        'alpha': 90,
        'beta': 90,
        'gamma': 120,
        'vertice_connectivity': 4,
        'edge_connectivity': 2,
        'vertices': [
            {'position': [1/4, 3*np.sqrt(3)/12, 0], 'angle': -15},
            {'position': [0.5, 0, 0], 'angle': 45},
            {'position': [-1/4, 3*np.sqrt(3)/12, 0], 'angle': 15}
            ],
        'edges': []
        },
    'FXT_A': {
        'a': 2,
        'b': 2,
        'c': 3.6,
        'alpha': 90,
        'beta': 90,
        'gamma': 120,
        'vertice_connectivity': 4,
        'edge_connectivity': 2,
        'vertices': [
            {'position': [1/4, 3*np.sqrt(3)/12, 0], 'angle': -15},
            {'position': [0.5, 0, 0], 'angle': 45},
            {'position': [-1/4, 3*np.sqrt(3)/12, 0], 'angle': 15}
            ],
        'edges': [
            {'position': [22/64, 7*np.sqrt(3)/64, 0], 'angle': -30},
            {'position': [85/128, 7*np.sqrt(3)/64, 0], 'angle': 30},
            {'position': [4/8, 35*np.sqrt(3)/128, 0], 'angle': 90},
            {'position': [0, 29*np.sqrt(3)/128, 0], 'angle': 90},
            {'position': [21/128, 25*np.sqrt(3)/64, 0], 'angle': -30},
            {'position': [-21/128, 25*np.sqrt(3)/64, 0], 'angle': 30},
            ]
        },
    'DIA': {
        'a': 1,
        'b': 1,
        'c': 1,
        'alpha': 60,
        'beta': 60,
        'gamma': 60,
        'lattice': [[0, 1, 1], [1, 0, 1], [1, 1, 0]],
        'vertice_connectivity': 4,
        'edge_connectivity': 4,
        'vertices': [
            {'position': [0, 0, 0], 'angle': 55, 'align_v': [1, 1, 1]},
            {'position': [1/4, 1/4, 1/4], 'angle': -55, 'align_v': [-1, -1, -1]},
            ],
        'edges': []
        },
    'DIA_A': {
        'a': 1,
        'b': 1,
        'c': 1,
        'alpha': 60,
        'beta': 60,
        'gamma': 60,
        'lattice': [[0, 1, 1], [1, 0, 1], [1, 1, 0]],
        'vertice_connectivity': 4,
        'edge_connectivity': 2,
        'vertices': [
            {'position': [0, 0, 0], 'angle': -7.5, 'align_v': [1, 1, 1]},
            {'position': [1/4, 1/4, 1/4], 'angle': 7.5, 'align_v': [-1, -1, -1]},
            ],
        'edges': [
            {'position': [1/8, 1/8, 1/8], 'angle': -8.8, 'align_v': [1, 1, 1]},
            {'position': [1/8, 3/8, 3/8], 'angle': 16, 'align_v': [-1/4, 1/4, 1/4]},
            {'position': [3/8, 1/8, 3/8], 'angle': -78, 'align_v': [1/4, -1/4, 1/4]},
            {'position': [3/8, 3/8, 1/8], 'angle': 16, 'align_v': [1/4, 1/4, -1/4]},
            ]
        },
    'BOR': {
        'a': 1,
        'b': 1,
        'c': 1,
        'alpha': 90,
        'beta': 90,
        'gamma': 90,
        'lattice': [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        'vertice_connectivity': 4,
        'edge_connectivity': 3,
        'vertices': [
            {'position': [1/2, 0, 0], 'angle': -10, 'align_v': [-1, 1, 1]},
            {'position': [0, 1/2, 0], 'angle': 40, 'align_v': [1, -1, 1]},
            {'position': [0, 0, 1/2], 'angle': -45, 'align_v': [1, 1, -1]},
            ],
        'edges': [
            {'position': [1/6, 1/6, 1/6], 'angle': 30, 'align_v': [1, 1, 1]},
            {'position': [1/6, 5/6, 5/6], 'angle': 0, 'align_v': [1/2, 1, 1]},
            {'position': [5/6, 1/6, 5/6], 'angle': 0, 'align_v': [1, 1/2, 1]},
            {'position': [5/6, 5/6, 1/6], 'angle': 10, 'align_v': [1, 1, 1]},
            ]
        }
    }
