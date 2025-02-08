# -*- coding: utf-8 -*-
# Created by Felipe Lopes de Oliveira
# Distributed under the terms of the MIT License.

"""
The dictionary containing the definitions nets for a Framework buiding. For each net, the the vertices and edges
positions are defined in cartesian coordinates divided by the lattice parameters. This definition is similar to
the alat defition used in QuantumESPRESSO, however instead of using the same alat for all dimentions, it is used
the lattice parameters a, b and c. This definition may not be the best one, however it allows a much more
intuitive definition of the vertices and edges positions.
The angles are defined in degrees.
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
        'a': 1,
        'b': 1,
        'c': 3.6,
        'alpha': 90,
        'beta': 90,
        'gamma': 90,
        'vertice_connectivity': 4,
        'edge_connectivity': 0,
        'vertices': [
            {'position': [0, 0, 0], 'angle': 45},
            {'position': [1/2, 1/2, 0], 'angle': 45}
            ],
        'edges': []
        },
    'SQL_A': {
        'a': 1,
        'b': 1,
        'c': 3.6,
        'alpha': 90,
        'beta': 90,
        'gamma': 90,
        'vertice_connectivity': 4,
        'edge_connectivity': 2,
        'vertices': [
            {'position': [0, 0, 0], 'angle': -45},
            {'position': [1/2, 1/2, 0], 'angle': -45}
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
        'edge_connectivity': 0,
        'vertices': [
            {'position': [0.5, 0, 0], 'angle': 0},
            {'position': [1/4, np.sqrt(3)/4, 0], 'angle': -60},
            {'position': [-1/4, np.sqrt(3)/4, 0], 'angle': 60}
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
            {'position': [0.5, 0, 0], 'angle': 0},
            {'position': [1/4, np.sqrt(3)/4, 0], 'angle': -60},
            {'position': [-1/4, np.sqrt(3)/4, 0], 'angle': 60}
            ],
        'edges': [
            {'position': [0.341475, 0.19715, 0], 'angle': -30},
            {'position': [0.6758525, 0.19715, 0], 'angle': 30},
            {'position': [0.5, 0.471724, 0], 'angle': 90},
            {'position': [0, 0.394301, 0], 'angle': 90},
            {'position': [0.1558525, 0.6668875, 0], 'angle': -30},
            {'position': [-0.158525, 0.668875, 0], 'angle': 30},
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
        },
    'LON': {
        'a': 1,
        'b': 1,
        'c': 2 * np.sqrt(2) / np.sqrt(3),
        'alpha': 90,
        'beta': 90,
        'gamma': 120,
        'lattice': [[1, 0, 0],
                    [-1/2, np.sqrt(3) / 2, 0],
                    [0, 0, 2 * np.sqrt(2) / np.sqrt(3)]],
        'vertice_connectivity': 4,
        'edge_connectivity': 4,
        'vertices': [
            {'position': [0.00000000, 0.57735027, 0.10206206], 'angle': 180, 'align_v': [0, 0, 1]},
            {'position': [0.50000000, 0.28867513, 0.91855856], 'angle': 0, 'align_v': [0, 0, 1]},
            {'position': [0.00000000, 0.57735027, 0.71443444], 'angle': 0, 'align_v': [0, 0, -1]},
            {'position': [0.50000000, 0.28867513, 1.53093094], 'angle': 180, 'align_v': [0, 0, -1]},
            ],
        'edges': []
        },
    'LON_A': {
        'a': 1,
        'b': 1,
        'c': 2 * np.sqrt(2) / np.sqrt(3),
        'alpha': 90,
        'beta': 90,
        'gamma': 120,
        'lattice': [[1, 0, 0], [-1/2, np.sqrt(3) / 2, 0], [0, 0, 2 * np.sqrt(2) / np.sqrt(3)]],
        'vertice_connectivity': 4,
        'edge_connectivity': 2,
        'vertices': [
            {'position': [0.00000000, 0.57735027, 0.10206206], 'angle': 180, 'align_v': [0, 0, 1]},
            {'position': [0.50000000, 0.28867513, 0.91855856], 'angle': 0, 'align_v': [0, 0, 1]},
            {'position': [0.00000000, 0.57735027, 0.71443444], 'angle': -20, 'align_v': [0, 0, -1]},
            {'position': [0.50000000, 0.28867513, 1.53093094], 'angle': 160, 'align_v': [0, 0, -1]},
            ],
        'edges': [
            {
                'position': [-0.00000,  0.57735,  0.40825],
                'angle': 90.00, 'align_v': [0, 0, -1]
            },
            {
                'position': [0.50000, 0.28867, 1.22474],
                'angle': 90.00, 'align_v': [0, 0, 1]
            },
            {
                'position': [0.25000, 0.43301, 0.00000],
                'angle': 61.87, 'align_v': [-0.8165, 0.4714, 0.33333]
            },
            {
                'position': [0.50000, 0.00000, 0.00000],
                'angle': 160.53, 'align_v': [-0.0, -0.94281, 0.33333]
            },
            {
                'position': [-0.25000, 0.43301, 0.00000],
                'angle': 118.13, 'align_v': [-0.8165, -0.47141, -0.33333]
            },
            {
                'position': [0.25000, 0.43301, 0.81650],
                'angle': 61.87, 'align_v': [-0.8165, 0.4714, -0.33333]
            },
            {
                'position': [0.50000, 0.00000, 0.81650],
                'angle': 160.53, 'align_v':  [-0.0, -0.94281, -0.33333]
            },
            {
                'position': [-0.25000, 0.43301, 0.81650],
                'angle': 118.13, 'align_v': [-0.8165, -0.47141, 0.33333]
            },
        ]
        },
    }
