# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020

@author: lipel
"""

import os
import numpy as np
from scipy.spatial import distance
try:
    from pymatgen.io.cif import CifParser
except Exception:
    print('WARNING: Could no import CifParser from pymatgen the conversion from cif to xyz and COF generation may not work properlly')
    cif_parser_imported = False
import simplejson

def elements_dict(prop='mass'):
    '''Returns a dictionary containing the elements symbol and its selected property.
    Parameters
    ----------
    prop : string
        The desired property. Can be:
        - mass: Atomi mass
        - full_name : Full name of all elements
        - atomic_number: Atomic number
        - polarizability: Polarizability
        - pauling: Enetronegativity on Pauling scale
        - thermo_electronecativity: Eletronecativity on the thermodynamic scale
        - mulliken: Eletronegativity on the Mulliken scale
        - covalent_radius: Covalent radius
        - atomic_radius: Atomic radius
        - element_symbols: Elements labels
    Returns
    -------
    prop_dic : dictionary
        A dictionary containing the elements symbol and its respective property.
    '''

    mass = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.0107, 'N': 14.0067,
            'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
            'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
            'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045, 'Fe': 55.845,
            'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64,
            'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585,
            'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
            'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6,
            'I': 126.90447, 'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116,
            'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
            'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04,
            'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
            'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804,
            'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0254, 'Ac': 227.0278,
            'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 'Am': 243.0614,
            'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0951,
            'No': 259.1009, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271, 'Bh': 270, 'Hs': 269, 'Mt': 278,
            'Ds': 281, 'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294,
            'X': 0.0, 'Q': 0.0, 'R': 0.0, 'R1': 0.0, 'R2': 0.0, 'R3': 0.0, 'R4': 0.0, 'R5': 0.0, 'R6': 0.0, 'R7': 0.0,
            'R8': 0.0, 'R9': 0.0}

    full_name = {'H': 'Hydrogen', 'He': 'Helium', 'Li': 'Lithium', 'Be': 'Beryllium', 'B': 'Boron',
                 'C': 'Carbon', 'N': 'Nitrogen', 'O': 'Oxygen', 'F': 'Fluorine', 'Ne': 'Neon', 'Na': 'Sodium',
                 'Mg': 'Magnesium', 'Al': 'Aluminium', 'Si': 'Silicon', 'P': 'Phosphorus', 'S': 'Sulfur',
                 'Cl': 'Chlorine', 'Ar': 'Argon', 'K': 'Potassium', 'Ca': 'Calcium', 'Sc': 'Scandium',
                 'Ti': 'Titanium', 'V': 'Vanadium', 'Cr': 'Chromium', 'Mn': 'Manganese', 'Fe': 'Iron', 'Co': 'Cobalt',
                 'Ni': 'Nickel', 'Cu': 'Copper', 'Zn': 'Zinc', 'Ga': 'Gallium', 'Ge': 'Germanium', 'As': 'Arsenic',
                 'Se': 'Selenium', 'Br': 'Bromine', 'Kr': 'Krypton', 'Rb': 'Rubidium', 'Sr': 'Strontium',
                 'Y': 'Yttrium', 'Zr': 'Zirconium', 'Nb': 'Niobium', 'Mo': 'Molybdenum', 'Tc': 'Technetium',
                 'Ru': 'Ruthenium', 'Rh': 'Rhodium', 'Pd': 'Palladium', 'Ag': 'Silver', 'Cd': 'Cadmium',
                 'In': 'Indium', 'Sn': 'Tin', 'Sb': 'Antimony', 'Te': 'Tellurium', 'I': 'Iodine', 'Xe': 'Xenon',
                 'Cs': 'Caesium', 'Ba': 'Barium', 'La': 'Lanthanum', 'Ce': 'Cerium', 'Pr': 'Praseodymium',
                 'Nd': 'Neodymium', 'Pm': 'Promethium', 'Sm': 'Samarium', 'Eu': 'Europium', 'Gd': 'Gadolinium',
                 'Tb': 'Terbium', 'Dy': 'Dysprosium', 'Ho': 'Holmium', 'Er': 'Erbium', 'Tm': 'Thulium',
                 'Yb': 'Ytterbium', 'Lu': 'Lutetium', 'Hf': 'Hafnium', 'Ta': 'Tantalum', 'W': 'Tungsten',
                 'Re': 'Rhenium', 'Os': 'Osmium', 'Ir': 'Iridium', 'Pt': 'Platinum', 'Au': 'Gold', 'Hg': 'Mercury',
                 'Tl': 'Thallium', 'Pb': 'Lead', 'Bi': 'Bismuth', 'Po': 'Polonium', 'At': 'Astatine', 'Rn': 'Radon',
                 'Fr': 'Francium', 'Ra': 'Radium', 'Ac': 'Actinium', 'Th': 'Thorium', 'Pa': 'Protactinium', 'U': 'Uranium',
                 'Np': 'Neptunium', 'Pu': 'Plutonium', 'Am': 'Americium', 'Cm': 'Curium', 'Bk': 'Berkelium', 'Cf': 'Californium',
                 'Es': 'Einsteinium', 'Fm': 'Fermium', 'Md': 'Mendelevium', 'No': 'Nobelium', 'Lr': 'Lawrencium', 'Rf': 'Rutherfordium',
                 'Db': 'Dubnium', 'Sg': 'Seaborgium', 'Bh': 'Bohrium', 'Hs': 'Hassium', 'Mt': 'Meitnerium', 'Ds': 'Darmstadtium',
                 'Rg': 'Roentgenium', 'Cn': 'Copernicium', 'Nh': 'Nihonium', 'Fl': 'Flerovium', 'Mc': 'Moscovium',
                 'Lv': 'Livermorium', 'Ts': 'Tennessine', 'Og': 'Oganesson'}

    atomic_number = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11,
                     'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21,
                     'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31,
                     'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41,
                     'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51,
                     'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61,
                     'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71,
                     'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81,
                     'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91,
                     'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101,
                     'No': 102, 'Lr': 103, 'Rf': 104, 'Ha': 105, 'Sg': 106, 'Hs': 107, 'Bh': 108, 'Mt': 109, 'Uun': 110,
                     'X': 0.0, 'Q': 0.0, 'R': 0.0, 'R1': 0.0, 'R2': 0.0, 'R3': 0.0, 'R4': 0.0, 'R5': 0.0, 'R6': 0.0, 'R7': 0.0,
                     'R8': 0.0, 'R9': 0.0}

    # https://www.tandfonline.com/doi/full/10.1080/00268976.2018.1535143
    polarizability = {'H': 4.51, 'He': 1.38, 'Li': 164, 'Be': 37.7, 'B': 20.5, 'C': 11.3, 'N': 7.4,
                      'O': 5.3, 'F': 3.74, 'Ne': 2.66, 'Na': 163, 'Mg': 71.2, 'Al': 58, 'Si': 37.3,
                      'P': 25, 'S': 19.4, 'Cl': 14.6, 'Ar': 11.1, 'K': 290, 'Ca': 161, 'Sc': 97,
                      'Ti': 100, 'V': 87, 'Cr': 83, 'Mn': 68, 'Fe': 62, 'Co': 55, 'Ni': 49, 'Cu': 74,
                      'Zn': 38.7, 'Ga': 50, 'Ge': 40, 'As': 30, 'Se': 29, 'Br': 21, 'Kr': 16.8,
                      'Rb': 320, 'Sr': 197, 'Y': 162, 'Zr': 112, 'Nb': 98, 'Mo': 87, 'Tc': 79,
                      'Ru': 72, 'Rh': 66, 'Pd': 26.1, 'Ag': 55, 'Cd': 46, 'In': 65, 'Sn': 53, 'Sb': 43,
                      'Te': 38, 'I': 32.9, 'Xe': 27.3, 'Cs': 401, 'Ba': 272, 'La': 215, 'Ce': 205,
                      'Pr': 216, 'Nd': 208, 'Pm': 200, 'Sm': 192, 'Eu': 184, 'Gd': 158, 'Tb': 170,
                      'Dy': 165, 'Ho': 156, 'Er': 150, 'Tm': 144, 'Yb': 139, 'Lu': 137, 'Hf': 103,
                      'Ta': 74, 'W': 68, 'Re': 62, 'Os': 57, 'Ir': 54, 'Pt': 48, 'Au': 36, 'Hg': 33.9,
                      'Tl': 50, 'Pb': 47, 'Bi': 48, 'Po': 44, 'At': 42, 'Rn': 35, 'Fr': 318, 'Ra': 246,
                      'Ac': 203, 'Th': 217, 'Pa': 154, 'U': 129, 'Np': 151, 'Pu': 132, 'Am': 131,
                      'Cm': 144, 'Bk': 125, 'Cf': 122, 'Es': 118, 'Fm': 113, 'Md': 109, 'No': 110,
                      'Lr': 320, 'Rf': 112, 'Db': 42, 'Sg': 40, 'Bh': 38, 'Hs': 36, 'Mt': 34, 'Ds': 32,
                      'Rg': 32, 'Cn': 28, 'Nh': 29, 'Fl': 31, 'Mc': 70, 'Lv': 67, 'Ts': 62, 'Og': 58,
                      'X': 0.0, 'Q': 0.0, 'R': 0.0, 'R1': 0.0, 'R2': 0.0, 'R3': 0.0, 'R4': 0.0, 'R5': 0.0, 'R6': 0.0, 'R7': 0.0,
                      'R8': 0.0, 'R9': 0.0}

    pauling = {'H': 2.2, 'Li': 0.98, 'Na': 0.93, 'K': 0.82, 'Rb': 0.82, 'Cs': 0.79, 'Fr': 0.7,
               'Be': 1.57, 'Mg': 1.31, 'Ca': 1.0, 'Sr': 0.95, 'Ba': 0.89, 'Ra': 0.9, 'Sc': 1.36,
               'Ti': 1.54, 'V': 1.63, 'Cr': 1.66, 'Mn': 1.55, 'Fe': 1.83, 'Co': 1.88, 'Ni': 1.91,
               'Cu': 1.9, 'Zn': 1.65, 'Y': 1.22, 'Zr': 1.33, 'Nb': 1.6, 'Mo': 2.16, 'Tc': 1.9,
               'Ru': 2.2, 'Rh': 2.28, 'Pd': 2.2, 'Ag': 1.93, 'Cd': 1.69, 'Hf': 1.3, 'Ta': 1.5,
               'W': 2.36, 'Re': 1.9, 'Os': 2.2, 'Ir': 2.2, 'Pt': 2.28, 'Au': 2.54, 'Hg': 2.0,
               'B': 2.04, 'Al': 1.61, 'Ga': 1.81, 'In': 1.78, 'Tl': 1.62, 'C': 2.55, 'Si': 1.9,
               'Ge': 2.01, 'Sn': 1.96, 'Pb': 2.33, 'N': 3.04, 'P': 2.19, 'As': 2.18, 'Sb': 2.05,
               'Bi': 2.02, 'O': 3.44, 'S': 2.58, 'Se': 2.55, 'Te': 2.1, 'Po': 2.0, 'F': 3.98,
               'Cl': 3.16, 'Br': 2.96, 'I': 2.66, 'At': 2.2, 'La': 1.1, 'Ce': 1.12, 'Pr': 1.13,
               'Nd': 1.14, 'Pm': 1.13, 'Sm': 1.17, 'Eu': 1.2, 'Gd': 1.2, 'Tb': 1.1, 'Dy': 1.22,
               'Ho': 1.23, 'Er': 1.24, 'Tm': 1.25, 'Yb': 1.1, 'Lu': 1.27, 'Th': 1.3, 'U': 1.38,
               'Kr': 3.23, 'Xe': 3.02, 'Rn': 2.81, 'X': 0.0, 'Q': 0.0, 'R': 0.0, 'R1': 0.0,
               'R2': 0.0, 'R3': 0.0, 'R4': 0.0, 'R5': 0.0, 'R6': 0.0, 'R7': 0.0, 'R8': 0.0, 'R9': 0.0}

    # https://www.nature.com/articles/s41467-021-22429-0#Sec5
    thermo_electronecativity = {'H': 3.04, 'Li': 2.17, 'Na': 2.15, 'K': 2.07, 'Rb': 2.07, 'Cs': 1.97,
                                'Fr': 2.01, 'Be': 2.42, 'Mg': 2.39, 'Ca': 2.2, 'Sr': 2.13, 'Ba': 2.02,
                                'Sc': 2.35, 'Ti': 2.23, 'V': 2.08, 'Cr': 2.12, 'Mn': 2.2, 'Fe': 2.32,
                                'Co': 2.34, 'Ni': 2.32, 'Cu': 2.86, 'Zn': 2.26, 'Y': 2.52, 'Zr': 2.05,
                                'Nb': 2.59, 'Mo': 2.47, 'Tc': 2.82, 'Ru': 2.68, 'Rh': 2.65, 'Pd': 2.7,
                                'Ag': 2.88, 'Cd': 2.36, 'Hf': 2.01, 'Ta': 2.32, 'W': 2.42, 'Re': 2.59,
                                'Os': 2.72, 'Ir': 2.79, 'Pt': 2.98, 'Au': 2.81, 'Hg': 2.92, 'B': 3.04,
                                'Al': 2.52, 'Ga': 2.43, 'In': 2.29, 'Tl': 2.26, 'C': 3.15, 'Si': 2.82,
                                'Ge': 2.79, 'Sn': 2.68, 'Pb': 2.62, 'N': 3.56, 'P': 3.16, 'As': 3.15,
                                'Sb': 3.05, 'O': 3.78, 'S': 3.44, 'Se': 3.37, 'Te': 3.14, 'F': 4.0,
                                'Cl': 3.56, 'Br': 3.45, 'I': 3.2, 'La': 2.49, 'Ce': 2.61, 'Pr': 2.24,
                                'Nd': 2.11, 'Sm': 1.9, 'Eu': 1.81, 'Gd': 2.4, 'Tb': 2.29, 'Dy': 2.07,
                                'Ho': 2.12, 'Er': 2.02, 'Tm': 2.03, 'Yb': 1.78, 'Lu': 2.68, 'Th': 2.62,
                                'U': 2.45, 'X': 0.0, 'Q': 0.0, 'R': 0.0, 'R1': 0.0, 'R2': 0.0, 'R3': 0.0,
                                'R4': 0.0, 'R5': 0.0, 'R6': 0.0, 'R7': 0.0, 'R8': 0.0, 'R9': 0.0}

    mulliken = {'H': 7.18, 'Li': 3.0, 'Na': 2.84, 'K': 2.42, 'Rb': 2.33, 'Cs': 2.18, 'Fr': 2.21, 'Be': 4.41,
                'Mg': 3.62, 'Ca': 3.07, 'Sr': 2.87, 'Ba': 2.68, 'Ra': 2.69, 'Sc': 3.37, 'Ti': 3.45, 'V': 3.64,
                'Cr': 3.72, 'Mn': 3.46, 'Fe': 4.03, 'Co': 4.27, 'Ni': 4.4, 'Cu': 4.48, 'Zn': 4.4, 'Y': 3.26,
                'Zr': 3.53, 'Nb': 3.84, 'Mo': 3.92, 'Tc': 3.91, 'Ru': 4.2, 'Rh': 4.3, 'Pd': 4.45, 'Ag': 4.44,
                'Cd': 4.14, 'Hf': 3.5, 'Ta': 4.1, 'W': 4.4, 'Re': 3.97, 'Os': 4.89, 'Ir': 5.34, 'Pt': 5.57,
                'Au': 5.77, 'Hg': 4.97, 'B': 4.29, 'Al': 3.21, 'Ga': 3.21, 'In': 3.09, 'Tl': 3.24, 'C': 6.26,
                'Si': 4.77, 'Ge': 4.57, 'Sn': 4.23, 'Pb': 3.89, 'N': 7.23, 'P': 5.62, 'As': 5.31, 'Sb': 4.85,
                'Bi': 4.11, 'O': 7.54, 'S': 6.22, 'Se': 5.89, 'Te': 5.49, 'Po': 4.91, 'F': 10.41, 'Cl': 8.29,
                'Br': 7.59, 'I': 6.76, 'At': 5.87, 'La': 3.06, 'Ce': 3.05, 'Pr': 3.21, 'Nd': 3.72, 'Pm': 2.86,
                'Sm': 2.9, 'Eu': 2.89, 'Gd': 3.14, 'Tb': 3.51, 'Dy': 3.15, 'Ho': 3.18, 'Er': 3.21, 'Tm': 3.61,
                'Yb': 3.12, 'Lu': 2.89, 'Th': 3.63, 'U': 3.36, 'He': 12.29, 'Ne': 10.78, 'Ar': 7.88, 'Kr': 7.0,
                'Xe': 6.07, 'Rn': 5.37, 'X': 0.0, 'Q': 0.0, 'R': 0.0, 'R1': 0.0, 'R2': 0.0, 'R3': 0.0,
                'R4': 0.0, 'R5': 0.0, 'R6': 0.0, 'R7': 0.0, 'R8': 0.0, 'R9': 0.0}

    covalent_radius = {'H': 0.32, 'He': 0.93, 'Li': 1.23, 'Be': 0.9, 'B': 0.82, 'C': 0.77, 'N': 0.75, 'O': 0.73,
                       'F': 0.72, 'Ne': 0.71, 'Na': 1.54, 'Mg': 1.36, 'Al': 1.18, 'Si': 1.11, 'P': 1.06,
                       'S': 1.02, 'Cl': 0.99, 'Ar': 0.98, 'K': 2.03, 'Ca': 1.74, 'Sc': 1.44, 'Ti': 1.32,
                       'V': 1.22, 'Cr': 1.18, 'Mn': 1.17, 'Fe': 1.17, 'Co': 1.16, 'Ni': 1.15, 'Cu': 1.17,
                       'Zn': 1.25, 'Ga': 1.26, 'Ge': 1.22, 'As': 1.2, 'Se': 1.16, 'Br': 1.14, 'Kr': 1.12,
                       'Rb': 2.16, 'Sr': 1.91, 'Y': 1.62, 'Zr': 1.45, 'Nb': 1.34, 'Mo': 1.3, 'Tc': 1.27,
                       'Ru': 1.25, 'Rh': 1.25, 'Pd': 1.28, 'Ag': 1.34, 'Cd': 1.48, 'In': 1.44, 'Sn': 1.41,
                       'Sb': 1.4, 'Te': 1.36, 'I': 1.33, 'Xe': 1.31, 'Cs': 2.35, 'Ba': 1.98, 'La': 1.69,
                       'Ce': 1.65, 'Pr': 1.65, 'Nd': 1.64, 'Pm': 1.63, 'Sm': 1.62, 'Eu': 1.85, 'Gd': 1.61,
                       'Tb': 1.59, 'Dy': 1.59, 'Ho': 1.58, 'Er': 1.57, 'Tm': 1.56, 'Yb': 1.74, 'Lu': 1.56,
                       'Hf': 1.44, 'Ta': 1.34, 'W': 1.3, 'Re': 1.28, 'Os': 1.26, 'Ir': 1.27, 'Pt': 1.3,
                       'Au': 1.34, 'Hg': 1.49, 'Tl': 1.48, 'Pb': 1.47, 'Bi': 1.46, 'Po': 1.46, 'At': 1.45,
                       'Rn': 2.0, 'Fr': 2.0, 'Ra': 2.0, 'Ac': 2.0, 'Th': 1.65, 'Pa': 2.0, 'U': 1.42, 'Np': 2.0,
                       'Pu': 2.0, 'Am': 2.0, 'Cm': 2.0, 'Bk': 2.0, 'Cf': 2.0, 'Es': 2.0, 'Fm': 2.0, 'Md': 2.0,
                       'No': 2.0, 'Lr': 2.0, 'Rf': 2.0, 'Ha': 2.0, 'Sg': 2.0, 'Hs': 2.0, 'Bh': 2.0, 'Mt': 2.0,
                       'X': 0.0, 'Q': 0.0, 'R': 0.0, 'R1': 0.0, 'R2': 0.0, 'R3': 0.0, 'R4': 0.0, 'R5': 0.0, 'R6': 0.0, 'R7': 0.0,
                       'R8': 0.0, 'R9': 0.0}

    atomic_radius = {'H': 0.79, 'He': 0.49, 'Li': 2.05, 'Be': 1.4, 'B': 1.17, 'C': 0.91, 'N': 0.75, 'O': 0.65,
                     'F': 0.57, 'Ne': 0.51, 'Na': 2.23, 'Mg': 1.72, 'Al': 1.82, 'Si': 1.46, 'P': 1.23, 'S': 1.09,
                     'Cl': 0.97, 'Ar': 0.88, 'K': 2.77, 'Ca': 2.23, 'Sc': 2.09, 'Ti': 2, 'V': 1.92, 'Cr': 1.85,
                     'Mn': 1.79, 'Fe': 1.72, 'Co': 1.67, 'Ni': 1.62, 'Cu': 1.57, 'Zn': 1.53, 'Ga': 1.81,
                     'Ge': 1.52, 'As': 1.33, 'Se': 1.22, 'Br': 1.12, 'Kr': 1.03, 'Rb': 2.98, 'Sr': 2.45, 'Y': 2.27,
                     'Zr': 2.16, 'Nb': 2.08, 'Mo': 2.01, 'Tc': 1.95, 'Ru': 1.89, 'Rh': 1.83, 'Pd': 1.79, 'Ag': 1.75,
                     'Cd': 1.71, 'In': 2, 'Sn': 1.72, 'Sb': 1.53, 'Te': 1.42, 'I': 1.32, 'Xe': 1.24, 'Cs': 3.34,
                     'Ba': 2.78, 'La': 2.74, 'Ce': 2.7, 'Pr': 2.67, 'Nd': 2.64, 'Pm': 2.62, 'Sm': 2.59, 'Eu': 2.56,
                     'Gd': 2.54, 'Tb': 2.51, 'Dy': 2.49, 'Ho': 2.47, 'Er': 2.45, 'Tm': 2.42, 'Yb': 2.4, 'Lu': 2.25,
                     'Hf': 2.16, 'Ta': 2.09, 'W': 2.02, 'Re': 1.97, 'Os': 1.92, 'Ir': 1.87, 'Pt': 1.83, 'Au': 1.79,
                     'Hg': 1.76, 'Tl': 2.08, 'Pb': 1.81, 'Bi': 1.63, 'Po': 1.53, 'At': 1.43, 'Rn': 1.34, 'Fr': 2.0,
                     'Ra': 2.0, 'Ac': 2.0, 'Th': 2.0, 'Pa': 2.0, 'U': 2.0, 'Np': 2.0, 'Pu': 2.0, 'Am': 2.0, 'Cm': 2.0,
                     'Bk': 2.0, 'Cf': 2.0, 'Es': 2.0, 'Fm': 2.0, 'Md': 2.0, 'No': 2.0, 'Lr': 2.0, 'Rf': 2.0, 'Ha': 2.0,
                     'Sg': 2.0, 'Hs': 2.0, 'Bh': 2.0, 'Mt': 2.0, 'Uun': 2.0, 'X': 0.0, 'Q': 0.0, 'R': 0.0, 'R1': 0.0,
                     'R2': 0.0, 'R3': 0.0, 'R4': 0.0, 'R5': 0.0, 'R6': 0.0, 'R7': 0.0, 'R8': 0.0, 'R9': 0.0}

    element_symbols = [
    'H', 'He',  # Period 1
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',  # Period 2
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',  # Period 3
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',  # Period 4
    'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',  # Period 4
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',  # Period 5
    'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',  # Period 5
    'Cs', 'Ba',  # Period 6
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',  # Lanthanides
    'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',  # Lanthanides
    'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',  # Period 6
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',  # Period 6
    'Fr', 'Ra',  # Period 7
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',  # Actinides
    'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',  # Actinides
    'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn',  # Period 7
    'Uut', 'Fl', 'Uup', 'Lv', 'Uus', 'Uuo', # Period 7
    'X', 'Q', 'R', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9']  #Specific labels for the programm

    out = {'full_name':full_name, 'mass':mass, 'atomic_number':atomic_number, 'polarizability':polarizability, 'pauling':pauling, 
    'thermo_electronecativity':thermo_electronecativity, 'mulliken':mulliken, 'covalent_radius':covalent_radius,
    'atomic_radius':atomic_radius, 'element_symbols':element_symbols}

    return out[prop]

def unit_vector(x):
    """Return a unit vector in the same direction as x."""
    y = np.array(x, dtype='float')
    return y / np.linalg.norm(y)

def angle(v1, v2, unit='degree'):
    """
    Calculates the angle between two vectors v1 and v2.
    ----------
    v1 : array
        (N,1) matrix with N dimensions
    v2 : array
        (N,1) matrix with N dimensions
    unit : str
        Unit of the output. Could be 'degree', 'radians' or 'cos'.
    Returns
    -------
    angle : float
        Angle in the selected unit. 
    """
    unit_vector1 = unit_vector(v1)
    unit_vector2 = unit_vector(v2)

    dot_product = np.dot(unit_vector1, unit_vector2)

    if unit == 'degree':
        angle = np.arccos(dot_product) * 180. / np.pi
    if unit == 'radians':
        angle = np.arccos(dot_product)
    if unit == 'cos':
        angle = dot_product
    return angle

def rotation_matrix_from_vectors(vec1, vec2):
    '''
    Find the rotation matrix that aligns vec1 to vec2
    ----------
    vec1 : array
        (3,3) array 
    vec2 : array
        (3,3) array 
    Returns
    -------
    rotation_matrix : array
        A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    '''
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    if s != 0:
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
        return rotation_matrix
    else:
        return np.identity(3)

def translate_inside(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if matrix[i][j] >= 1:
                matrix[i][j] -= 1
    return matrix

def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    Parameters
    ----------
    V : array
        (N,D) matrix, where N is points and D is dimension.
    W : array
        (N,D) matrix, where N is points and D is dimension.
    Returns
    -------
    rmsd : float
        Root-mean-square deviation between the two vectors
    """
    diff = np.array(V) - np.array(W)
    N = len(V)
    return np.sqrt((diff * diff).sum() / N)

def cell_to_cellpar(cell, radians=False):
    """Returns the cell parameters [a, b, c, alpha, beta, gamma] 
    given a 3x3 cell matrix [v1, v2, v3]

    Angles are in degrees unless radian=True is used.
    """
    lengths = [np.linalg.norm(v) for v in cell]
    angles = []
    for i in range(3):
        j = i - 1
        k = i - 2
        ll = lengths[j] * lengths[k]
        if ll > 1e-16:
            x = np.dot(cell[j], cell[k]) / ll
            angle = 180.0 / np.pi * np.arccos(x)
        else:
            angle = 90.0
        angles.append(angle)
    if radians:
        angles = [angle * np.pi / 180 for angle in angles]
    return np.array(lengths + angles)

def cellpar_to_cell(cellpar, ab_normal=(0, 0, 1), a_direction=None):
    """Return a 3x3 cell matrix from cellpar=[a,b,c,alpha,beta,gamma].

    Angles must be in degrees.

    The returned cell is orientated such that a and b
    are normal to `ab_normal` and a is parallel to the projection of
    `a_direction` in the a-b plane.

    Default `a_direction` is (1,0,0), unless this is parallel to
    `ab_normal`, in which case default `a_direction` is (0,0,1).

    The returned cell has the vectors va, vb and vc along the rows. The
    cell will be oriented such that va and vb are normal to `ab_normal`
    and va will be along the projection of `a_direction` onto the a-b
    plane.
    Example:
    >>> cell = cellpar_to_cell([1, 2, 4, 10, 20, 30], (0, 1, 1), (1, 2, 3))
    >>> np.round(cell, 3)
    array([[ 0.816, -0.408,  0.408],
            [ 1.992, -0.13 ,  0.13 ],
            [ 3.859, -0.745,  0.745]])
    """
    if a_direction is None:
        if np.linalg.norm(np.cross(ab_normal, (1, 0, 0))) < 1e-5:
            a_direction = (0, 0, 1)
        else:
            a_direction = (1, 0, 0)

    # Define rotated X,Y,Z-system, with Z along ab_normal and X along
    # the projection of a_direction onto the normal plane of Z.
    ad = np.array(a_direction)
    Z = unit_vector(ab_normal)
    X = unit_vector(ad - np.dot(ad, Z) * Z)
    Y = np.cross(Z, X)

    # Express va, vb and vc in the X,Y,Z-system
    alpha, beta, gamma = 90., 90., 90.
    if isinstance(cellpar, (int, float)):
        a = b = c = cellpar
    elif len(cellpar) == 1:
        a = b = c = cellpar[0]
    elif len(cellpar) == 3:
        a, b, c = cellpar
    else:
        a, b, c, alpha, beta, gamma = cellpar

    # Handle orthorhombic cells separately to avoid rounding errors
    eps = 2 * np.spacing(90.0, dtype=np.float64)  # around 1.4e-14
    # alpha
    if abs(abs(alpha) - 90) < eps:
        cos_alpha = 0.0
    else:
        cos_alpha = np.cos(alpha * np.pi / 180.0)
    # beta
    if abs(abs(beta) - 90) < eps:
        cos_beta = 0.0
    else:
        cos_beta = np.cos(beta * np.pi / 180.0)
    # gamma
    if abs(gamma - 90) < eps:
        cos_gamma = 0.0
        sin_gamma = 1.0
    elif abs(gamma + 90) < eps:
        cos_gamma = 0.0
        sin_gamma = -1.0
    else:
        cos_gamma = np.cos(gamma * np.pi / 180.0)
        sin_gamma = np.sin(gamma * np.pi / 180.0)

    # Build the cell vectors
    va = a * np.array([1, 0, 0])
    vb = b * np.array([cos_gamma, sin_gamma, 0])
    cx = cos_beta
    cy = (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    cz_sqr = 1. - cx * cx - cy * cy
    assert cz_sqr >= 0
    cz = np.sqrt(cz_sqr)
    vc = c * np.array([cx, cy, cz])

    # Convert to the Cartesian x,y,z-system
    abc = np.vstack((va, vb, vc))
    T = np.vstack((X, Y, Z))
    cell = np.dot(abc, T)

    return cell

def get_fractional_to_cartesian_matrix(a, b, c, alpha, beta, gamma, angle_in_degrees=True):
    """
    Return the transformation matrix that converts fractional coordinates to
    cartesian coordinates.
    Parameters
    ----------
    a, b, c : float
        The lengths of the edges.
    alpha, gamma, beta : float
        The angles between the sides.
    angle_in_degrees : bool
        True if alpha, beta and gamma are expressed in degrees.
    Returns
    -------
    r : array_like
        The 3x3 rotation matrix. ``V_cart = np.dot(r, V_frac)``.
    """
    if angle_in_degrees:
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)
    cosa = np.cos(alpha)
    cosb = np.cos(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)
    volume = 1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
    volume = np.sqrt(volume)
    r = np.zeros((3, 3))
    r[0, 0] = a
    r[0, 1] = b * cosg
    r[0, 2] = c * cosb
    r[1, 1] = b * sing
    r[1, 2] = c * (cosa - cosb * cosg) / sing
    r[2, 2] = c * volume / sing
    return r

def get_cartesian_to_fractional_matrix(a, b, c, alpha, beta, gamma, angle_in_degrees=True):
    """
    Return the transformation matrix that converts cartesian coordinates to
    fractional coordinates.
    Parameters
    ----------
    a, b, c : float
        The lengths of the edges.
    alpha, gamma, beta : float
        The angles between the sides.
    angle_in_degrees : bool
        True if alpha, beta and gamma are expressed in degrees.
    Returns
    -------
    r : array_like
        The 3x3 rotation matrix. ``V_frac = np.dot(r, V_cart)``.
    """
    if angle_in_degrees:
        alpha = np.deg2rad(alpha)
        beta = np.deg2rad(beta)
        gamma = np.deg2rad(gamma)
    cosa = np.cos(alpha)
    cosb = np.cos(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)
    volume = 1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
    volume = np.sqrt(volume)
    r = np.zeros((3, 3))
    r[0, 0] = 1.0 / a
    r[0, 1] = -cosg / (a * sing)
    r[0, 2] = (cosa * cosg - cosb) / (a * volume * sing)
    r[1, 1] = 1.0 / (b * sing)
    r[1, 2] = (cosb * cosg - cosa) / (b * volume * sing)
    r[2, 2] = sing / (c * volume)
    return r

def get_reciprocal_vectors(cell):
    '''
    Get the reciprocal vectors of a cell given in cell parameters of cell vectors
    ----------
    cell : array
        (3,1) array for cell vectors or (6,1) array for cell parameters
    Returns
    -------
    b1 : array
        (3,1) array containing b_1 vector in the reciprocal space
    b2 : array
        (3,1) array containing b_2 vector in the reciprocal space
    b3 : array
        (3,1) array containing b_3 vector in the reciprocal space
    '''
    if len(cell) == 3:
        v1, v2, v3 = cell
    if len(cell) == 6:
        v1, v2, v3 = cellpar_to_cell(cell)

    vol = np.dot(v1, np.cross(v2, v3))
    b1 = 2*np.pi*np.cross(v2, v3)/vol
    b2 = 2*np.pi*np.cross(v3, v1)/vol
    b3 = 2*np.pi*np.cross(v1, v2)/vol

    return b1, b2, b3

def get_kgrid(cell, distance=0.3):
    '''Get the k-points grid in the reciprocal space with a given distance for a 
    cell given in cell parameters of cell vectors.
    ----------
    cell : array
        (3,1) array for cell vectors or (6,1) array for cell parameters
    distance : float
        distance between the points in the reciprocal space
    Returns
    -------
    kx : int
        Number of points in the x direction on reciprocal space
    ky : int
        Number of points in the y direction on reciprocal space
    kz : int
        Number of points in the z direction on reciprocal space
    '''
    b1, b2, b3 = get_reciprocal_vectors(cell)
    b = np.array([np.linalg.norm(b1), np.linalg.norm(b2), np.linalg.norm(b3)])
    kx = np.ceil(b[0]/distance)
    ky = np.ceil(b[1]/distance)
    kz = np.ceil(b[2]/distance)

    return kx, ky, kz

def create_CellBox(A, B, C, alpha, beta, gamma):
    """Creates the CellBox using the same expression as RASPA."""
    tempd = np.cos(alpha) - np.cos(gamma) * np.cos(beta) / np.sin(gamma)
    ax = A
    ay = 0
    az = 0
    bx = B * np.cos(gamma)
    by = B * np.sin(gamma)
    bz = 0 
    cx = C * np.cos(beta)
    cy = C * tempd
    cz = C * np.sqrt(1 - np.cos(beta) ** 2 - tempd ** 2 )
    
    CellBox = np.array([[ax, ay, az], [bx, by, bz], [cx, cy, cz]])
    
    return CellBox

def calculate_UnitCells(cell, cutoff):
    '''
    Calculate the number of unit cell repetitions so that all supercell lengths are larger than
    twice the interaction potential cut-off radius. 
    
    RASPA considers the perpendicular directions the directions perpendicular to the `ab`, `bc`,
    and `ca` planes. Thus, the directions depend on who the crystallographic vectors `a`, `b`, 
    and `c` are and the length in the perpendicular directions would be the projections 
    of the crystallographic vectors on the vectors `a x b`, `b x c`, and `c x a`. 
    (here `x` means cross product)
    ----------
    cell_matrix : array
        (3,3) cell vectors or (6,1)
    Returns
    -------
    SuperCell
        (3,1) list containg the number of repiting units in `x`, `y`, `z` directions. 
    '''

    # Make sure that the cell is in the format of cell matrix
    if len(cell) == 6:
        CellBox = cellpar_to_cell(cell)
    if len(cell) == 3:
        CellBox = cell
    
    # Pre-calculate the cross products
    axb = np.cross(CellBox[0], CellBox[1])
    bxc = np.cross(CellBox[1], CellBox[2])
    cxa = np.cross(CellBox[2], CellBox[0])
    
    # Calculates the cell volume
    V = np.dot(np.cross(CellBox[0], CellBox[1]), CellBox[2])

    # Calculate perpendicular widths
    cx = V / np.linalg.norm(bxc)
    cy = V / np.linalg.norm(cxa)
    cz = V / np.linalg.norm(axb)

    # Calculate UnitCells array
    SuperCell = np.ceil(2.0 * cutoff / np.array([cx, cy, cz])).astype(int)

    return SuperCell

def cellpar_to_lammpsbox(a, b, c, alpha, beta, gamma, angle_in_degrees=True):
    """
    Return the box parameters lx, ly, lz, xy, xz, yz for LAMMPS data input.
    Parameters
    ----------
    a, b, c : float
        The lengths of the edges.
    alpha, gamma, beta : float
        The angles between the sides.
    angle_in_degrees : bool
        True if alpha, beta and gamma are expressed in degrees.
    Returns
    -------
    r : array_like
        The 1x6 array with the box parameters 'lx', 'ly', 'lz', 'xy', 'xz', 'yz'.
    """
    if angle_in_degrees:
        alpha = alpha*(np.pi/180)
        beta = beta*(np.pi/180)
        gamma = gamma*(np.pi/180)

    lx = a
    xy = b * np.cos(gamma)
    xz = c * np.cos(beta)
    ly = np.sqrt( b**2 - xy **2)
    yz = (b * c * np.cos(alpha) - xy * xz) / ly
    lz = np.sqrt(c**2 - xz**2 - yz**2)
    
    return np.array([lx, ly, lz, xy, xz, yz])

def find_index(element, e_list):
    '''Finds the index of a given element in a list
    ----------
    element : string
        String containing the label of the element in e_list
    e_list : list
        List with the atom labels
    Returns
    ----------
    i : int
        The index of element in the e_list
    '''
    for i in range(len(e_list)):
        if np.array_equal(e_list[i], element):
            return i

def change_X_atoms(atom_labels, atom_pos, bond_atom):
    ''' Changes the X atom for the desired bond_atom or remove it if bond_atom == 'R'.
    ----------
    atom_labels : list
        List containing the atom labels
    atom_pos : list 
        List containing the atom position
    Returns
    ----------
    labels : list
        List containing the processed atom labels
    pos : list 
        List containing the processed atom position
    '''
    label, pos = [], []

    for i in range(len(atom_labels)):
        if atom_labels[i] == 'X' and bond_atom != 'R':
            label += [bond_atom]
            pos += [atom_pos[i]]
        if atom_labels[i] != 'X':
            label += [atom_labels[i]]
            pos += [atom_pos[i]]

    return label, pos 

def find_bond_atom(cof_name):
    '''Finds the type of atom that the program heve to substitute X based on the building blocks'''
    bb1, bb2, net, stacking = cof_name.split('-')
    conect_1 = bb1.split('_')[2]
    conect_2 = bb2.split('_')[2]

    bond_dict = {'NH2':'N',
                 'CONHNH2':'N',
                 'BOH2': 'B',
                 'Cl': 'R',
                 'Br':'R'}

    for group in list(bond_dict.keys()):
        if group in [conect_1, conect_2]:
            return bond_dict[group]


def closest_atom(label_1, pos_1, labels, pos):

    list_labels = []
    list_pos = []

    for i in range(len(labels)):
        if labels[i] != label_1:
            list_labels += [labels[i]]
            list_pos += [pos[i]]

    closest_index = distance.cdist([pos_1], list_pos).argmin()

    return list_labels[closest_index], list_pos[closest_index], np.linalg.norm(pos_1 - list_pos[closest_index])

def closest_atom_struc(label_1, pos_1, labels, pos):
    list_labels = []
    list_pos = []
    for i in range(len(labels)):
        if labels[i] != label_1:
            if 'C' in labels[i]:
                list_labels += [labels[i]]
                list_pos += [pos[i]]
    closest_index = distance.cdist([pos_1], list_pos).argmin()

    return list_labels[closest_index], list_pos[closest_index], np.linalg.norm(pos_1-list_pos[closest_index])

def print_result(name, lattice, hall, space_group, space_number, symm_op):
    '''Print the results of the created structure'''
    print('{:<60s} {:^12s} {:<4s} {:^4s} #{:^5s}   {:^2} sym. op.'.format(name, lattice, hall.lstrip('-'), space_group, space_number, symm_op))

def print_comand(text, verbose, match):
    if verbose in match:
        print(text)

############# Reads and save files #####################

def save_csv(path, file_name, data, delimiter=',', head=False):
    """
    Saves a file in format `.csv`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the `csv` file. Does not neet to contain the `.csv` extention. 
    data : list
        Data to be saved.
    delimiter: str
        Delimiter of the columns. `,` is the default. 
    head : str
        Names of the columns.
    """
    file_name = file_name.split('.')[0] # Remove the extention if exists
    file_name = os.path.join(path, file_name + '.csv')

    file_temp = open(file_name, 'w')
    if head is not False:
        file_temp.write(head)
    for i in range(len(data)):
        file_temp.write(delimiter.join([str(j) for j in data[i]]) + '\n')

    file_temp.close()

def read_xyz_file(path, file_name):
    """
    Reads a file in format `.xyz` from the `path` given and returns a list containg the N atom labels and 
    a Nx3 array contaning the atoms coordinates.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the `xyz` file. Does not neet to contain the `.xyz` extention. 

    Returns
    -------
    atom_labels : list
        List of strings containing containg the N atom labels. 
    atom_pos : numpy array
        Nx3 array contaning the atoms coordinates
    """
    
    file_name = file_name.split('.')[0] # Remove the extention if exists

    if os.path.exists(os.path.join(path, file_name + '.xyz')):
        temp_file = open(os.path.join(path, file_name + '.xyz'), 'r').readlines()
        
        atoms = [i.split() for i in temp_file[2:]]

        atom_labels = [i[0] for i in atoms if len(i) > 1]
        atom_pos = np.array([[float(i[1]), float(i[2]), float(i[3])] for i in atoms if len(i) > 1])

        return atom_labels, atom_pos
    else:
        print(f'File {file_name} not found!')
        return None

def read_gjf_file(path, file_name):
    """
    Reads a file in format `.gjf` from the `path` given and returns a list containg the N atom labels and 
    a Nx3 array contaning the atoms coordinates.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the `gjf` file. Does not neet to contain the `.gjf` extention. 

    Returns
    -------
    atom_labels : list
        List of strings containing containg the N atom labels. 
    atom_pos : numpy array
        Nx3 array contaning the atoms coordinates
    """

    file_name = file_name.split('.')[0] # Remove the extention if exists

    if os.path.exists(os.path.join(path, file_name + '.gjf')):
    
        temp_file = open(os.path.join(path, file_name + '.gjf'), 'r').readlines()
        temp_file = [i.split() for i in temp_file if i != '\n']

        atoms = [i for i in temp_file if i[0] in elements_dict()]

        atom_labels = [i[0] for i in atoms]
        atom_pos = np.array([[float(i[1]), float(i[2]), float(i[3])] for i in atoms])

        return atom_labels, atom_pos
    else:
        print(f'File {file_name} not found!')
        return None

def read_cif(path, file_name):
    """
    Reads a file in format `.cif` from the `path` given and returns a list containg the N atom labels and 
    a Nx3 array contaning the atoms coordinates.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the `cif` file. Does not neet to contain the `.cif` extention. 

    Returns
    -------
    cell : numpy array
        3x3 array contaning the cell vectors.
    atom_labels : list
        List of strings containing containg the N atom labels. 
    atom_pos : numpy array
        Nx3 array contaning the atoms coordinates
    charges : list
        List of strings containing containg the N atom partial charges. 
    """

    file_name = file_name.split('.')[0] # Remove the extention if exists

    if os.path.exists(os.path.join(path, file_name + '.cif')):

        temp_file = open(os.path.join(path, file_name + '.cif'), 'r').readlines()
        cell = []
        atom_label = []
        atom_pos = []
        charges = []
        has_charges = False
        for i in temp_file:
            if 'cell_length_a' in i:
                cell += [float(i.split()[-1])]
            if 'cell_length_b' in i:
                cell += [float(i.split()[-1])]    
            if 'cell_length_c' in i:
                cell += [float(i.split()[-1])]  
            if 'cell_angle_alpha' in i:
                cell += [float(i.split()[-1])]  
            if '_cell_angle_beta' in i:
                cell += [float(i.split()[-1])]  
            if '_cell_angle_gamma' in i:
                cell += [float(i.split()[-1])]  
            if '_atom_site_charge' in i:
                has_charges = True

        for i in temp_file:
            line = i.split()
            if len(line) > 1 and line[0] in elements_dict().keys():
                atom_label += [line[0]]
                atom_pos += [[float(j) for j in line[2:5]]]
                if has_charges:
                    charges += [float(line[-1])]
        cell = cellpar_to_cell(cell)

        return cell, atom_label, atom_pos, charges 
    else:
        print(f'File {file_name} not found!')
        return None   

def save_xsf(path, file_name, cell, atom_label, atom_pos):
    """
    Save a file in format `.xsf` on the `path`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.xsf` extention. 
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters. 
    atom_label : list
        List of strings containing containg the N atom label. 
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 6:
        cell = cellpar_to_cell(cell)

    xsf_file = open(os.path.join(path, file_name + '.xsf'), 'w')
    xsf_file.write(' CRYSTAL\n')
    xsf_file.write('  PRIMVEC\n')

    for i in range(len(cell)):
        xsf_file.write(f'  {cell[i][0]:>5.9f}    {cell[i][1]:>5.9f}    {cell[i][2]:>5.9f}\n')

    xsf_file.write('   PRIMCOORD\n')
    xsf_file.write(f'           {len(atom_pos)}           1\n')

    for i in range(len(atom_pos)):
        xsf_file.write(f'{atom_label[i]}        {atom_pos[i][0]:>5.9f}    {atom_pos[i][1]:>5.9f}    {atom_pos[i][2]:>5.9f}\n')

    xsf_file.close()

def save_pqr(path, file_name, cell, atom_label, atom_pos, partial_charges=False):
    """
    Save a file in format `.pqr` on the `path`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.pqr` extention. 
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters. 
    atom_label : list
        List of strings containing containg the N atom partial charges. 
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    partial_charges: list
        List of strings containing containg the N atom partial charges. 
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 3:
        cell = cell_to_cellpar(cell)

    pqr_file = open(os.path.join(path, file_name + '.pqr'), 'w')
    pqr_file.write(f'TITLE       {file_name}  \n')
    pqr_file.write('REMARK   4\n')
    pqr_file.write(f'CRYST1{cell[0]:>9.3f}{cell[1]:>9.3f}{cell[2]:>9.3f}{cell[3]:>7.2f}{cell[4]:>7.2f}{cell[5]:>7.2f} P1\n')

    if partial_charges is not False:
        for i in range(len(atom_pos)):
            pqr_file.write(f'ATOM   {i+1:>4} {atom_label[i]:>2}   MOL A   0    {atom_pos[i][0]:>8.3f}{atom_pos[i][1]:>8.3f}{atom_pos[i][2]:>8.3f}{partial_charges:>8.5f}                {atom_label[i]}\n')
    if partial_charges is False:
        for i in range(len(atom_pos)):
            pqr_file.write(f'ATOM   {i+1:>4} {atom_label[i]:>2}   MOL A   0    {atom_pos[i][0]:>8.3f}{atom_pos[i][1]:>8.3f}{atom_pos[i][2]:>8.3f}                {atom_label[i]}\n')

    pqr_file.close()

def save_pdb(path, file_name, cell, atom_label, atom_pos):
    """
    Save a file in format `.pdb` on the `path`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.pdb` extention. 
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters. 
    atom_label : list
        List of strings containing containg the N atom partial charges. 
    atom_pos : list
        Nx3 array contaning the atoms coordinates in cartesian form.
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 3:
        cell = cell_to_cellpar(cell)

    pdb_file = open(os.path.join(path, file_name + '.pdb'), 'w')
    pdb_file.write(f'TITLE       {file_name}  \n')
    pdb_file.write('REMARK   pyCOFBuilder\n')
    pdb_file.write(f'CRYST1{cell[0]:>9.3f}{cell[1]:>9.3f}{cell[2]:>9.3f}{cell[3]:>7.2f}{cell[4]:>7.2f}{cell[5]:>7.2f} P1\n')

    for i in range(len(atom_pos)):
        pdb_file.write(f'ATOM   {i+1:>4} {atom_label[i]:>2}   MOL          {atom_pos[i][0]:>8.3f}{atom_pos[i][1]:>8.3f}{atom_pos[i][2]:>8.3f}  1.00  0.00           {atom_label[i]}\n')

    pdb_file.close()

def save_gjf(path, file_name, atom_labels, atom_pos, text='opt pm6'):
    """
    Save a file in format `.gjf` on the `path`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.gjf` extention. 
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters. 
    atom_label : list
        List of strings containing containg the N atom partial charges. 
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    text : str
        Parameters for Gaussian calculations.
    """

    file_name = file_name.split('.')[0]

    temp_file = open(os.path.join(path, file_name), 'w')
    temp_file.write(f'%chk={file_name}.chk \n')
    temp_file.write(f'# {text}\n')
    temp_file.write('\n')
    temp_file.write(f'{file_name}\n')
    temp_file.write('\n')
    temp_file.write('0 1 \n')

    for i in range(len(atom_labels)):
        temp_file.write('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atom_labels[i], atom_pos[i][0], atom_pos[i][1], atom_pos[i][2]))

    temp_file.write('\n')
    temp_file.write('\n')
    temp_file.close()

def save_xyz(path, file_name, atom_labels, atom_pos, cell=None):
    """
    Save a file in format `.xyz` on the `path`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.xyz` extention. 
    atom_label : list
        List of strings containing containg the N atom partial charges. 
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters. 
    """

    file_name = file_name.split('.')[0]

    temp_file = open(os.path.join(path, file_name + '.xyz'), 'w')
    temp_file.write(f'{len(atom_labels)}\n')

    if cell is None:
        temp_file.write(f'{file_name}\n')
    else:
        if len(cell) == 3:
            cell = cell_to_cellpar(cell)
        temp_file.write(f'{cell[0]}  {cell[1]}  {cell[2]}  {cell[3]}  {cell[4]}  {cell[5]}\n')

    for i in range(len(atom_labels)):
        temp_file.write('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atom_labels[i], atom_pos[i][0], atom_pos[i][1], atom_pos[i][2]))

    temp_file.close()

def save_json(path, file_name, cell, atom_labels, atom_pos):
    """
    Save a file in format `.json` on the `path`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.cif` extention. 
    atom_label : list
        List of strings containing containg the N atom partial charges. 
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters. 
    """

    file_name = file_name.split('.')[0]

    cof_json = create_COF_json(file_name)

    if len(cell) == 3:
        cell_par = cell_to_cellpar(np.array(cell)).tolist()
        cell_par =  [round(i, 10) for i in cell_par]

    if len(cell) == 6:
        cell_par = cell
        cell = cellpar_to_cell(cell_par).tolist()

    cof_json['system']['geo_opt'] = False
    
    cof_json['geometry']['cell_matrix'] = cell
    cof_json['geometry']['cell_parameters'] = cell_par
    cof_json['geometry']['atom_labels'] = atom_labels
    cof_json['geometry']['atom_pos'] = atom_pos

    write_json(path, file_name, cof_json)

def save_cif(path, file_name, cell, atom_labels, atom_pos, partial_charges=False, frac_coords=True ):
    """
    Save a file in format `.cif` on the `path`.

    Parameters
    ----------
    path : str
        Path to the file.
    file_name : str
        Name of the file. Does not neet to contain the `.cif` extention. 
    atom_label : list
        List of strings containing containg the N atom partial charges. 
    atom_pos : list
        Nx3 array contaning the atoms coordinates.
    cell : numpy array
        Can be a 3x3 array contaning the cell vectors or a list with the 6 cell parameters. 
    """

    file_name = file_name.split('.')[0]

    if len(cell) == 3:
        a, b, c, alpha, beta, gamma = cell_to_cellpar(cell)
    if len(cell) == 6:
        a, b, c, alpha, beta, gamma = cell

    r = get_cartesian_to_fractional_matrix(a, b, c, alpha, beta, gamma)

    cif_file = open(os.path.join(path, file_name + '.cif'), 'w')

    cif_file.write(f'data_{file_name}\n')
    cif_file.write(f'_chemical_name_common                  \'{file_name}\'\n')
    cif_file.write(f'_cell_length_a                         {a:>10.6f}\n')
    cif_file.write(f'_cell_length_b                         {b:>10.6f}\n')
    cif_file.write(f'_cell_length_c                         {c:>10.6f}\n')
    cif_file.write(f'_cell_angle_alpha                      {alpha:>6.2f}\n')
    cif_file.write(f'_cell_angle_beta                       {beta:>6.2f}\n')
    cif_file.write(f'_cell_angle_gamma                      {gamma:>6.2f}\n')
    cif_file.write('_space_group_name_H-M_alt              \'P 1\'\n')
    cif_file.write('_space_group_IT_number                 1\n')
    cif_file.write('\n')
    cif_file.write('loop_\n')
    cif_file.write('_symmetry_equiv_pos_as_xyz\n')
    cif_file.write('   \'x, y, z\'\n')
    cif_file.write('\n')
    cif_file.write('loop_\n')
    cif_file.write('   _atom_site_label\n')
    cif_file.write('   _atom_site_type_symbol\n')
    cif_file.write('   _atom_site_fract_x\n')
    cif_file.write('   _atom_site_fract_y\n')
    cif_file.write('   _atom_site_fract_z\n')
    if partial_charges is not False:
        cif_file.write('   _atom_site_charge\n')

    if frac_coords is False:
        atom_pos = [np.dot(r, [i[0], i[1], i[2]]) for i in atom_pos]

    for i in range(len(atom_pos)):
        u, v, w = atom_pos[i][0], atom_pos[i][1], atom_pos[i][2]
        if partial_charges is not False:
            cif_file.write(f'{atom_labels[i]}    {atom_labels[i]} {u:>15.9f} {v:>15.9f} {w:>15.9f} {partial_charges[i]:>10.5f}\n')
        else:
            cif_file.write(f'{atom_labels[i]}    {atom_labels[i]} {u:>15.9f} {v:>15.9f} {w:>15.9f}\n')

    cif_file.close()
     
def convert_json_2_cif(origin_path, file_name, destiny_path, charge_type='None'):
    """
    Convert a file in format `.json` to `.cif`.

    Parameters
    ----------
    origin_path : str
        Path to the '.json' file.
    file_name : str
        Name of the file. Does not neet to contain the `.json` extention. 
    destiny_path : str
        path where the `.cif` file will be saved.
    """

    framework_JSON = read_json(origin_path, file_name)

    cell = framework_JSON['geometry']['cell_matrix']
    atom_labels = framework_JSON['geometry']['atom_labels']
    atom_pos = framework_JSON['geometry']['atom_pos']

    if charge_type + '_charges' in list(framework_JSON['system'].keys()):
        partial_charges = framework_JSON['geometry'][charge_type + '_charges']
    else:
        partial_charges = False

    save_cif(destiny_path, file_name, cell, atom_labels, atom_pos, partial_charges, frac_coords=False)

def convert_gjf_2_xyz(path, file_name):

    file_name = file_name.split('.')[0]

    atom_labels, atom_pos = read_gjf_file(path, file_name + '.gjf')

    save_xyz(path, file_name + '.xyz', atom_labels, atom_pos)

def convert_xyz_2_gjf(path, file_name):

    file_name = file_name.split('.')[0]

    atom_labels, atom_pos = read_xyz_file(path, file_name + '.xyz')

    save_xyz(path, file_name + '.gjf', atom_labels, atom_pos)

def convert_cif_2_xyz(path, file_name, supercell=[1, 1, 1]):

    file_name = file_name.split('.')[0]

    if cif_parser_imported is not False:

        structure = CifParser(os.path.join(path, file_name + '.cif')).get_structures(primitive=True)[0]

        structure.make_supercell([[supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, supercell[2]]])

        dict_sctructure = structure.as_dict()
        a, b, c = dict_sctructure['lattice']['a'], dict_sctructure['lattice']['b'], dict_sctructure['lattice']['c']
        alpha = round(dict_sctructure['lattice']['alpha'])
        beta = round(dict_sctructure['lattice']['beta'])
        gamma = round(dict_sctructure['lattice']['gamma'])

        atom_labels = [i['label'] for i in dict_sctructure['sites']]

        atom_pos = [i['xyz'] for i in dict_sctructure['sites']]

    if cif_parser_imported is False:
        cell, atom_labels, atom_pos, charges = read_cif(path, file_name)
        a, b, c, alpha, beta, gamma = cell

    temp_file = open(os.path.join(path, file_name + '.xyz'), 'w')
    temp_file.write(f'{len(atom_labels)} \n')

    temp_file.write(f'{a}  {b}  {c}  {alpha}  {beta}  {gamma}\n')

    for i in range(len(atom_labels)):
        temp_file.write('{:<5s}{:>15.7f}{:>15.7f}{:>15.7f}\n'.format(atom_labels[i], atom_pos[i][0], atom_pos[i][1], atom_pos[i][2]))

    temp_file.close()


########################### JSON related ##########################  


def write_json(path, name, COF_json):

    name = name.split('.')[0]

    if os.path.exists(path) is not True:
        os.mkdir(path)
    
    save_path = os.path.join(path, name + '.json')
    
    with open(save_path, 'w', encoding='utf-8') as f:
        simplejson.dump(COF_json, f, ensure_ascii=False, separators=(',', ':'), indent=2, ignore_nan=True)
        
def read_json(path, cof_name):
    
    cof_path = os.path.join(path, cof_name + '.json')
    
    with open(cof_path, 'r') as r:
        json_object = simplejson.loads(r.read())
    
    return json_object
    
def create_COF_json(name):

    system_info = 'Informations about the system such as name, if it is optimized and other relevant information.'
    geometry_info = 'Informations about the geometry: cell parameters, cell matrix, atomic positions, partial charges, bond orders, simmetry information'
    optimization_info = 'Information about the optimization process such as level of calculations, optimization schema and optimization steps.'
    adsorption_info = 'Information about the adsorption simulation experiments on RASPA2'
    textural_info = 'Information about the textural calculations of the structure such as specific area, pore volume, void fraction.'
    spectrum_info = 'Information about spectra simulation like DRX, FTIR, ssNMR, UV-VIS, Band dispersion, Phonon dispersion...'
    experimental_info = 'Experimental data DRX, FTIR, ssNMR, UV-VIS...'

    COF_json = {'system':{'description':system_info,
                           'name':name,
                            'geo_opt':False,
                            'execution_times_seconds':{}},
                'geometry':{'description':geometry_info},
                'optimization':{'description':optimization_info},
                'adsorption':{'description':adsorption_info},
                'textural':{'description':textural_info},
                'spectrum':{'description':spectrum_info},
                'experimental':{'description':experimental_info}
                }

    return COF_json
