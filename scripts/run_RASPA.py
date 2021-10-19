# -*- coding: utf-8 -*-
"""
Created on Mon Oct 4 12:51:19 2020
@author: Felipe Lopes de Oliveira
"""
import tools as Tools
import os
from textwrap import dedent
from tqdm import tqdm

def get_pseudoatoms(molecule):
    '''Returns the pseudoatoms of a given molecule. If the molecule is not in the supported list will returns `None`. 
    Parameters
    ----------
    molecule : string
        Molecule name. Could be CO2, N2, O2, or H2O.
    Returns
    ----------
    pseudoatoms : list
        List containing the strings with the pseudotoms.
    '''

    pseudoatoms_dict = {'CO2': ['C_co2', 'O_co2'],
                        'N2': ['N_n2'],
                        'O2': ['O_o2'],
                        'H2O': ['Ow', 'Hw', 'Lw']}

    if molecule in list(pseudoatoms_dict.keys()):
        return pseudoatoms_dict[molecule]
    else:
        return None


def create_grid_inp(path, FrameworkName, 
                    GasComposition={'CO2':1.0}, 
                    Unit_Cell=[2,2,8],
                    SpacingVDWGrid=0.05,
                    SpacingCoulombGrid=0.05):

    GRID_DICT = {'FrameworkName': FrameworkName,
                 'UnitCells': ' '.join(map(str, Unit_Cell)),
                 'SpacingVDWGrid': SpacingVDWGrid,
                 'SpacingCoulombGrid': SpacingCoulombGrid,
                 'NumberOfGrids': 0,
                 'GridTypes': ''}   

    for gas in list(GasComposition.keys()):
        pseudo_atoms = get_pseudoatoms(gas)
        GRID_DICT['NumberOfGrids'] += len(pseudo_atoms)           # int
        GRID_DICT['GridTypes'] += ' '.join(pseudo_atoms) + ' '    # string

    Grid_InputFile = dedent("""\
        SimulationType                      MakeGrid                                            

        ForceField                          DREIDING-UFF-TRAPPE                                   

        Framework                           0         
        FrameworkName                       {FrameworkName}                                       
        UseChargesFromCIFFile               yes                                     
        UnitCells                           {UnitCells}  

        NumberOfGrids                       {NumberOfGrids}
        GridTypes                           {GridTypes}
        SpacingVDWGrid                      {SpacingVDWGrid}
        SpacingCoulombGrid                  {SpacingCoulombGrid}                     

        """).format(**GRID_DICT)

    # Write input to file
    with open(os.path.join(path, 'simulation-MakeGrid.input'), 'w') as f:
        f.write(Grid_InputFile)


def create_raspa_inp(path, FrameworkName, HeliumVoidFraction=0.0, ExternalTemperature=77, 
                     Ext_Press=[1000,10000,100000], Unit_Cell=[2,2,8], GasComposition={'CO2':1.0},
                     Grid=True, Movies=True, WriteMoviesEvery=100, VTK=True):
    '''Create the RASPA GCMC simulation input.
    Parameters
    ----------
    path : string
        Path where the file will be saved.
    FrameworkName : string
        Name of the structure. Must be the same name in the `.cif` file. 
    HeliumVoidFraction : float
        Helium Void Fraction of the structure. Is needed to calculate the excess adsorption. The default is 0.0.
    ExternalTemperature : float
        External temperature in Kelvin. Currently only one value is possible.
    Ext_Press : int or list
        External pressure in Pascal. Can be a integer value or a list of integers. 
    Unit_Cell : list
        Size of the supercell used for the calculation. 
    GasComposition : dict
        A dictionary containing the gas composition of the simulation. e.g. {'CO2':1.0} or {'CO2':0.2, 'N2':0.8}
    Grid : `True` or `False`
        Use grid calculation to speed up simulation of isotherms. Only helps if more than 1 point of pressure is used. 
    Movies : `True` or `False`
        Salve a `.pdb` file containg the atomic posioton at every `WriteMoviesEver` 
    WriteMoviesEvery : int
        Frequency to save the movie. 
    VTK : `True` or `False`
        Save the adsorption density in a `VTK` file.
    '''

    # Calculation parameters dictionary
    CALC_DICT = {'NumberOfCycles': 10000,                          # int
                'NumberOfInitializationCycles': 1000,              # int
                'PrintEvery': 100,                                 # int
                'PrintPropertiesEvery': 0,                         # int
                'CutOffVDW': 12.8,                                 # float
                'CutOffChargeCharge': 12.8,                        # float
                'CutOffChargeBondDipole': 12.8,                    # float
                'CutOffBondDipoleBondDipole': 12.8,                # float
                'EwaldPrecision': 1.0e-6,                          # float
                'FrameworkName': FrameworkName,                    # string
                'HeliumVoidFraction': HeliumVoidFraction,          # float
                'ExternalTemperature': ExternalTemperature,        # int
                'ExternalPressure': ' '.join(map(str, Ext_Press)), # float or csv
                'UseChargesFromCif': 'yes',                        # yes / no
                'UnitCells': ' '.join(map(str, Unit_Cell))}        # int int int

    GRID_DICT = {'Use':Grid,
                'SpacingVDWGrid': 0.075,                            # float
                'SpacingCoulombGrid': 0.075,                        # float
                'UseTabularGrid': 'yes',                           # yes / no
                'NumberOfGrids': 0,                                # int
                'GridTypes': ''}                                   # string

    MOVIES_DICT = {'Use': Movies,
                        'WriteMoviesEvery': WriteMoviesEvery }     # int

    VTK_DENSITY_DICT = {'Use': VTK,
                        'DensityProfile3DVTKGridPoints': '100 100 100',  # int int int
                        'WriteDensityProfile3DVTKGridEvery': 100}        # int

    # Create file header as string
    GCMC_InputFile = dedent("""\
    SimulationType                      MonteCarlo
    NumberOfCycles                      {NumberOfCycles}                        
    NumberOfInitializationCycles        {NumberOfInitializationCycles}          
    PrintEvery                          {PrintEvery}                            
    PrintPropertiesEvery                {PrintPropertiesEvery} 
    PrintForcefieldToOutput             no          
    PrintPseudoAtomsToOutput            no      
    PrintMoleculeDefinitionToOutput     no 

    ContinueAfterCrash                  yes                                     
    WriteBinaryRestartFileEvery         500                                    

    ForceField                          DREIDING-UFF-TRAPPE                                   
    CutOffVDW                           {CutOffVDW}                          
    CutOffChargeCharge                  {CutOffChargeCharge}
    CutOffChargeBondDipole              {CutOffChargeBondDipole}
    CutOffBondDipoleBondDipole          {CutOffBondDipoleBondDipole}
    ChargeMethod                        Ewald          
    EwaldPrecision                      {EwaldPrecision}  

    Framework                           0         
    FrameworkName                       {FrameworkName}       
    HeliumVoidFraction                  {HeliumVoidFraction}                   
    ExternalTemperature                 {ExternalTemperature}                   
    ExternalPressure                    {ExternalPressure}                     
    UseChargesFromCIFFile               yes                                     
    UnitCells                           {UnitCells}                             

    """).format(**CALC_DICT)

    if GRID_DICT['Use'] is True:

        for gas in list(GasComposition.keys()):
            pseudo_atoms = get_pseudoatoms(gas)
            GRID_DICT['NumberOfGrids'] += len(pseudo_atoms)           # int
            GRID_DICT['GridTypes'] += ' '.join(pseudo_atoms) + ' '    # string

        GCMC_InputFile += dedent("""\
        NumberOfGrids                       {NumberOfGrids} 
        GridTypes                           {GridTypes}                           
        SpacingVDWGrid                      {SpacingVDWGrid}        
        SpacingCoulombGrid                  {SpacingCoulombGrid}  
        UseTabularGrid                      {UseTabularGrid}                        

        """).format(**GRID_DICT)

    if VTK_DENSITY_DICT['Use'] is True:
        GCMC_InputFile += dedent("""\
        ComputeDensityProfile3DVTKGrid      yes                                    
        DensityProfile3DVTKGridPoints       {DensityProfile3DVTKGridPoints}         
        WriteDensityProfile3DVTKGridEvery   {WriteDensityProfile3DVTKGridEvery}

        """).format(**VTK_DENSITY_DICT)

    if MOVIES_DICT['Use'] is True:
        GCMC_InputFile += dedent("""\
        Movies yes                                  
        WriteMoviesEvery       {WriteMoviesEvery}

        """).format(**MOVIES_DICT)              


    # Create component list as string
    for name, fraction in GasComposition.items():

        number_of_components = len(GasComposition)
        index_of_component = list(GasComposition).index(name)

        # Append component string block to input file
        GCMC_InputFile += dedent(f"""\
        Component {index_of_component} MoleculeName                  {name}
                    MolFraction                   {fraction}
                    MoleculeDefinition            TraPPE
                    SwapProbability               0.5
                    TranslationProbability        0.5
                    RotationProbability           0.5
                    IdentityChangeProbability     0.5
                        NumberOfIdentityChanges     {number_of_components}
                        IdentityChangesList         {' '.join(map(str, range(number_of_components)))}
                    CreateNumberOfMolecules       0

        """)

    #Write input to file
    with open(os.path.join(path, 'simulation-MonteCarlo.input'), 'w') as f:
        f.write(GCMC_InputFile)


def low_pressure_isotherm(DataBasePath, OutPath, FrameworkName, ExternalTemperature=298, 
                          GasComposition={'N2':0.8, 'CO2':0.2}, ExternalPressures='default', 
                          grid=False, charge_type='EQeq'):
    '''Create the RASPA GCMC simulation input for low pressure (from 0 to 1 bar) isotherm.
    Parameters
    ----------
    DataBasePath : str
        Path containg the COFs Database JSON files.
    OutPath : string
        Path where the file will be saved.
    FrameworkName : string
        Name of the structure. Must be the same name in the `.cif` file. 
    ExternalTemperature : float
        External temperature in Kelvin. Currently only one value is possible.
    GasComposition : dict
        A dictionary containing the gas composition of the simulation. e.g. {'CO2':1.0} or {'CO2':0.2, 'N2':0.8}
    ExternalPressures : int or list
        External pressure in Pascal. Can be a integer value or a list of integers. `default` uses a pre-set of pressures.
    '''

    Ex_press_list = [1000,10000,25000,50000,101250]

    COF_JSON = Tools.read_json(DataBasePath, FrameworkName)

    cell = COF_JSON['geometry']['cell_matrix']
    atom_labels = COF_JSON['geometry']['atom_labels']
    atom_pos = COF_JSON['geometry']['atom_pos']

    charge_dict = {'eqeq':'EQeq_charges',
                   'ddec':'DDEC_charges',
                   'lowdin': 'lowdin_charges',
                   'hirshfeld': 'hirshfeld_charges',
                   'mulliken': 'mulliken_charges'}

    charge_label = charge_dict[charge_type.lower()]

    if COF_JSON['system'][charge_label] is True:
        partial_charges = COF_JSON['geometry'][charge_label]
    else:
        partial_charges = False

    if COF_JSON['system']['textural'] is True:
        HeliumVoidFraction = COF_JSON['textural']['acessible_volume_fraction']
    else:
        HeliumVoidFraction = 0.0

    Tools.save_cif(OutPath, FrameworkName, cell, atom_labels, atom_pos, partial_charges, frac_coords=False)

    Unit_Cell = Tools.calculate_UnitCells(cell, 12.8)

    if grid is True:
        create_grid_inp(OutPath, FrameworkName, GasComposition=GasComposition)

    create_raspa_inp(OutPath, FrameworkName, HeliumVoidFraction, ExternalTemperature, 
                     GasComposition=GasComposition, Ext_Press=Ex_press_list, Unit_Cell=Unit_Cell, 
                     Grid=grid, Movies=False, VTK=False)

def single_pressure(DataBasePath, OutPath, FrameworkName, ExternalTemperature=298, 
                    GasComposition={'N2':0.8, 'CO2':0.2}, ExternalPressures='default', 
                    grid=False, charge_type='EQeq'):
    '''Create the RASPA GCMC simulation input for low pressure (from 0 to 1 bar) isotherm.
    Parameters
    ----------
    DataBasePath : str
        Path containg the COFs Database JSON files.
    OutPath : string
        Path where the file will be saved.
    FrameworkName : string
        Name of the structure. Must be the same name in the `.cif` file. 
    ExternalTemperature : float
        External temperature in Kelvin. Currently only one value is possible.
    GasComposition : dict
        A dictionary containing the gas composition of the simulation. e.g. {'CO2':1.0} or {'CO2':0.2, 'N2':0.8}
    ExternalPressures : int or list
        External pressure in Pascal. Can be a integer value or a list of integers. `default` uses a pre-set of pressures.
    '''
    if ExternalPressures == 'default':
        Ex_press = [101250]
    else:
        Ex_press = ExternalPressures

    COF_JSON = Tools.read_json(DataBasePath, FrameworkName)

    cell = COF_JSON['geometry']['cell_matrix']
    atom_labels = COF_JSON['geometry']['atom_labels']
    atom_pos = COF_JSON['geometry']['atom_pos']

    charge_dict = {'eqeq':'eqeq_charge',
                   'ddec':'DDEC_charges',
                   'lowdin': 'lowdin_charges',
                   'hirshfeld': 'hirshfeld_charges',
                   'mulliken': 'mulliken_charges'}

    charge_label = charge_dict[charge_type.lower()]
    
    if COF_JSON['system'][charge_label] is True:
        partial_charges = COF_JSON['geometry'][charge_label + 's']
    else:
        partial_charges = False

    Tools.save_cif(OutPath, FrameworkName, cell, atom_labels, atom_pos, partial_charges, frac_coords=False)

    if COF_JSON['system']['textural'] is True:
        HeliumVoidFraction = COF_JSON['textural']['acessible_volume_fraction']
    else:
        HeliumVoidFraction = 0.0

    Unit_Cell = Tools.calculate_UnitCells(cell, 12.8)

    if grid is True:
        create_grid_inp(OutPath, FrameworkName, GasComposition=GasComposition)

    create_raspa_inp(OutPath, FrameworkName, HeliumVoidFraction, ExternalTemperature, 
                     GasComposition=GasComposition, Ext_Press=Ex_press, Unit_Cell=Unit_Cell, 
                     Grid=grid, Movies=False, VTK=False)


DataBasePath = '/home/felipelopes/PythonCodes/pyCOFBuilder/COF_Database'
OutPath = '/home/felipelopes/Simulations/Teste2'

#FrameworkName = 'C3_BENZ_CHO_H-C2_3BPD_NH2_H_H_H-HCB_A-AA'
ROOT_PATH = os.getcwd()

######################################
COF_DATABASE_PATH = '/home/felipe/Simulations/Size_Matter/'
OUT_PATH = '/home/felipe/Simulations/Size_EQeq'
######################################

RASPA_CMD = '/home/felipe/Programas/RASPA/simulations/bin/simulate'

COF_list = [i for i in os.listdir(COF_DATABASE_PATH) if '.json' in i]

for cof in tqdm(COF_list):

    FrameworkName = cof.rstrip('.json')
    try:
        os.mkdir(os.path.join(OUT_PATH, FrameworkName))
    except:
        None
    try:
        os.mkdir(os.path.join(OUT_PATH, FrameworkName, 'RASPA'))
    except:
        None
    try:
        os.mkdir(os.path.join(OUT_PATH, FrameworkName, 'RASPA', '298'))
    except:
        None
    try:
        os.mkdir(os.path.join(OUT_PATH, FrameworkName, 'RASPA', '298', 'CO2'))
    except:
        None

    os.chdir(os.path.join(OUT_PATH, FrameworkName, 'RASPA', '298', 'CO2'))

    single_pressure(COF_DATABASE_PATH, os.path.join(OUT_PATH, FrameworkName, 'RASPA', '298', 'CO2'), 
                   FrameworkName, ExternalTemperature=298, GasComposition={'CO2':1.0}, 
                   ExternalPressures='default', grid=False)

    os.system(f'{RASPA_CMD} -i simulation-MonteCarlo.input >/dev/null')

    os.chdir(ROOT_PATH)