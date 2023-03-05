# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020

@author: Felipe Lopes de Oliveira
"""

import os
import glob
import time
from tqdm import tqdm

import pycofbuilder.tools as Tools
from pycofbuilder.reticulum import Reticulum
from pycofbuilder.building_block import Building_Block


_ROOTDIR = os.path.abspath(os.path.dirname(__file__))


def build(cof_name=None,
          save_format='json',
          print_result=True,
          supercell=[1, 1, 1],
          save_dir=None,
          verbosity=False):
    '''Build a COF with a given name

    Parameters
    ----------
    cof_name : str
        Name of the COF to be build.
    save_format : list
        List containg the formats to save the file.
        Can be `json`, `cif`, `xyz`, `turbomole`, `vasp`, `xsf`, `pdb`.
    print_result : bool
        Boolean to print in the screen or not the result of the creation.
    supercell : list
        List containg the units or repetition for supercell creation.
        Default: [1,1,1]
    save_dir : str
        Path to the directory where the file will be saved.
    verbosity : bool
        Boolean to print in the screen or not the result of the creation.
    '''
    bond_atom = Tools.find_bond_atom(cof_name)

    if isinstance(save_format, str):
        save_format = [save_format]

    Ret = Reticulum(out_dir=save_dir, verbosity=verbosity)

    save_dict = {'qe': Ret.save_qe,
                 'xyz': Ret.save_xyz,
                 'cif': Ret.save_cif,
                 'turbomole': Ret.save_turbomole,
                 'vasp': Ret.save_vasp,
                 'json': Ret.save_json,
                 'xsf': Ret.save_xsf,
                 'pdb': Ret.save_pdb
                 }

    # Check if save_format elements is in save_dict
    error_msg = f'Save format not recognized. Available formats: {list(save_dict.keys())}'
    assert all([i in save_dict for i in save_format]), error_msg

    net_build_dict = {'HCB': Ret._create_hcb_structure,
                      'HCB_A': Ret._create_hcb_a_structure,
                      'SQL': Ret._create_sql_structure,
                      'SQL_A': Ret._create_sql_a_structure,
                      'KGM': Ret._create_kgm_structure,
                      'KGM_A': Ret._create_kgm_a_structure,
                      'KGD': Ret._create_kgd_structure,
                      'HXL_A': Ret._create_hxl_a_structure}

    bb1, bb2, net, stacking = cof_name.split('-')

    # Check if net is in net_build_dict
    error_msg = f'Net not recognized. Available nets: {list(net_build_dict.keys())}'
    assert net in net_build_dict, error_msg

    # Try to build the net
    try:
        simm_data = net_build_dict[net](bb1,
                                        bb2,
                                        stacking=stacking,
                                        bond_atom=bond_atom,
                                        print_result=print_result)

        # Save in all formats requesteds
        for i in save_format:
            save_dict[i](supercell)

        return [True, simm_data]

    except Exception as exception:
        print(exception)
        return [False, f'{bb1}-{bb2}-{net}-{stacking}']


def build_COFs_list(cofs_list,
                    save_format='cif',
                    supercell=[1, 1, 1],
                    save_dir='.',
                    print_result=False):
    '''Build all the COF structures with a given name on a list

    Parameters
    ----------
    cofs_list : str
        List with the names of the COFs to be build.
    save_format : list
        List containg the formats to save the file.
        Can be `json`, `cif`, `xyz`, `turbomole`, `vasp`, `xsf`, `pdb`.
    print_result : bool
        Boolean to print in the screen or not the result of the creation.
    supercell : list
        List containg the units or repetition for supercell creation.
        Default: [1,1,1]
    save_dir : str
        Path to the directory where the file will be saved.
    '''

    failed_list = []
    sucess_list = []

    val = input(
        f'{len(cofs_list)} COFs will be created. Do you want to proceed? Type [y] to continue.\n')

    if val == 'y':
        t_i = time.time()
        for cof in tqdm(cofs_list):
            succes, name = build(cof,
                                 save_format=save_format,
                                 print_result=print_result,
                                 save_dir=save_dir,
                                 supercell=supercell)
            if succes is True:
                sucess_list += [name]
            if succes is False:
                failed_list += [name]

        print(f'{"COF Name:^60"}|    Lattice    | Point Group | N° of symmetry op. |')
        for s in sucess_list:
            Tools.print_result(*s)

        print('{} sucessful. {} failled ({:.2f} % success rate)'.format(
            len(sucess_list),
            len(failed_list),
            100*len(sucess_list)/(len(failed_list) + len(sucess_list))
            )
            )
        print(f'Enlapsed time: {time.time() - t_i:.3f} s \n')
        if len(failed_list) > 0:
            print('Failed list:')
            for i in failed_list:
                print(i)
    else:
        print('Exiting...')


def create_all_C2(nucleos=None, radicais=None, conectores=None):
    '''Creates a set of C2 symmetry building block based on your choice of a nucleo,
    a radical group, and a type of connector group.

    Be warned that the building blocks created only had one radical group at the position R1.

    The creation of blocks with more than one group or in specific positions must be done manually.

    For exemple, the code below will create a `C2` building block based on a `Antracene` core
    with a `NH2` connection group and a `OH` group in the position `R2`:

    >>> BB = Building_Block()

    >>> BB.create_C4_BB('ANTR', 'NH2', *['H', 'OH'])

    ----------
    nucleos : list
        List containing the desired cores for the creation of blocks.
        Ex.: ['2BPD', '3BPD', 'ANTR', 'BBTZ', 'BENZ', 'BPNY', 'BPYB', 'BTPH',
              'DBTP', 'DHPI', 'DHSI', 'DPBY', 'DPDA', 'DPEL', 'DPEY', 'HDZ,
              'NAPT', 'PYRN', 'TPNY', 'TTPH']
    conectores : list
        List containing the connector groups desired for the creation of blocks.
        Ex.: ['NH2', 'CHO', 'BOH2', 'OH2', 'Cl', 'Br', 'NHNH2']
    radicais : list
        List containing the desired radical groups for the creation of blocks.
        Ex.: ['CH3', 'CHO', 'CN', 'COOH', 'F', 'H', 'NC2H', 'NH2',
        'NO2', 'O', 'OEt', 'OH', 'OMe', 'SO2H', 'tBut']
    '''

    if nucleos is None:
        nucleos = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOTDIR, 'data', 'nucleo', 'C2'))]
    if radicais is None:
        radicais = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOTDIR, 'data', 'radical'))]
    if conectores is None:
        conectores = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOTDIR, 'data', 'conector'))]

    for n in nucleos:
        for c in conectores:
            for r in radicais:
                BB = Building_Block()
                BB.create_C2_BB(n, c, r)
                print(BB.name, 'created')
                BB.save()


def create_all_C3(nucleos=None, radicais=None, conectores=None):
    '''Creates a set of C3 symmetry building block based on your choice of a nucleo, a radical
    group, and a type of connector group.

    Be warned that the building blocks created only had one radical group at the position R1.

    The creation of blocks with more than one group or in specific positions must be done manually.

    For exemple, the code below will create a `C4` building block based on a `Triphenilene`
    core with a `NH2` connection group and a `OH` group in the position `R2`:

    >>> BB = Building_Block()

    >>> BB.create_C4_BB('TPNY', 'NH2', *['H', 'OH'])

    ----------
    nucleos : list
        List containing the desired cores for the creation of blocks.
        Ex.: ['BENZ', 'DICZ', 'TPAM', 'TPBZ', 'TPNY', 'TPOB', 'TPTA', 'TPTZ']
    conectores : list
        List containing the connector groups desired for the creation of blocks.
        Ex.: ['NH2', 'CHO', 'BOH2', 'OH2', 'Cl', 'Br', 'NHNH2']
    radicais : list
        List containing the desired radical groups for the creation of blocks.
        Ex.: ['CH3', 'CHO', 'CN', 'COOH', 'F', 'H', 'NC2H', 'NH2',
        'NO2', 'O', 'OEt', 'OH', 'OMe', 'SO2H', 'tBut']
    '''
    if nucleos is None:
        nucleos = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOTDIR, 'data', 'nucleo', 'C3'))]
    if radicais is None:
        radicais = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOTDIR, 'data', 'radical'))]
    if conectores is None:
        conectores = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOTDIR, 'data', 'conector'))]

    for n in nucleos:
        for c in conectores:
            for r in radicais:
                BB = Building_Block()
                BB.create_C3_BB(n, c, r)
                print(BB.name, 'created')
                BB.save()


def create_all_C4(nucleos=None, conectores=None, radicais=None):
    '''
    Creates a set of C4 symmetry building block based on your choice of a core,W
    a radical group, and a type of connector group.

    Be warned that the building blocks created only had one radical group at the position R1.

    The creation of blocks with more than one group or in specific positions must be done manually.

    For exemple, the code below will create a `C4` building block based on a `Porphirin` core with
    a `NH2` connection group and a `OH` group in the position `R2`:

    >>> BB = Building_Block()

    >>> BB.create_C4_BB('PORP', 'NH2', ['H', 'OH', 'H'])

    ----------
    nucleos : list
        List containing the desired cores for the creation of blocks.
        Ex.: ['BENZ', 'PHPR', 'PORP', 'PYRN']
    conectores : list
        List containing the connector groups desired for the creation of blocks.W
        Ex.: ['NH2', 'CHO', 'BOH2', 'OH2', 'Cl', 'Br', 'NHNH2']
    radicais : list
        List containing the desired radical groups for the creation of blocks.
        Ex.: ['CH3', 'CHO', 'CN', 'COOH', 'F', 'H', 'NC2H', 'NH2', 'NO2',
        'O', 'OEt', 'OH', 'OMe', 'SO2H', 'tBut']
    '''

    if nucleos is None:
        nucleos = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOTDIR, 'data', 'nucleo', 'C4'))]
    if radicais is None:
        radicais = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOTDIR, 'data', 'radical'))]
    if conectores is None:
        conectores = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOTDIR, 'data', 'conector'))]

    for n in nucleos:
        for c in conectores:
            for r in radicais:
                print(n, c, r)
                BB = Building_Block()
                BB.create_C4_BB(n, c, r)
                print(BB.name, 'created')
                BB.save()


def build_all_available_COFs(stacking='AA',
                             save_format='cif'):
    '''Build all available COFs in the out directory.

    Parameters
    ----------
    stacking : str
        Stacking type.
        Default: 'AA'
    save_format : str
        Format to save the file.
        Default: 'cif'
    '''

    BB = Building_Block()

    lista_amina_2 = BB.get_bipodal_NH2()
    lista_amina_3 = BB.get_tripodal_NH2()
    lista_amina_4 = BB.get_tetrapodal_squared_NH2()
    lista_oh_2 = BB.get_bipodal_OH2()
    lista_oh_3 = BB.get_tripodal_OH2()
    lista_oh_4 = BB.get_tetrapodal_squared_OH2()
    lista_aldeido_2 = BB.get_bipodal_CHO()
    lista_aldeido_3 = BB.get_tripodal_CHO()
    lista_aldeido_4 = BB.get_tetrapodal_squared_CHO()
    lista_b_2 = BB.get_bipodal_BOH2()
    lista_b_3 = BB.get_tripodal_BOH2()
    lista_b_4 = BB.get_tetrapodal_squared_BOH2()

    cofs_list = []

    failed_list = []
    sucess_list = []

    for file_a in lista_aldeido_3:
        for file_b in lista_amina_2:
            cofs_list += [f'{file_a}-{file_b}-HCB_A-{stacking}']

    for file_a in lista_aldeido_4:
        for file_b in lista_amina_2:
            cofs_list += [f'{file_a}-{file_b}-SQL_A-{stacking}']

    for file_a in lista_aldeido_3:
        for file_b in lista_amina_3:
            cofs_list += [f'{file_a}-{file_b}-HCB-{stacking}']

    for file_a in lista_amina_3:
        for file_b in lista_aldeido_2:
            cofs_list += [f'{file_a}-{file_b}-HCB_A-{stacking}']

    for file_a in lista_amina_4:
        for file_b in lista_aldeido_2:
            cofs_list += [f'{file_a}-{file_b}-SQL_A-{stacking}']

    for file_a in lista_b_3:
        for file_b in lista_oh_2:
            cofs_list += [f'{file_a}-{file_b}-HCB_A-{stacking}']

    for file_a in lista_b_3:
        for file_b in lista_oh_3:
            cofs_list += [f'{file_a}-{file_b}-HCB-{stacking}']

    for file_a in lista_oh_3:
        for file_b in lista_b_2:
            cofs_list += [f'{file_a}-{file_b}-HCB_A-{stacking}']

    for file_a in lista_oh_4:
        for file_b in lista_b_2:
            cofs_list += [f'{file_a}-{file_b}-SQL_A-{stacking}']

    for file_a in lista_b_4:
        for file_b in lista_oh_2:
            cofs_list += [f'{file_a}-{file_b}-SQL_A-{stacking}']

    val = input(
        f'{len(cofs_list)} COFs will be created. Do you want to proceed? Type [y] to continue.\n')

    if val == 'y':
        t_i = time.time()
        for cof in tqdm(cofs_list):
            succes, name = build(cof,
                                 save_format=save_format,
                                 print_result=False)
            if succes is True:
                sucess_list += [name]
            if succes is False:
                failed_list += [name]

        print('                      COF Name                              ',
              '    Lattice    ',
              ' Point Group ',
              ' N° of symmetry op.',
              sep='|')
        print('-'*35)
        for s in sucess_list:
            Tools.print_result(*s)

        print('{:i} sucessful. {:i} failled ({:.2f} % success rate)'.format(
            len(sucess_list),
            len(failed_list),
            100*len(sucess_list)/(len(failed_list) + len(sucess_list))))

        print(f'Enlapsed time: {time.time() - t_i:.3f} s \n')
        if len(failed_list) > 0:
            print('Failed list:')
            for i in failed_list:
                print(i)
    else:
        print('Exiting...')


def clean_out_dir():
    # Loop Through all output files and deleting them one by one
    for file in glob.glob(os.path.join(os.getcwd(), 'out', '*')):
        os.remove(file)
        print(f'Deleted {file}')
