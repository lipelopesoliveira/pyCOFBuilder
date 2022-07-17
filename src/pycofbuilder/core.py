# -*- coding: utf-8 -*-
import os
import glob
import time
from tqdm import tqdm

import pycofbuilder.tools as Tools
from pycofbuilder.reticulum import Reticulum
from pycofbuilder.building_block import Building_Block

_ROOT = os.path.abspath(os.path.dirname(__file__))


def build(cof_name=None,
          save_format=['json'],
          lib='bb_lib',
          print_result=True,
          supercell=[1, 1, 2],
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
    lib : str
        Library of building block.
    print_result : bool
        Boolean to print in the screen or not the result of the creation.
    supercell : list
        List containg the units or repetition for supercell creation. Default: [1,1,2]
    '''
    bond_atom = Tools.find_bond_atom(cof_name)

    qe = False
    xyz = False
    cif = False
    turbomole = False
    vasp = False
    json = False
    xsf = False
    pdb = False

    if type(save_format) is not list:
        save_format = [save_format]

    for i in save_format:
        if i == 'qe':
            qe = True
        if i == 'xyz':
            xyz = True
        if i == 'cif':
            cif = True
        if i == 'turbomole':
            turbomole = True
        if i == 'vasp':
            vasp = True
        if i == 'json':
            json = True
        if i == 'xsf':
            xsf = True
        if i == 'pdb':
            pdb = True

    bb1, bb2, net, stacking = cof_name.split('-')
    print(bb1, bb2, net, stacking, bond_atom)

    if net == 'HCB':
        try:
            Ret = Reticulum(bb_lib=lib, out_dir=save_dir, verbosity=verbosity)
            simm_data = Ret.create_hcb_structure(bb1,
                                                 bb2,
                                                 stacking=stacking,
                                                 bond_atom=bond_atom,
                                                 print_result=print_result)
            if cif is True:
                Ret.save_cif(supercell)
            if xyz is True:
                if stacking == 'AA':
                    Ret.save_xyz(supercell)
                else:
                    Ret.save_xyz()
            if qe is True:
                Ret.save_qe(supercell)
            if json is True:
                Ret.save_json(supercell)
            if xsf is True:
                Ret.save_xsf(supercell)
            if pdb is True:
                Ret.save_pdb(supercell)
            if turbomole is True:
                Ret.save_turbomole(supercell)
            if vasp is True:
                Ret.save_vasp()
            return [True, simm_data]
        except Exception:
            return [False, f'{bb1}-{bb2}-{net}-{stacking}']

    if net == 'HCB_A':
        try:
            Ret = Reticulum(out_dir=save_dir, verbosity=verbosity)
            simm_data = Ret.create_hcb_a_structure(bb1,
                                                   bb2,
                                                   stacking=stacking,
                                                   bond_atom=bond_atom,
                                                   print_result=print_result)
            if cif is True:
                Ret.save_cif()
            if xyz is True:
                if stacking == 'AA':
                    Ret.save_xyz(supercell)
                else:
                    Ret.save_xyz()
            if qe is True:
                Ret.save_qe()
            if json is True:
                Ret.save_json(supercell)
            if xsf is True:
                Ret.save_xsf()
            if pdb is True:
                Ret.save_pdb()
            if turbomole is True:
                Ret.save_turbomole()
            if vasp is True:
                Ret.save_vasp()
            return [True, simm_data]
        except Exception:
            return [False, f'{bb1}-{bb2}-{net}-{stacking}']

    if net == 'SQL':
        try:
            Ret = Reticulum(bb_lib=lib, out_dir=save_dir, verbosity=verbosity)
            simm_data = Ret.create_sql_structure(
                bb1, bb2, stacking=stacking, bond_atom=bond_atom, print_result=print_result)
            if cif is True:
                Ret.save_cif()
            if xyz is True:
                if stacking == 'AA':
                    Ret.save_xyz(supercell)
                else:
                    Ret.save_xyz()
            if qe is True:
                Ret.save_qe()
            if json is True:
                Ret.save_json(supercell)
            if xsf is True:
                Ret.save_xsf()
            if pdb is True:
                Ret.save_pdb()
            if turbomole is True:
                Ret.save_turbomole()
            if vasp is True:
                Ret.save_vasp()
            return [True, simm_data]
        except Exception:
            return [False, f'{bb1}-{bb2}-{net}-{stacking}']

    if net == 'SQL_A':
        try:
            Ret = Reticulum(bb_lib=lib, out_dir=save_dir, verbosity=verbosity)
            simm_data = Ret.create_sql_a_structure(
                bb1, bb2, stacking=stacking, bond_atom=bond_atom, print_result=print_result)
            if cif is True:
                Ret.save_cif()
            if xyz is True:
                if stacking == 'AA':
                    Ret.save_xyz(supercell)
                else:
                    Ret.save_xyz()
            if qe is True:
                Ret.save_qe()
            if json is True:
                Ret.save_json(supercell)
            if xsf is True:
                Ret.save_xsf()
            if pdb is True:
                Ret.save_pdb()
            if turbomole is True:
                Ret.save_turbomole()
            if vasp is True:
                Ret.save_vasp()
            return [True, simm_data]
        except Exception:
            return [False, f'{bb1}-{bb2}-{net}-{stacking}']


def build_all_available_COFs(lib='bb_lib', stacking='AA', qe=False, xyz=False, cif=True, turbomole=False, vasp=False, json=True):

    save_f = []

    for i in [[qe, 'qe'],
              [xyz, 'xyz'],
              [cif, 'cif'],
              [turbomole, 'turbomole'],
              [vasp, 'vasp'],
              [json, 'json']]:
        if i[0] is True:
            save_f += [i[1]]

    BB = Building_Block(lib=lib)

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
            succes, name = build(cof, save_format=save_f, print_result=False)
            if succes is True:
                sucess_list += [name]
            if succes is False:
                failed_list += [name]

        print('                      COF Name                              |    Lattice    | Point Group | N° of symmetry op. |')
        for s in sucess_list:
            Tools.print_result(*s)

        print(f'{len(sucess_list)} sucessful. {len(failed_list)} failled ({100*len(sucess_list)/(len(failed_list) + len(sucess_list)):.2f} % success rate)')
        print(f'Enlapsed time: {time.time() - t_i:.3f} s \n')
        if len(failed_list) > 0:
            print('Failed list:')
            for i in failed_list:
                print(i)
    else:
        print('Exiting...')
        
def build_COFs_list(cofs_list, save_format=['json'], lib='bb_lib', supercell=[1, 1, 2]):

    failed_list = []
    sucess_list = []

    val = input(
        f'{len(cofs_list)} COFs will be created. Do you want to proceed? Type [y] to continue.\n')

    if val == 'y':
        t_i = time.time()
        for cof in tqdm(cofs_list):
            succes, name = build(cof, save_format=save_format, print_result=False, supercell=supercell)
            if succes is True:
                sucess_list += [name]
            if succes is False:
                failed_list += [name]

        print('                      COF Name                              |    Lattice    | Point Group | N° of symmetry op. |')
        for s in sucess_list:
            Tools.print_result(*s)

        print(f'{len(sucess_list)} sucessful. {len(failed_list)} failled ({100*len(sucess_list)/(len(failed_list) + len(sucess_list)):.2f} % success rate)')
        print(f'Enlapsed time: {time.time() - t_i:.3f} s \n')
        if len(failed_list) > 0:
            print('Failed list:')
            for i in failed_list:
                print(i)
    else:
        print('Exiting...')


def create_all_C2(nucleos=None, radicais=None, conectores=None):
    '''Creates a set of C2 symmetry building block based on your choice of a nucleo, a radical group, and a type of connector group. 

    Be warned that the building blocks created only had one radical group at the position R1.  

    The creation of blocks with more than one group or in specific positions must be done manually.

    For exemple, the code below will create a `C2` building block based on a `Antracene` core with a `NH2` connection group and a `OH` group in the position `R2`:

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

    if nucleos == None:
        nucleos = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOT, 'data', 'nucleo', 'C2'))]
    if radicais == None:
        radicais = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOT, 'data', 'radical'))]
    if conectores == None:
        conectores = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOT, 'data', 'conector'))]

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
    if nucleos == None:
        nucleos = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOT, 'data', 'nucleo', 'C3'))]
    if radicais == None:
        radicais = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOT, 'data', 'radical'))]
    if conectores == None:
        conectores = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOT, 'data', 'conector'))]

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

    if nucleos == None:
        nucleos = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOT, 'data', 'nucleo', 'C4'))]
    if radicais == None:
        radicais = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOT, 'data', 'radical'))]
    if conectores == None:
        conectores = [i.rstrip('.gjf') for i in os.listdir(
            os.path.join(_ROOT, 'data', 'conector'))]

    for n in nucleos:
        for c in conectores:
            for r in radicais:
                print(n, c, r)
                BB = Building_Block()
                BB.create_C4_BB(n, c, r)
                print(BB.name, 'created')
                BB.save()

def clean_bb_dir():
    # Loop Through all building block files and deleting them one by one
    for file in glob.glob(os.path.join(_ROOT, 'data', 'bb_lib', '*')):
        os.remove(file)
        print(f'Deleted {file}')


def clean_out_dir():
    # Loop Through all output files and deleting them one by one
    for file in glob.glob(os.path.join(os.getcwd(), 'out', '*')):
        os.remove(file)
        print(f'Deleted {file}')

def clean():
    '''Clean both the building block dir and the out dir'''
    val = input(
        'This action will delete all building blocks and COFs created, do you want to proceed? Type [y] to continue.\n')
    if val == 'y':
        clean_bb_dir()
        clean_out_dir()
    else:
        print(f'{val} pressed. Nothing to be done.')
