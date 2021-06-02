# -*- coding: utf-8 -*-
import os
import glob

import pycofbuilder.tools as Tools
from pycofbuilder.reticulum import Reticulum
from pycofbuilder.building_block import Building_Block


def build(cof_name=None, save_format=['cif'], bond_atom='N', lib='bb_lib'):
    '''Create a COF with a given name'''

    qe = False
    xyz = False
    cif = False
    turbomole = False
    vasp = False

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
    
    bb1, bb2, net, stacking = cof_name.split('-')

    if net == 'HCB':
        try:
            Ret = Reticulum(bb_lib=lib)
            Ret.create_hcb_structure(bb1, bb2, stack=stacking, bond_atom=bond_atom)
            if cif is True:
                Ret.save_cif()
            if xyz is True:
                if stacking == 'AA':
                    Ret.save_xyz(supercell=[1, 1, 2])
                else:
                    Ret.save_xyz()
            if qe is True:
                Ret.save_qe()
            if turbomole is True:
                Ret.save_turbomole()
            if vasp is True:
                Ret.save_vasp()
            return  [True, f'{bb1}-{bb2}-{net}-{stacking}']
        except Exception:
            return [False, f'{bb1}-{bb2}-{net}-{stacking}']

    if net == 'HCB_A':
            Ret = Reticulum(bb_lib=lib)
            Ret.create_hcb_a_structure(bb1, bb2, stack=stacking, bond_atom=bond_atom)
            if cif is True:
                Ret.save_cif()
            if xyz is True:
                if stacking == 'AA':
                    Ret.save_xyz(supercell=[1, 1, 2])
                else:
                    Ret.save_xyz()
            if qe is True:
                Ret.save_qe()
            if turbomole is True:
                Ret.save_turbomole()
            if vasp is True:
                Ret.save_vasp()
            return  [True, f'{bb1}-{bb2}-{net}-{stacking}']


def build_all_available_COFs(lib='bb_lib', stacking='AA', qe=False, xyz=False, cif=True, turbomole=False, vasp=False):

    save_f = []

    for i in [[qe, 'qe'], [xyz, 'xyz'], [cif, 'cif'], [turbomole, 'turbomole'], [vasp, 'vasp']]:
        if i[0] is True:
            save_f += [i[1]]

    BB = Building_Block(lib=lib)

    lista_amina_2 = BB.get_bipodal_NH2()
    lista_amina_3 = BB.get_tripodal_NH2()
    lista_amina_4 = BB.get_tetrapodal_squared_NH2()
    lista_oh_2 = BB.get_bipodal_OH2()
    lista_oh_3 = BB.get_tripodal_OH2()
    lista_aldeido_2 = BB.get_bipodal_CHO()
    lista_aldeido_3 = BB.get_tripodal_CHO()
    lista_aldeido_4 = BB.get_tetrapodal_squared_CHO()
    lista_b_3 = BB.get_tripodal_BOH2()
    lista_b_2 = BB.get_bipodal_BOH2()
    lista_b_4 = BB.get_tetrapodal_squared_BOH2()

    cofs_list = []

    failed_list = []
    sucess_list = []

    for file_a in lista_aldeido_3:
        for file_b in lista_amina_2:
            cofs_list += [f'{file_a}-{file_b}-HCB_A-{stacking}']

    for file_a in lista_aldeido_3:
        for file_b in lista_amina_3:
            cofs_list += [f'{file_a}-{file_b}-HCB-{stacking}']

    for file_a in lista_amina_3:
        for file_b in lista_aldeido_2:
            cofs_list += [f'{file_a}-{file_b}-HCB_A-{stacking}']

    '''for file_b in lista_b_2:
        try:
            Ret = Reticulum(bb_lib=lib)
            Ret.create_hcb_a_structure('BDBA_1', file_b, stack=stacking, bond_atom='B')
            if cif is True:
                Ret.save_cif()
            if xyz is True:
                if stacking == 'AA':
                    Ret.save_xyz(supercell=[1, 1, 2])
                else:
                    Ret.save_xyz()
            if qe is True:
                Ret.save_qe()
            if turbomole is True:
                Ret.save_turbomole()
            if vasp is True:
                Ret.save_vasp()
            sucess_list += [f'BDBA_1_{file_b}']
        except Exception:
            failed_list += [f'BDBA_1_{file_b}']

    Ret = Reticulum(bb_lib=lib)
    Ret.create_hcb_a_structure('BDBA_1', 'BDBA_2', stack='AB1', bond_atom='B')
    if cif is True:
        Ret.save_cif()
    if xyz is True:
        Ret.save_xyz(supercell=[1, 1, 2])
    if qe is True:
        Ret.save_qe()
    if turbomole is True:
        Ret.save_turbomole()
    if vasp is True:
        Ret.save_vasp()

    for file_a in lista_b_4:
        for file_b in lista_oh_2:
            try:
                Ret = Reticulum(bb_lib=lib)
                Ret.create_sql_a_structure(file_a, file_b, stack=stacking, bond_atom='B')
                if cif is True:
                    Ret.save_cif()
                if xyz is True:
                    if stacking == 'AA':
                        Ret.save_xyz(supercell=[1, 1, 2])
                    else:
                        Ret.save_xyz()
                if qe is True:
                    Ret.save_qe()
                if turbomole is True:
                    Ret.save_turbomole()
                if vasp is True:
                    Ret.save_vasp()
                sucess_list += [f'{file_a}_{file_b}']
            except Exception:
                failed_list += [f'{file_a}_{file_b}']

    for file_a in lista_oh_3:
        for file_b in lista_b_2:
            try:
                Ret = Reticulum(bb_lib=lib)
                Ret.create_hcb_a_structure(file_a, file_b, stack=stacking, bond_atom='B')
                if cif is True:
                    Ret.save_cif()
                if xyz is True:
                    if stacking == 'AA':
                        Ret.save_xyz(supercell=[1, 1, 2])
                    else:
                        Ret.save_xyz()
                if qe is True:
                    Ret.save_qe()
                if turbomole is True:
                    Ret.save_turbomole()
                if vasp is True:
                    Ret.save_vasp()
                sucess_list += [f'{file_a}_{file_b}']
            except Exception:
                failed_list += [f'{file_a}_{file_b}']

    for file_a in lista_b_3:
        for file_b in lista_oh_2:
            try:
                Ret = Reticulum(bb_lib=lib)
                Ret.create_hcb_a_structure(file_a, file_b, stack=stacking, bond_atom='B')
                if cif is True:
                    Ret.save_cif()
                if xyz is True:
                    if stacking == 'AA':
                        Ret.save_xyz(supercell=[1, 1, 2])
                    else:
                        Ret.save_xyz()
                if qe is True:
                    Ret.save_qe()
                if turbomole is True:
                    Ret.save_turbomole()
                if vasp is True:
                    Ret.save_vasp()
                sucess_list += [f'{file_a}_{file_b}']
            except Exception:
                failed_list += [f'{file_a}_{file_b}']

    for file_a in lista_oh_3:
        for file_b in lista_b_3:
            try:
                Ret = Reticulum(bb_lib=lib)
                Ret.create_hcb_structure(file_a, file_b, stack=stacking, bond_atom='B')
                if cif is True:
                    Ret.save_cif()
                if xyz is True:
                    if stacking == 'AA':
                        Ret.save_xyz(supercell=[1, 1, 2])
                    else:
                        Ret.save_xyz()
                if qe is True:
                    Ret.save_qe()
                if turbomole is True:
                    Ret.save_turbomole()
                if vasp is True:
                    Ret.save_vasp()
                sucess_list += [f'{file_a}_{file_b}']
            except Exception:
                failed_list += [f'{file_a}_{file_b}']

    for file_a in lista_amina_4:
        for file_b in lista_aldeido_2:
            try:
                Ret = Reticulum(bb_lib=lib)
                Ret.create_sql_a_structure(file_a, file_b, stack=stacking)
                if cif is True:
                    Ret.save_cif()
                if xyz is True:
                    if stacking == 'AA':
                        Ret.save_xyz(supercell=[1, 1, 2])
                    else:
                        Ret.save_xyz()
                if qe is True:
                    Ret.save_qe()
                if turbomole is True:
                    Ret.save_turbomole()
                if vasp is True:
                    Ret.save_vasp()
                sucess_list += [f'{file_a}_{file_b}']
            except Exception:
                failed_list += [f'{file_a}_{file_b}']

    for file_a in lista_aldeido_4:
        for file_b in lista_amina_2:
            try:
                Ret = Reticulum(bb_lib=lib)
                Ret.create_sql_a_structure(file_a, file_b, stack=stacking)
                if cif is True:
                    Ret.save_cif()
                if xyz is True:
                    if stacking == 'AA':
                        Ret.save_xyz(supercell=[1, 1, 2])
                    else:
                        Ret.save_xyz()
                if qe is True:
                    Ret.save_qe()
                if turbomole is True:
                    Ret.save_turbomole()
                if vasp is True:
                    Ret.save_vasp()
                sucess_list += [f'{file_a}_{file_b}']
            except Exception:
                failed_list += [f'{file_a}_{file_b}']

    for file_a in lista_aldeido_4:
        for file_b in lista_amina_4:
            try:
                Ret = Reticulum(bb_lib=lib)
                Ret.create_sql_structure(file_a, file_b, stack=stacking)
                if cif is True:
                    Ret.save_cif()
                if xyz is True:
                    if stacking == 'AA':
                        Ret.save_xyz(supercell=[1, 1, 2])
                    else:
                        Ret.save_xyz()
                if qe is True:
                    Ret.save_qe()
                if turbomole is True:
                    Ret.save_turbomole()
                if vasp is True:
                    Ret.save_vasp()
                sucess_list += [f'{file_a}_{file_b}']
            except Exception:
                failed_list += [f'{file_a}_{file_b}']'''

    val = input(f'{len(cofs_list)} COFs will be created. Do you want o proceed? Type [y] to continue.\n')

    if val == 'y':
        print('                      COF Name                              |    Lattice    | Point Group | N° of symmetry op. |')
        for cof in cofs_list:
            succes, name = build(cof, save_format=save_f)
            if succes is True:
                sucess_list += [name]
            if succes is False:
                failed_list += [name]


        print(f'{len(sucess_list)} sucessful. {len(failed_list)} failled ({100*len(sucess_list)/(len(failed_list) + len(sucess_list))} % success rate)\n')
        if len(failed_list) > 0:
            print('Failed list:')
            for i in failed_list:
                print(i)
    else:
        print('Exiting...')


def creat_all_C2():
    nucleos = [i.rstrip('.gjf') for i in os.listdir(os.path.join('data', 'nucleo', 'C2'))]
    radicais = [i.rstrip('.gjf') for i in os.listdir(os.path.join('data', 'radical'))]
    conectores = [i.rstrip('.gjf') for i in os.listdir(os.path.join('data', 'conector'))]

    nucleos = ['BENZ', 'NAPT', 'BPNY', 'ANTR', 'TPNY', 'DPBY', 'PYRN', 'BPYB', 'DPEY']
    radicais = ['H']
    conectores = ['NH2']
    for n in nucleos:
        for c in conectores:
            for r1 in radicais:
                BB = Building_Block()
                BB.create_C2_BB(n, c, r1)
                print(BB.name, 'created')
                BB.save()


def creat_all_C3():

    nucleos = [i.rstrip('.gjf') for i in os.listdir(os.path.join('data', 'nucleo', 'C3'))]
    radicais = [i.rstrip('.gjf') for i in os.listdir(os.path.join('data', 'radical'))]
    conectores = [i.rstrip('.gjf') for i in os.listdir(os.path.join('data', 'conector'))]

    nucleos = ['BENZ', 'TPBZ', 'TPOB', 'DICZ']
    radicais = ['H']
    conectores = ['NH2', 'CHO']
    for n in nucleos:
        for c in conectores:
            for r1 in radicais:
                BB = Building_Block()
                BB.create_C3_BB(n, c, r1)
                print(BB.name, 'created')
                BB.save()


def creat_all_C4():

    nucleos = [i.rstrip('.gjf') for i in os.listdir(os.path.join('data', 'nucleo', 'C4'))]
    radicais = [i.rstrip('.gjf') for i in os.listdir(os.path.join('data', 'radical'))]
    conectores = [i.rstrip('.gjf') for i in os.listdir(os.path.join('data', 'conector'))]

    radicais = ['H']
    conectores = ['CHO']
    for n in nucleos:
        for c in conectores:
            for r1 in radicais:
                print(n, c, r1)
                BB = Building_Block()
                BB.create_C4_BB(n, c, r1)
                print(BB.name, 'created')
                BB.save()

def clean_bb_list():
    #Loop Through the folder projects all files and deleting them one by one
    for file in glob.glob(os.path.join('data', 'bb_lib', '*')):
        os.remove(file)
        print("Deleted " + str(file))

