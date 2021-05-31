# -*- coding: utf-8 -*-


def create_all_available_COFs(lib='built', stacking='AA', qe=True, xyz=False, cif=True, turbomole=False, vasp=False):

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

    failed_list = []
    sucess_list = []

    for file_a in lista_aldeido_3:
        for file_b in lista_amina_2:
            try:
                Ret = Reticulum(bb_lib=lib)
                Ret.create_hcb_a_structure(file_a, file_b, stack=stacking)
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

    for file_a in lista_aldeido_3:
        for file_b in lista_amina_3:
            try:
                Ret = Reticulum(bb_lib=lib)
                Ret.create_hcb_structure(file_a, file_b, stack=stacking)
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

    for file_a in lista_amina_3:
        for file_b in lista_aldeido_2:
            try:
                Ret = Reticulum(bb_lib=lib)
                Ret.create_hcb_a_structure(file_a, file_b, stack=stacking)
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

    for file_b in lista_b_2:
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

    '''Ret = Reticulum(bb_lib=lib)
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
        Ret.save_vasp()'''

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
                failed_list += [f'{file_a}_{file_b}']

    print(f'{len(sucess_list)} sucessful. {len(failed_list)} failled ({100*len(failed_list)/(len(failed_list) + len(sucess_list))} %)\n')
    if len(failed_list) > 0:
        print('Failed list:',failed_list)

'''
def creat_all_C2():

    nucleos = [i.rstrip('.gjf') for i in os.listdir(os.path.join(os.getcwd(), 'Nucleo', 'C2'))]
    conectores = [i.rstrip('.gjf') for i in os.listdir(os.path.join(os.getcwd(), 'Conectores'))]
    radicais = [i.rstrip('.gjf') for i in os.listdir(os.path.join(os.getcwd(), 'Radicais'))]
    conectores = ['CHO', 'NH2']
    for n in nucleos:
        for c in conectores:
            for r1 in radicais:
                BB = Building_Block()
                BB.create_C2_BB(n, c, r1)
                print(BB.name, 'created')
                BB.save_BB()


def creat_all_C3():

    nucleos = [i.rstrip('.gjf') for i in os.listdir(os.path.join(os.getcwd(), 'Nucleo', 'C3'))]
    conectores = [i.rstrip('.gjf') for i in os.listdir(os.path.join(os.getcwd(), 'Conectores'))]
    radicais = [i.rstrip('.gjf') for i in os.listdir(os.path.join(os.getcwd(), 'Radicais'))]
    conectores = ['CHO', 'NH2']
    for n in nucleos:
        for c in conectores:
            for r1 in radicais:
                BB = Building_Block()
                BB.create_C3_BB(n, c, r1)
                print(BB.name, 'created')
                BB.save_BB()'''


def creat_all_C2():

    nucleos = ['BENZ', 'NAPT', 'BPNY', 'ANTR', 'TPNY', 'DPBY', 'PYRN', 'BPYB', 'DPEY'] #[i.rstrip('.gjf') for i in os.listdir(os.path.join(os.getcwd(), 'Nucleo', 'C2'))]
    radicais = [i.rstrip('.gjf') for i in os.listdir(os.path.join(os.getcwd(), 'Radicais'))]
    conectores = [i.rstrip('.gjf') for i in os.listdir(os.path.join(os.getcwd(), 'Conectores'))]
    radicais = ['H']
    conectores = ['NH2']
    for n in nucleos:
        for c in conectores:
            for r1 in radicais:
                BB = Building_Block()
                BB.create_C2_BB(n, c, r1)
                print(BB.name, 'created')
                BB.save_BB()


def creat_all_C3():

    nucleos = ['BENZ', 'TPBZ', 'TPOB', 'DICZ'] #[i.rstrip('.gjf') for i in os.listdir(os.path.join(os.getcwd(), 'Nucleo', 'C3'))]
    radicais = ['H'] #[i.rstrip('.gjf') for i in os.listdir(os.path.join(os.getcwd(), 'Radicais'))]
    conectores = ['NH2', 'CHO'] #[i.rstrip('.gjf') for i in os.listdir(os.path.join(os.getcwd(), 'Conectores'))]
    for n in nucleos:
        for c in conectores:
            for r1 in radicais:
                BB = Building_Block()
                BB.create_C3_BB(n, c, r1)
                print(BB.name, 'created')
                BB.save_BB()


def creat_all_C4():

    nucleos = [i.rstrip('.gjf') for i in os.listdir(os.path.join(os.getcwd(), 'Nucleo', 'C4'))]
    radicais = [i.rstrip('.gjf') for i in os.listdir(os.path.join(os.getcwd(), 'Radicais'))]
    radicais = ['H']
    conectores = ['CHO']
    for n in nucleos:
        for c in conectores:
            for r1 in radicais:
                print(n, c, r1)
                BB = Building_Block()
                BB.create_C4_BB(n, c, r1)
                print(BB.name, 'created')
                BB.save_BB()

