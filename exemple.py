# -*- coding: utf-8 -*-

import pycofbuilder as COF

COF_list = ['C3_TPBZ_NH2_H-C2_BENZ_CHO_H-HCB_A-AA',
            'C3_TPBZ_NH2_H-C2_BPNY_CHO_H-HCB_A-AA',
            'C3_TPBZ_NH2_H-C2_TIDA_CHO_H-HCB_A-AA',
            'C3_TPBZ_NH2_H-C2_BENZ_CHO_OMe-HCB_A-AA',
            'C3_TPTZ_NH2_H-C3_TPTZ_CHO_H-HCB-AA',
            'C3_TPTZ_NH2_H-C3_BENZ_CHO_H-HCB-AA',
            'C3_TPBZ_NH2_H-C2_BENZ_CHO_CH3-HCB_A-AA',
            'C3_TPTZ_NH2_H-C2_TTPH_CHO_H-HCB_A-AA',
            'C3_STAR_(OH)2_H-C2_BENZ_B(OH)2_H-HCB_A-AA',
            'C3_STAR_(OH)2_H-C2_PYRN_B(OH)2_H-HCB_A-AA',
            'C3_STAR_(OH)2_H-C2_BPNY_B(OH)2_H-HCB_A-AA',
            'C3_DBA1_(OH)2_H-C2_BENZ_B(OH)2_H-HCB_A-AA',
            'C3_DBA2_(OH)2_H-C2_BENZ_B(OH)2_H-HCB_A-AA',
            'C3_BENZ_CHO_OH-C2_BENZ_NH2_H-HCB_A-AA',
            'C3_BENZ_CHO_OH-C2_BPNY_NH2_H-HCB_A-AA',
            'C3_BENZ_CHO_OH-C2_BPNY_NH2_OMe-HCB_A-AA',
            'C3_BENZ_CHO_OH-C3_TPBZ_NH2_H-HCB-AA',
            'C3_TPNY_(OH)2_H-C2_BENZ_B(OH)2_H-HCB_A-AA',
            'C3_BENZ_B(OH)2_H-C3_TPNY_(OH)2_H-HCB-AA',
            'C3_TPBZ_B(OH)2_H-C3_TPNY_(OH)2_H-HCB-AA',
            'C3_TPNY_(OH)2_H-C2_BPNY_B(OH)2_H-HCB_A-AA',
            'C3_BENZ_CHO_H-C2_BENZ_CONHNH2_OEt-HCB_A-AA',
            'C3_TPBZ_CHO_H-C2_BENZ_CONHNH2_OEt-HCB_A-AA',
            'C3_TPBZ_NH2_H-C2_BENZ_CHO_OH-HCB_A-AA',
            'C3_TPBZ_NH2_H-C2_BENZ_CHO_OMe-HCB_A-AA',
            'C3_TPAM_NH2_H-C2_INTO_O_H-HCB_A-AA',
            'C3_TPBZ_NH2_H-C2_INTO_O_H-HCB_A-AA',
            'C3_TBBZ_NH2_H-C2_INTO_O_H-HCB_A-AA',
            'C3_TPBZ_NH2_H-C2_PTCD_O_H-HCB_A-AA',
            'C3_TPBZ_NH2_H-C2_INTO_O_H-HCB_A-AA',
            'C3_TPAM_NH2_H-C2_INTO_O_H-HCB_A-AA',
            'C3_TPTZ_NH2_H-C2_INTO_O_H-HCB_A-AA',
            'C3_BENZ_CHO_OH-C2_BENZ_NH2_H-HCB_A-AA',
            'C3_BENZ_CHO_OH-C2_ANTR_NH2_O-HCB_A-AA']

for cof in COF_list:
    COF.build(cof, supercell=[1,1,1], save_format=['cif'])

#COF.build('C3_BENZ_CHO_OH-C2_HDZ_NH2-HCB_A-AA', supercell=[1,1,1])
#COF.build('C3_BENZ_CHO_H-C2_HDZ_NH2-HCB_A-AA', supercell=[1,1,1])
#COF.build('C3_BENZ_CHO_NH2-C2_BENZ_NH2-HCB_A-AA', supercell=[1,1,1])


#COF.clean()

#COF.create_all_C2(nucleos=['BENZ', 'TPNY', 'ANTR', 'BPNY', 'BPYB', 'DHSI', 'DPBY', 'DPEL', 'DPEY', 'NAPT', 'PYRN'], conectores=['NH2'], radicais=['CH3', 'CHO', 'CN', 'COOH', 'F', 'H', 'NC2H', 'NH2', 'NO2', 'OEt', 'OH', 'OMe', 'tBut'])

#COF.create_all_C2(nucleos=['HDZ', 'BENZ'], conectores=['NH2'], radicais=['H', 'NH2', 'NO2', 'OH', 'OMe'])


#COF.create_all_C3(nucleos=['BENZ', 'TPBZ'], conectores=['CHO'], radicais=['H'])

#COF.build_all_available_COFs()




