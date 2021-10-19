# -*- coding: utf-8 -*-

import pycofbuilder as COF

COF.clean()

COF.create_all_C2(nucleos=['BENZ', 'TPNY', 'ANTR', 'BPNY', 'BPYB', 'DHSI', 'DPBY', 'DPEL', 'DPEY', 'NAPT', 'PYRN'], conectores=['NH2'], radicais=['CH3', 'CHO', 'CN', 'COOH', 'F', 'H', 'NC2H', 'NH2', 'NO2', 'OEt', 'OH', 'OMe', 'tBut'])

#COF.create_all_C2(nucleos=['HDZ', 'BENZ'], conectores=['NH2'], radicais=['H', 'NH2', 'NO2', 'OH', 'OMe'])


COF.create_all_C3(nucleos=['BENZ', 'TPBZ'], conectores=['CHO'], radicais=['H'])

COF.build_all_available_COFs()




