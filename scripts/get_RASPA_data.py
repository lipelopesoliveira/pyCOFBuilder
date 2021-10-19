# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020
@author: lipel
"""
import os
import numpy as np
import traceback


def read_raspa_single_ads(path, gas, temp):
    raspa_out_path = os.path.join(path, temp, gas, 'Output', 'System_0')
    if os.path.exists(raspa_out_path):
        out_files = os.listdir(raspa_out_path)
        ads_dict = {'gas': gas,
                    'temperature': temp,
                    'pressure': [],
                    'mol/unit_cell':[],
                    'sd_mol/unit_cell':[],
                    'mol/kg':[],
                    'sd_mol/kg':[],
                    'mg/g':[],
                    'sd_mg/g':[],
                    'cm3/g':[],
                    'sd_cm3/g':[],
                    'cm3/cm3':[],
                    'sd_cm3/cm3':[],
                    'enthalpy':[],
                    'sd_enthalpy':[],
                    'heat_capacity_J/mol/K':[],
                    'sd_heat_capacity_J/mol/K':[],
                    'heat_capacity_cal/mol/K':[],
                    'sd_heat_capacity_cal/mol/K':[]
                    }

        for file in out_files:
            tmp = open(os.path.join(raspa_out_path, file), 'r').readlines()
            if 'Simulation finished' in tmp[-3]:
                for line in tmp:
                    if 'External temperature' in line: 
                        ads_dict['temperature'] = float(line.split()[-2])
                    if 'External Pressure' in line: 
                        ads_dict['pressure'] += [float(line.split()[-2])]
                    if '[J/mol/K] +/-  ' in line:
                        ads_dict['heat_capacity_J/mol/K'] += [float(line.split()[1])]
                        ads_dict['sd_heat_capacity_J/mol/K'] += [float(line.split()[-2])]
                    if '[cal/mol/K] +/- ' in line:
                        ads_dict['heat_capacity_cal/mol/K'] += [float(line.split()[1])]
                        ads_dict['sd_heat_capacity_cal/mol/K'] += [float(line.split()[-2])]
                    if '[KJ/MOL]' in line:
                        ads_dict['enthalpy'] += [float(line.split()[0])]
                        ads_dict['sd_enthalpy'] += [float(line.split()[-2])]
                    if 'Average loading absolute [molecules/unit cell]' in line:
                        ads_dict['mol/unit_cell'] += [float(line.split()[-4])]
                        ads_dict['sd_mol/unit_cell'] += [float(line.split()[-2])]
                    if 'Average loading absolute [mol/kg framework]' in line:
                        ads_dict['mol/kg'] += [float(line.split()[-4])]
                        ads_dict['sd_mol/kg'] += [float(line.split()[-2])]
                    if 'Average loading absolute [milligram/gram framework]' in line:
                        ads_dict['mg/g'] += [float(line.split()[-4])]
                        ads_dict['sd_mg/g'] += [float(line.split()[-2])]
                    if 'Average loading absolute [cm^3 (STP)/gr framework]' in line:
                        ads_dict['cm3/g'] += [float(line.split()[-4])]
                        ads_dict['sd_cm3/g'] += [float(line.split()[-2])]
                    if 'Average loading absolute [cm^3 (STP)/cm^3 framework]' in line:
                        ads_dict['cm3/cm3'] += [float(line.split()[-4])]
                        ads_dict['sd_cm3/cm3'] += [float(line.split()[-2])]

        return ads_dict
    else:
        print(raspa_out_path, 'do not exist! Skipping this structure.')
                
