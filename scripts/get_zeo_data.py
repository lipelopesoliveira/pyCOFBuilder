# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:31:19 2020
@author: lipel
"""
import os
import numpy as np

def get_zeo_results(path, file_name):

    zeo_dict = {}

    if os.path.exists(os.path.join(path, file_name + '.chan')):

        chan_file = open(os.path.join(path, file_name + '.chan'), 'r').readlines()
        zeo_dict['number_of_channels'] = int(chan_file[0].split()[1])
        zeo_dict['channel_dimentionality'] = int(chan_file[0].split()[-1])
        zeo_dict['channel_size'] = [float(i) for i in chan_file[1].split()[2:]]
    
    if os.path.exists(os.path.join(path, file_name + '.vol')):
        vol_file = open(os.path.join(path, file_name)+ '.vol', 'r').readlines()
        zeo_dict['density'] = float(vol_file[0].split()[5])
        zeo_dict['acessible_volume_fraction'] = float(vol_file[0].split()[9])
        zeo_dict['pore_volume_cm3/g'] = float(vol_file[0].split()[11])
    
    if os.path.exists(os.path.join(path, file_name + '.sa')):
        sa_file = open(os.path.join(path, file_name) + '.sa', 'r').readlines()
        zeo_dict['specific_area'] = {'A2': float(sa_file[0].split()[7]),
                                     'm2/cm3': float(sa_file[0].split()[9]),
                                     'm2/g': float(sa_file[0].split()[11])} 

    if os.path.exists(os.path.join(path, file_name + '.psd_histo')):
        temp_file = open(os.path.join(path, file_name) + '.psd_histo', 'r').readlines()

        histo_data = [i.split() for i in temp_file[11:]]
        bins, count, cumm, deriv = np.transpose(histo_data)

        bins = [float(i) for i in bins]
        count = [float(i) for i in count]
        cumm = [float(i) for i in cumm]
        deriv = [float(i) for i in deriv]

        zeo_dict['psd'] = {'bins': bins, 'count':count}

    return zeo_dict
	