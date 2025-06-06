#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

import os, sys
import numpy as np
import pandas as pd
from datetime import timedelta
from NANUK_lkf_tools  import *
import pickle
import calendar

#----  nanuk12_driver_LKF_detect.py --------------------------
#
# Driver that loops through a series of files (dates) and that
# calls the funtion NANUK_lkf_detect.
#
#------------------------------------------------------------

myhost = os.uname()[1]

print('  *** Host: '+myhost)

#----- INPUT -----
cregflag = 1 # 1: output includes vorticity, 2: no vorticity
CONF = 'NANUK12'

EXP='CPL00'

if   myhost=='frazilo':
    main_dir = '/data/laurent/'+CONF
elif myhost in ['luitel']:
    main_dir = '/data/gcm_setup/'+CONF
else:
    print(' UNKNOWN HOST: '+myhost+' !\n')
    

main_dir_grid = main_dir+'/'+CONF+'.L31-I'

main_dir_exp  = main_dir+'/'+CONF+'-'+EXP+'_BBM'

store_main_dirTP = '/data/gcm_setup/tmp/'+CONF

kvalue = 7 # value for kernel
produce_plot = True
pack_ice_mask = False
SDATE = '19970202'
EDATE = '19970206'
FREQ = '3h'
DTGAP = '48h'
suffix = 'icemod'

#----- check kernel value ----------------
# kvalue should be odd: see Nils' email (4 nov 2022)

if kvalue % 2 == 0:
    print("kernel value should be an odd integer")
    exit()
else:
    print("kernel value = ")
    print(kvalue)

#----- define paths and file name --------

grid_path=os.path.join(main_dir_grid+'/coordinates_'+CONF+'.nc')


if (pack_ice_mask):
    store_main_dir=store_main_dirTP+'/LKF_diag_pack'
else:
    store_main_dir=store_main_dirTP+'/LKF_diag'

store_path=os.path.join(store_main_dir+'/'+EXP+'/detectedLKFs/')

list_dates=list(pd.date_range(SDATE,EDATE, freq=DTGAP))


print( "list_dates = ", list_dates )


for i in range(len(list_dates)) :

    
    #date0 = (list_dates[i] + timedelta(days=-0)).strftime('%Y%m%d%H')
    date0 = (list_dates[i] + timedelta(days=-0)).strftime('%Y%m%d')
    print(' * date0 = ',date0)
    data_path = main_dir_exp+'/'+CONF+'-'+EXP+'_'+FREQ+'_'+date0+'_'+suffix+'.nc4'
    print(' * data_path = ',data_path)
    #data_path=os.path.join(main_dir+'/'+EXP+'/hourly/'+date0+suffix+'.nc')
    fileout=date0+'_'+EXP
    #print(fileout)
    print('\n *** Going for detection of '+fileout+':')
    NANUK_lkf_detect(date0, CONF, cregflag, grid_path, data_path, store_path, fileout, kvalue, produce_plot, pack_ice_mask)
    print('\n '+fileout+' detecttion done!\n\n')


    
print('Detection done for experiment:')
print(CONF+'-'+EXP)


exit(0)

#for i in range(len(list_dates)) :
#    date0 = (list_dates[i] + timedelta(days=-0)).strftime('%Y%m%d%H')
#    data_path=os.path.join(main_dir+'/'+EXP+'/hourly/'+date0+suffix+'.nc')
#    fileout=date0 + '_' + EXP
#    print(fileout)
#    NANUK_lkf_detect(date0, CONF, cregflag, grid_path, data_path, store_path, fileout, kvalue, produce_plot, pack_ice_mask)

data_path = main_dir_exp+'/'+CONF+'-'+EXP+'_'+FREQ+'_'+SDATE+'_'+EDATE+'_'+suffix+'.nc4'

list_dates=list(pd.date_range(SDATE,SDATE, freq=FREQ))
date0 = (list_dates[0] + timedelta(days=-0)).strftime('%Y%m%d%H')

print('*** date0 =', date0)

fileout=date0+'_'+CONF+'_'+EXP

