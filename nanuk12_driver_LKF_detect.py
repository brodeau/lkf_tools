#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
##################################################################

import os,sys
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

#----- INPUT -----
cregflag = 1 # 1: output includes vorticity, 2: no vorticity
CONF = 'NANUK12'

EXP='CPL00'

main_dir = '/data/gcm_setup/'+CONF

main_dir_grid = main_dir+'/'+CONF+'.L31-I'

main_dir_exp  = main_dir+'/'+CONF+'-'+EXP+'_BBM'

store_main_dirTP = '/data/gcm_setup/tmp/'+CONF

kvalue = 7 # value for kernel
produce_plot = True
pack_ice_mask = False
SDATE = '19970201'
EDATE = '19970202'
FREQ = '3h'
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

#list_dates=list(pd.date_range(SDATE,EDATE, freq=FREQ))

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

NANUK_lkf_detect(date0, CONF, cregflag, grid_path, data_path, store_path, fileout, kvalue, produce_plot, pack_ice_mask)

print('Detection done for experiment:')
print(CONF+'-'+EXP)
