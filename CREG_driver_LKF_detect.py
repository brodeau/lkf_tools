import os,sys
import numpy as np
import pandas as pd
from datetime import timedelta
from CREG_lkf_tools  import *
import pickle
import calendar

#----  CREG_driver_LKF_detect -------------------------------
#
# Driver that loops through a series of files (dates) and that 
# calls the funtion CREG_lkf_detect.
#
#------------------------------------------------------------

#----- INPUT -----
#ni = 528 ; creg025
#nj = 735 ;
#ni = 1580 ; creg12
#nj = 2198 ;
cregflag=1 # used to customize code for CREG applications.
creggrid='creg025' # creg025 or creg12
EXP='run12fb'
main_dir='/home/jfl001/data/runsLemieux_et_al_2022/'
main_dir_grid='/home/socn000/data/eccc-ppp5/env_rhel-8-icelake-64/datafiles/constants/oce/repository/master/CONCEPTS/'
store_main_dir='/home/jfl001/data/Lemieux2022/LKF_diag'
kvalue=7 # value for kernel
produce_plot=False
FREQ='24H'
SDATE='20050217'
EDATE='20050217'
suffix='_000'

#----- define paths and file name --------

store_path=os.path.join(store_main_dir+'/'+creggrid+'/'+EXP+'/')

if (creggrid == 'creg025'):
    grid_path=os.path.join(main_dir_grid+'/creg025pe/grid/coordinates_CREG025_LIM.nc')
elif (creggrid == 'creg12'):
    grid_path=os.path.join(main_dir_grid+'/creg012pe/grid/coordinates_CREG12_ext.nc')
else:
    print ("Wrong choice of grid")


list_dates=list(pd.date_range(SDATE,EDATE, freq=FREQ))

for i in range(len(list_dates)) :
    date0 = (list_dates[i] + timedelta(days=-0)).strftime('%Y%m%d%H')
    data_path=os.path.join(main_dir+creggrid+'/'+EXP+'/netcdf/'+date0+suffix+'.nc')
    fileout=date0 + suffix + '_' + EXP
    CREG_lkf_detect(date0, creggrid, grid_path, data_path, store_path, fileout, kvalue, produce_plot)

print('Detection done for experiment:')
print(EXP)
