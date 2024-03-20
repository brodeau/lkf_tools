import os,sys
import numpy as np
import pandas as pd
from datetime import timedelta
from CREG_lkf_tools  import *
import pickle
import calendar

#----  CREG_driver_LKF_get_angles ------------------------------
#
# finds pairs of intersecting LKFs, calc angles and identify 
# conjugate fault lines 
# 
# note: there is no condition applied here for distance to 
#       land. This could be applied later for plotting. 
#
#------------------------------------------------------------

#----- INPUT -----
#ni = 528 ; creg025
#nj = 735 ;
#ni = 1580 ; creg12
#nj = 2198 ;
creggrid='creg12' # creg025 or creg12
EXP='eg1p5_ef1p5'
main_dir='/home/jfl001/data/Lemieux_et_al_plast_pot/LKF_diag'
store_main_dir='/home/jfl001/data/Lemieux_et_al_plast_pot/LKF_diag'
SDATE='20060301'
EDATE='20060301'
FREQ='24H'
suffix='0000_iceh_inst'

#-----------------------------------------

store_path=os.path.join(main_dir+'/'+EXP+'/Angle_grid')
if not os.path.isdir(store_path):
    os.mkdir(store_path)

list_dates=list(pd.date_range(SDATE,EDATE, freq=FREQ))

for i in range(len(list_dates)) :
    date0 = (list_dates[i] + timedelta(days=-0)).strftime('%Y%m%d%H')
    filein='lkf_' + date0 + '_' + EXP + '_001.npy'
    tpdir=date0 + '_' + EXP
    path_filein=os.path.join(main_dir+'/'+EXP+'/'+tpdir+'/'+filein)
    print(path_filein)
    fileout=os.path.join(store_path + '/' + date0 + '_anggrid_' + EXP + '.py')
    print(fileout)
    CREG_lkf_angles_with_grid(date0,creggrid,path_filein,fileout)

print('CREG_driver_LKF_angles_with_grid is done')
print(EXP)
