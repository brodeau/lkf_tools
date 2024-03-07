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

SDATE='20061221'
EDATE='20061221'
FREQ='24H'

#-----------------------------------------

store_path=os.path.join(main_dir+'/'+EXP+'/')

list_dates=list(pd.date_range(SDATE,EDATE, freq=FREQ))

for i in range(len(list_dates)) :
    date0 = (list_dates[i] + timedelta(days=-0)).strftime('%Y%m%d%H')
    filein='lkf_' + date0 + '_' + EXP + '_001.npy'
    tpdir=date0 + '_' + EXP
    path_filein=os.path.join(main_dir+'/'+EXP+'/'+tpdir+'/'+filein)
    print(path_filein)
    #lkfs = np.load(path_filein,allow_pickle=True)
    #fileout=date0 + '_pairs_' + EXP
    CREG_lkf_pairs_and_angles(date0,creggrid,path_filein)

print('CREG_driver_LKF_pairs_and_angles is done')
print(EXP)
