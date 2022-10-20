import os,sys
import numpy as np
import pandas as pd
from datetime import timedelta
from CREG_lkf_tools  import *
import pickle
import calendar

#----  CREG_driver_LKF_analysis -----------------------------
#
# Driver that loops through a series of files (dates) and that 
# calls the funtion CREG_lkf_analysis that calculates the LKF
# half widths.
#
#------------------------------------------------------------

#----- INPUT -----
#ni = 528 ; creg025
#nj = 735 ;
#ni = 1580 ; creg12
#nj = 2198 ;
creggrid='creg025' # creg025 or creg12
EXP='run6f'
main_dir='/home/jfl001/data/Lemieux2022/LKF_diag'
main_dirnc='/home/jfl001/data/runsLemieux_et_al_2022/'
dir_util='/home/jfl001/Lemieux2022/UTIL'
dsearch=5 # +- dsearch cells around one LKF cell (dist is capped if searching too far!!!)
frac=0.5 # half width is defined as eps_tot < frac*LKFepsmax 
mindist=150.0 # LKF point is analysed if dist from land > mindist (km)

FREQ='24H'
SDATE='20050329'
EDATE='20050329'
suffix='_000'
#-----------------------------------------

list_dates=list(pd.date_range(SDATE,EDATE, freq=FREQ))

for i in range(len(list_dates)) :
    date0 = (list_dates[i] + timedelta(days=-0)).strftime('%Y%m%d%H')
    date0ext=date0 + '_000'
    filein='lkf_' + date0ext + '_' + EXP + '_001.npy'
    fileout='lkf_' + date0ext + '_' + EXP + '_a.npy' # a for analysed
    tpdir=date0ext + '_' + EXP
    path_filein=os.path.join(main_dir+'/'+creggrid+'/'+EXP+'/'+tpdir+'/'+filein)
    path_fileout=os.path.join(main_dir+'/'+creggrid+'/'+EXP+'/'+tpdir+'/'+fileout)
    data_path=os.path.join(main_dirnc+creggrid+'/'+EXP+'/netcdf/'+date0ext+'.nc')
    path_filedist=os.path.join(dir_util +'/dist_'+creggrid+'.pkl')

    CREG_lkf_analysis(date0,creggrid,path_filedist,path_filein,path_fileout,data_path,dsearch,frac,mindist)
