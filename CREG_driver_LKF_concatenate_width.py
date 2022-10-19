import os,sys
import numpy as np
import pandas as pd
from datetime import timedelta
from CREG_lkf_tools import CREG_lkf_width
import pickle
import calendar

#----- INPUT -----
creggrid='creg025' # creg025 or creg12
EXP='run6f'
main_dir='/home/jfl001/data/Lemieux2022/LKF_diag'

FREQ='24H'
SDATE='20050329'
EDATE='20050329'
suffix='_000'
#-----------------------------------------

fileout1='hwidth1_lkf_'+SDATE+'_'+EDATE+'.npy'
path_fileout1=os.path.join(main_dir+'/'+creggrid+'/'+EXP+'/WIDTH/'+fileout1)
fileout2='hwidth2_lkf_'+SDATE+'_'+EDATE+'.npy'
path_fileout2=os.path.join(main_dir+'/'+creggrid+'/'+EXP+'/WIDTH/'+fileout2)

list_dates=list(pd.date_range(SDATE,EDATE, freq=FREQ))

hwidth1=[]
hwidth2=[]
tpvect=[]
for i in range(len(list_dates)):
    date0 = (list_dates[i] + timedelta(days=-0)).strftime('%Y%m%d%H')
    date0ext=date0 + '_000'
    filein='lkf_' + date0ext + '_' + EXP + '_a.npy' # a for analysed
    tpdir=date0ext + '_' + EXP
    
    tpvect=[]
    path_filein=os.path.join(main_dir+'/'+creggrid+'/'+EXP+'/'+tpdir+'/'+filein)
    tpvect=CREG_lkf_concatenate_width (date0, path_filein, hwidth=1)
    hwidth1.append(tpvect)

    tpvect=[]
    path_filein=os.path.join(main_dir+'/'+creggrid+'/'+EXP+'/'+tpdir+'/'+filein)
    tpvect=CREG_lkf_concatenate_width (date0, path_filein, hwidth=2)
    hwidth2.append(tpvect)

np.save(path_fileout1,hwidth1,allow_pickle=True)
np.save(path_fileout2,hwidth2,allow_pickle=True)
