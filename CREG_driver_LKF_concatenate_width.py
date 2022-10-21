import os,sys
import numpy as np
import pandas as pd
from datetime import timedelta
from CREG_lkf_tools import CREG_lkf_concatenate_width
import pickle
import calendar

#----  CREG_driver_concatenate_width ------------------------
#
# Driver that loops through a series of files (dates) and that 
# calls the funtion CREG_lkf_concatenate_width. That function
# concatenate the half widths for one LKF file (one date) and
# the concatenation for multiple files (mutiple dates) is done
# here in the driver.
#
#------------------------------------------------------------

#----- INPUT -----
creggrid='creg025' # creg025 or creg12
EXP='run8fb'
main_dir='/home/jfl001/data/Lemieux2022/LKF_diag'

FREQ='24H'
SDATE='20050201'
EDATE='20050228'
suffix='_000'
#-----------------------------------------

widthdir=os.path.join(main_dir+'/'+creggrid+'/'+EXP+'/WIDTH/')

if not os.path.isdir(widthdir):
    os.makedirs(widthdir)

fileout1='hwidth1_lkf_'+SDATE+'_'+EDATE+'.npy'
path_fileout1=os.path.join(widthdir+fileout1)
fileout2='hwidth2_lkf_'+SDATE+'_'+EDATE+'.npy'
path_fileout2=os.path.join(widthdir+fileout2)

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
    hwidth1.extend(tpvect)

    tpvect=[]
    path_filein=os.path.join(main_dir+'/'+creggrid+'/'+EXP+'/'+tpdir+'/'+filein)
    tpvect=CREG_lkf_concatenate_width (date0, path_filein, hwidth=2)
    hwidth2.extend(tpvect)

np.save(path_fileout1,hwidth1,allow_pickle=True)
np.save(path_fileout2,hwidth2,allow_pickle=True)
