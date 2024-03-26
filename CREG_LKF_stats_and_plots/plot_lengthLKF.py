import os,sys
import numpy as np
import pandas as pd
#import matplotlib as plt
import matplotlib.pyplot as plt
from datetime import timedelta
#import pickle
import calendar

EXP1='eg1p0_ef1p5'
label1='e_g=1.0'
EXP2='eg1p5_ef1p5'
label2='e_g=1.5'
EXP3='eg2p25_ef1p5'
label3='e_g=2.25'
main_dir='/home/jfl001/data/Lemieux_et_al_plast_pot/LKF_diag'
zdir='Length'
SDATE='20050101'
EDATE='20050531'
FREQ='24H'
addlabel='length'

#-----------------------------------------

list_dates=list(pd.date_range(SDATE,EDATE, freq=FREQ))

mlength1=[] # mean length
tlength1=[] # total length
mlength2=[] # mean length
tlength2=[] # total length
mlength3=[] # mean length
tlength3=[] # total length
for i in range(len(list_dates)) :
    date0 = (list_dates[i] + timedelta(days=-0)).strftime('%Y%m%d%H')
    filein1=date0+'_'+addlabel+'_'+EXP1+'.py'
    path_filein1=os.path.join(main_dir+'/'+EXP1+'/'+zdir+'/'+filein1)
    df1 = pd.read_csv(path_filein1)
    mlength1.append(df1['length'].mean())
    tlength1.append(df1['length'].sum())

    filein2=date0+'_'+addlabel+'_'+EXP2+'.py'
    path_filein2=os.path.join(main_dir+'/'+EXP2+'/'+zdir+'/'+filein2)
    df2 = pd.read_csv(path_filein2)
    mlength2.append(df2['length'].mean())
    tlength2.append(df2['length'].sum())

    filein3=date0+'_'+addlabel+'_'+EXP3+'.py'
    path_filein3=os.path.join(main_dir+'/'+EXP3+'/'+zdir+'/'+filein3)
    df3 = pd.read_csv(path_filein3)
    mlength3.append(df3['length'].mean())
    tlength3.append(df3['length'].sum())


plt.figure(1)
plt.plot(tlength1)
plt.plot(tlength2)
plt.plot(tlength3)
plt.legend([label1, label2, label3], loc ="lower right")


#plt.xlabel('angle', fontsize=14)
plt.ylabel('Total length of LKFs (km)', fontsize=14)
plt.xlabel('days', fontsize=14)
#plt.figure(2)

plt.figure(2)
plt.plot(mlength1)
plt.plot(mlength2)
plt.plot(mlength3)
plt.legend([label1, label2, label3], loc ="upper right")
plt.ylabel('Mean length of LKFs (km)', fontsize=14)
plt.xlabel('days', fontsize=14)

plt.show()
