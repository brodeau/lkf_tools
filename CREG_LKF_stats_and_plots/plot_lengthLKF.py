import os,sys
import numpy as np
import pandas as pd
#import matplotlib as plt
import matplotlib.pyplot as plt
from datetime import timedelta
#import pickle
import calendar

EXP1='exp1'
label1='exp1'
EXP2='control'
label2='control'
#EXP3='eg2p25_ef1p5'
#label3='e_g=2.25'
main_dir='/home/jfl001/data/LKF_rips_analysis'
zdir='Length'
SDATE='20200101'
EDATE='20201230'
FREQ='168H'
addlabel='length'

#-----------------------------------------

list_dates=list(pd.date_range(SDATE,EDATE, freq=FREQ))

mlength1=[] # mean length
tlength1=[] # total length
mlength2=[] # mean length
tlength2=[] # total length
#mlength3=[] # mean length
#tlength3=[] # total length
datesl=[] # lits of dates
for i in range(len(list_dates)) :
    date0 = (list_dates[i] + timedelta(days=-0)).strftime('%Y%m%d%H')
    filein1=date0+'_'+addlabel+'_'+EXP1+'.py'
    path_filein1=os.path.join(main_dir+'/'+EXP1+'/'+zdir+'/'+filein1)
    df1 = pd.read_csv(path_filein1)
    mlength1.append(df1['length'].mean())
    tlength1.append(df1['length'].sum())
    datesl.append(list_dates[i])

    filein2=date0+'_'+addlabel+'_'+EXP2+'.py'
    path_filein2=os.path.join(main_dir+'/'+EXP2+'/'+zdir+'/'+filein2)
    df2 = pd.read_csv(path_filein2)
    mlength2.append(df2['length'].mean())
    tlength2.append(df2['length'].sum())

#    filein3=date0+'_'+addlabel+'_'+EXP3+'.py'
#    path_filein3=os.path.join(main_dir+'/'+EXP3+'/'+zdir+'/'+filein3)
#    df3 = pd.read_csv(path_filein3)
#    mlength3.append(df3['length'].mean())
#    tlength3.append(df3['length'].sum())


#--- create the panda dataframes for plotting

dfpm = pd.DataFrame(datesl, columns=['date'])
dfpm.insert(1, 'mean_length1', mlength1)
dfpm.insert(2, 'mean_length2', mlength2)

dfpt = pd.DataFrame(datesl, columns=['date'])
dfpt.insert(1, 'total_length1', tlength1)
dfpt.insert(2, 'total_length2', tlength2)


ax = dfpm.plot(x = 'date', y =['mean_length1', 'mean_length2'], color =['dodgerblue', 'orange'])
ax.legend([EXP1, EXP2])
ax.set_xlabel("", fontsize='2')
ax.set_ylabel("Mean length of LKFs (km)", fontsize=14)
plt.savefig('Mean_length_LKFs_2020.png')

ax = dfpt.plot(x = 'date', y =['total_length1', 'total_length2'], color =['dodgerblue', 'orange'])
ax.legend([EXP1, EXP2])
ax.set_xlabel("", fontsize='2')
ax.set_ylabel("Total length of LKFs (km)", fontsize=14)
plt.savefig('Total_length_LKFs_2020.png')

#plt.show()
