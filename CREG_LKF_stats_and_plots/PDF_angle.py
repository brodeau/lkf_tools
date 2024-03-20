import os,sys
import numpy as np
import pandas as pd
#import matplotlib as plt
import matplotlib.pyplot as plt
from datetime import timedelta
#import pickle
import calendar

EXP='eg1p5_ef1p5'
main_dir='/home/jfl001/data/Lemieux_et_al_plast_pot/LKF_diag'
zdir='Angle_grid'
SDATE='20060301'
EDATE='20060331'
FREQ='24H'
addlabel='anggrid'

#-----------------------------------------

list_dates=list(pd.date_range(SDATE,EDATE, freq=FREQ))

for i in range(len(list_dates)) :
    date0 = (list_dates[i] + timedelta(days=-0)).strftime('%Y%m%d%H')
    filein=date0+'_'+addlabel+'_'+EXP+'.py'
    path_filein=os.path.join(main_dir+'/'+EXP+'/'+zdir+'/'+filein)
    print(path_filein)
    if i==0: # date0 is SDATE
        df1 = pd.read_csv(path_filein)
    else:
        df2 = pd.read_csv(path_filein)
        df1 = pd.concat([df1, df2])

#pangle=df.x_angle
#arr=pangle.to_numpy()
#print(pangle)
#print(arr.shape)

#--- define bins ---
nbbins=90
delta=1.0
mybins=np.zeros(nbbins+1)
binc=np.zeros(nbbins)
for b in range(nbbins+1):
    mybins[b]=b*delta

for b in range(nbbins):
    binc[b]=0.5*(mybins[b]+mybins[b+1])


plt.figure(1)
counts, bins, bars = plt.hist(df1['x_angle'], bins=mybins, color = "dodgerblue", ec="dodgerblue")
plt.xlabel('angle', fontsize=14)
plt.figure(2)
counts, bins, bars = plt.hist(df1['y_angle'], bins=mybins, color = "dodgerblue", ec="dodgerblue")
plt.xlabel('angle', fontsize=14)
plt.figure(3)
counts, bins, bars = plt.hist(df2['x_angle'], bins=mybins, color = "dodgerblue", ec="dodgerblue")
plt.xlabel('angle', fontsize=14)
plt.figure(4)
counts, bins, bars = plt.hist(df2['y_angle'], bins=mybins, color = "dodgerblue", ec="dodgerblue")
plt.xlabel('angle', fontsize=14)
#plt.bar(binc,counts, width = 0.2)
#plt.ylabel('PDF', fontsize=14)
#plt.xlim([0, 90]) 
#plt.ylim([0, 0.2]) 
#plt.savefig(path_fileout)
plt.show()
