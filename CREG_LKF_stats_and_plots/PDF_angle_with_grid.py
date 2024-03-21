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
SDATE='20060101'
EDATE='20060531'
FREQ='24H'
addlabel='anggrid'

#-----------------------------------------

list_dates=list(pd.date_range(SDATE,EDATE, freq=FREQ))

for i in range(len(list_dates)) :
    date0 = (list_dates[i] + timedelta(days=-0)).strftime('%Y%m%d%H')
    filein=date0+'_'+addlabel+'_'+EXP+'.py'
    path_filein=os.path.join(main_dir+'/'+EXP+'/'+zdir+'/'+filein)
    if i==0: # date0 is SDATE
        df1 = pd.read_csv(path_filein)
    else:
        df2 = pd.read_csv(path_filein)
        df1 = pd.concat([df1, df2])

#--- define bins 0-90 deg ---
nbbins=90
delta=1.0
mybins=np.zeros(nbbins+1)
binc=np.zeros(nbbins)
for b in range(nbbins+1):
    mybins[b]=b*delta

for b in range(nbbins):
    binc[b]=0.5*(mybins[b]+mybins[b+1])

#--- define bins 0-45 deg ---
nbbins=45
delta=1.0
myotherbins=np.zeros(nbbins+1)
otherbinc=np.zeros(nbbins)
for b in range(nbbins+1):
    myotherbins[b]=b*delta

for b in range(nbbins):
    otherbinc[b]=0.5*(mybins[b]+mybins[b+1])

#--- calc mean values ---
mean_x_angle=df1['x_angle'].mean()
print('mean angle with x axis', mean_x_angle)
mean_y_angle=df1['y_angle'].mean()
print('mean angle with y axis', mean_y_angle)
mean_min_angle=df1['min_angle'].mean()
print('mean min angle with x or y axis', mean_min_angle)

plt.figure(1)
counts, bins, bars = plt.hist(df1['x_angle'], bins=mybins, density=True, color = "dodgerblue", ec="dodgerblue")
plt.xlabel('angle', fontsize=14)
plt.ylabel('PDF (angle with x axis)', fontsize=14)
plt.figure(2)
counts, bins, bars = plt.hist(df1['y_angle'], bins=mybins, density=True, color = "dodgerblue", ec="dodgerblue")
plt.xlabel('angle', fontsize=14)
plt.ylabel('PDF (angle with y axis)', fontsize=14)
plt.figure(3)
counts, bins, bars = plt.hist(df1['min_angle'], bins=myotherbins, density=True, color = "dodgerblue", ec="dodgerblue")
plt.xlabel('Minimum angle with x or y axis', fontsize=14)
plt.ylabel('PDF', fontsize=14)

plt.show()
