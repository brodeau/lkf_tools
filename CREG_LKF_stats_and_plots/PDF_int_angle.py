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
zdir='Int_Angle'
SDATE='20050101'
EDATE='20050111'
FREQ='24H'
addlabel='intpairs'

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

#--- calc mean values ---
#mean_x_angle=df1['x_angle'].mean()
#print('mean angle with x axis', mean_x_angle)

nrow=len(df1)
#print('nrow=',nrow)

#print(df1['clean_int'])
#for n in range(nrow):
#    print(df1['clean_int'][n])
#    aa=1

for row in df1.iterrows():
    print(row['conj_pair'])


#plt.figure(1)
#counts, bins, bars = plt.hist(df1['x_angle'], bins=mybins, density=True, color = "dodgerblue", ec="dodgerblue")
#plt.xlabel('angle', fontsize=14)
#plt.ylabel('PDF (angle with x axis)', fontsize=14)

#plt.show()
