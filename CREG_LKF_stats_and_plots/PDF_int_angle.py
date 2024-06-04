import os,sys
import numpy as np
import pandas as pd
#import matplotlib as plt
import matplotlib.pyplot as plt
from datetime import timedelta
#import pickle
import calendar

#EXP='run_eg1p0_ef1p5'
#EXP='run_eg1p5_ef1p5'
#EXP='run_eg2p25_ef1p5'
#EXP='run_eg1p16_ef1p75'
#EXP='run_eg1p75_ef1p75'
EXP='run_eg2p63_ef1p75'
#EXP='run_eg1p33_ef2p0'
#EXP='run_eg2p0_ef2p0'
#EXP='run_eg3p0_ef2p0'


main_dir='/home/jfl001/data/Lemieux_et_al_plast_pot/LKF_diag'
zdir='Int_Angle'
year='2005'
SDATE=year+'0101'
EDATE=year+'0531'
FREQ='24H'
addlabel='intpairs'
dlabel='10'

nbmin=10

#-----------------------------------------

list_dates=list(pd.date_range(SDATE,EDATE, freq=FREQ))

for i in range(len(list_dates)) :
    date0 = (list_dates[i] + timedelta(days=-0)).strftime('%Y%m%d%H')
    filein=date0+'_'+addlabel+'_'+EXP+'_delta'+dlabel+'.py'
    path_filein=os.path.join(main_dir+'/'+EXP+'/'+zdir+'/'+filein)
    if i==0: # date0 is SDATE
        df1 = pd.read_csv(path_filein)
    else:
        df2 = pd.read_csv(path_filein)
        df1 = pd.concat([df1, df2])

#--- define bins 0-90 deg ---
nbbins=36
delta=2.5
mybins=np.zeros(nbbins+1)
binc=np.zeros(nbbins)
for b in range(nbbins+1):
    mybins[b]=b*delta

for b in range(nbbins):
    binc[b]=0.5*(mybins[b]+mybins[b+1])

nrow=len(df1)

# WATCHOUT ONLY ACUTE ANGLES!!!!!!!!!!!!!!!!!!!!!!!

all_angles=[]
conj_angles=[]
for index,row in df1.iterrows():
    if row.nb1 > nbmin and row.nb2 > nbmin:
        if row.clean_int:
            all_angles.append(row.int_angle)
            if row.conj_pair:
                conj_angles.append(row.int_angle)

nb_all=len(all_angles)
nb_conj=len(conj_angles)

print('total nb=', len(df1))
print('nb rejected=', len(df1)-nb_all)
print('nb all intersections=', nb_all)
print('nb conj pairs=', nb_conj)

mean_ang_int=np.mean(all_angles)
print('mean angle all intersection=', mean_ang_int)
mean_ang_conj=np.mean(conj_angles)
print('mean angle conjugate pairs=', mean_ang_conj)

plt.figure(1)
counts, bins, bars = plt.hist(all_angles, bins=mybins, density=True, color = "dodgerblue", ec="dodgerblue")
plt.xlabel('angle of intersection', fontsize=14)
plt.ylabel('PDF', fontsize=14)
plt.ylim(0,0.06)
plt.xlim(0,90)
fileout1='FIGS/PDF_int_angle_'+EXP+'_'+year+'_delta' + dlabel +'.png'
plt.savefig(fileout1)

plt.figure(2)
counts, bins, bars = plt.hist(conj_angles, bins=mybins, density=True, color = "dodgerblue", ec="dodgerblue")
plt.xlabel('angle of conjugate pair', fontsize=14)
plt.ylabel('PDF', fontsize=14)
plt.ylim(0,0.06)
plt.xlim(0,90)
fileout2='FIGS/PDF_conj_angle_'+EXP+'_'+year+'_delta' + dlabel +'.png'
plt.savefig(fileout2)

plt.show()
