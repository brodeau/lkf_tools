import os,sys
import numpy as np
import pandas as pd
#import matplotlib as plt
import matplotlib.pyplot as plt
#from datetime import timedelta
#import pickle
#import calendar

EXP='eg1p5_ef1p5'
main_dir='/home/jfl001/data/Lemieux_et_al_plast_pot/LKF_diag'
zdir='Angle_grid'
SDATE='2006030100'
#EDATE='2005022800'
addlabel='anggrid'

#-----------------------------------------

filename=SDATE+'_'+addlabel+'_'+EXP+'.py'
filein=os.path.join(main_dir+'/'+EXP+'/'+zdir+'/'+filename)
print(filein)

#path_filein=os.path.join(densitydir+filein)

#density = np.load(path_filein,allow_pickle=True)

flights = pd.read_csv(filein)
print(flights)
        
nbbins=90
delta=1.0
mybins=np.zeros(nbbins+1)
binc=np.zeros(nbbins)
for b in range(nbbins+1):
    mybins[b]=b*delta

for b in range(nbbins):
    binc[b]=0.5*(mybins[b]+mybins[b+1])

print(mybins)
print(binc)
