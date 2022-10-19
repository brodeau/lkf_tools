#load_ext autoreload
#autoreload 2

import numpy as np
#import xarray as xr
import os
from pathlib import Path
import matplotlib.pyplot as plt
#import cartopy.crs as ccrs

creggrid='creg025' # creg025 or creg12
EXP='run6f'
SDATE='20050329'
EDATE='20050329'
main_dir='/home/jfl001/data/Lemieux2022/LKF_diag'
mybins=np.linspace(0, 8, 100)

#----- define paths and file names --------

filein1='hwidth1_lkf_'+SDATE+'_'+EDATE+'.npy'
path_filein1=os.path.join(main_dir+'/'+creggrid+'/'+EXP+'/WIDTH/'+filein1)
filein2='hwidth2_lkf_'+SDATE+'_'+EDATE+'.npy'
path_filein2=os.path.join(main_dir+'/'+creggrid+'/'+EXP+'/WIDTH/'+filein2)

#----- open npy files -----

hwidth1 = np.load(path_filein1,allow_pickle=True)
hwidth2 = np.load(path_filein2,allow_pickle=True)

#n=0
#hconca1=np.zeros(1060)
#for hw in hwidth1:
#    for hwi in  hw:
#        hconca1[n]=hwi
#        n=n+1
        

#print(n)
#n=0
#hconca2=np.zeros(1060)
#for hw in hwidth2:
#    for hwi in  hw:
#        hconca2[n]=hwi
#        n=n+1


plt.hist(hwidth1, bins=mybins, color = "dodgerblue", ec="dodgerblue")
plt.hist(hwidth2, bins=mybins, alpha=0.3, color = "magenta", ec="magenta")
#plt.hist(Ddef2, bins=mybins, alpha=0.3, color = "orange", ec="orange", )
#plt.yscale('log')
#plt.xlabel(r'$\Delta \dot{\epsilon_{II}}$', fontsize=14)
#plt.ylabel('Nb of counts', fontsize=12)
#plt.legend([label1, label2], loc ="upper right", markerscale=1)
plt.savefig('hwidth2.png')



