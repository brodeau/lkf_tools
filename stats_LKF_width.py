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
mybins=np.linspace(0, 15, 100)

#----- define paths and file names --------

filein1='hwidth1_lkf_'+SDATE+'_'+EDATE+'.npy'
path_filein1=os.path.join(main_dir+'/'+creggrid+'/'+EXP+'/WIDTH/'+filein1)
filein2='hwidth2_lkf_'+SDATE+'_'+EDATE+'.npy'
path_filein2=os.path.join(main_dir+'/'+creggrid+'/'+EXP+'/WIDTH/'+filein2)

#----- open npy files -----

hwidth1 = np.load(path_filein1,allow_pickle=True)
hwidth2 = np.load(path_filein2,allow_pickle=True)

print(hwidth1.shape[0])
print(hwidth2.shape[0])
npoints=hwidth1.shape[0]
width=np.zeros(npoints)
width[:]=hwidth1[:]+hwidth2[:]
print('mean values for hw1, hw2 and width are:')
print(np.mean(hwidth1))
print(np.mean(hwidth2))
print(np.mean(width))

plt.hist(width, bins=mybins, color = "dodgerblue", ec="dodgerblue")
#plt.hist(hwidth2, bins=mybins, alpha=0.3, color = "magenta", ec="magenta")
#plt.hist(Ddef2, bins=mybins, alpha=0.3, color = "orange", ec="orange", )
#plt.yscale('log')
#plt.xlabel(r'$\Delta \dot{\epsilon_{II}}$', fontsize=14)
#plt.ylabel('Nb of counts', fontsize=12)
#plt.legend([label1, label2], loc ="upper right", markerscale=1)
plt.savefig('width.png')



