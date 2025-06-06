#load_ext autoreload
#autoreload 2

import numpy as np
import xarray as xr
import os
from pathlib import Path
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

#----- input -----

#ni = 528 ; creg025
#nj = 735 ;

#ni = 1580 ; creg12
#nj = 2198 ;

creggrid='creg025' # creg025 or creg12
EXP='run6f'
ddate='2005032900_000'
lkfchoice=62
main_dir='/home/jfl001/data/Lemieux2022/LKF_diag'
main_dirnc='/home/jfl001/data/runsLemieux_et_al_2022/'
dir_util='/home/jfl001/Lemieux2022/UTIL'
dsearch=5 # +- dsearch cells around one LKF cell (dist is capped if searching too far!!!)
frac=0.5 # half width is defined as eps_tot < frac*LKFepsmax 
mindist=150.0 # LKF point is analysed if dist from land > mindist (km)

#---- exemples interessants 2005032903_000 run6f kvalue=7----
# 46: ligne horizontale
# 50: longue ligne diagonale
# 62: longue ligne diagonale de faibles def
# 63: courte ligne diagonale parfaite de fortes def
# 69: ligne verticale

#---- begin -------------------------------

print('working on date:')
print(ddate)

#----- define paths and file names --------

filein='lkf_' + ddate + '_' + EXP + '_001.npy'
fileout='lkf_' + ddate + '_' + EXP + '_width.npy'
tpdir=ddate + '_' + EXP
path_filein=os.path.join(main_dir+'/'+creggrid+'/'+EXP+'/'+tpdir+'/'+filein)
path_fileout=os.path.join(main_dir+'/'+creggrid+'/'+EXP+'/'+tpdir+'/'+fileout)
data_path=os.path.join(main_dirnc+creggrid+'/'+EXP+'/netcdf/'+ddate+'.nc')
path_filedist=os.path.join(dir_util +'/dist_'+creggrid+'.pkl')

#----- open npy file -----

lkfs = np.load(path_filein,allow_pickle=True)
print(lkfs.shape)

#----- open nc file and calc tot def -----
# array is opened as j,i in python

creg_nc = xr.open_dataset(data_path)
div = creg_nc.divu[0,:,:]/100.0
shr = creg_nc.shear[0,:,:]/100.0
eps_tot = np.sqrt(div**2+shr**2)

#----- open distance to land file --------

dist = np.load(path_filedist,allow_pickle=True)

#----- shift indices ---------------------

# arrays (i,j) are read in python as (j,i)

# il y presentement un bug dans les sorties des j,i. 
# les indices ne correspondent pas aux indices de la grille native.
# dans ce code: jl,il indices des LKFs (avec bug) et j,i indices grille native
# jl=lkf[:,0] et il=lkf[:,1].
# Voir courriel de Nils du 4 oct 2022. Pour corriger les i,j je dois faire:
# i = lkf[:,0] + lkf_data.index_x[0][0]
# j = lkf[:,1] + lkf_data.index_y[0][0]

# pour creg025:
# lkf_data.index_x[0][0]=93
# lkf_data.index_y[0][0]=329
# 
# pour creg12:
# lkf_data.index_x[0][0]=278
# lkf_data.index_y[0][0]=985

# je pense que ce que Nils a écrit n'est pas ok. Ça devrait être:

# j = lkf[:,0] + lkf_data.index_y[0][0] - 1
# i = lkf[:,1] + lkf_data.index_x[0][0] - 1

if (creggrid == 'creg025'):
    jshift=329
    ishift=93
elif (creggrid == 'creg12'):
    jshift=985
    ishift=278
else:
    print ("Wrong choice of grid")

#---- analyse chosen lkf -----

l=0
for ilkf in lkfs:
#    print(l)
#    print(ilkf.shape)
#    if l == lkfchoice:
    nb=ilkf.shape[0] # nb of points in LKF i
    maxjl=np.max(ilkf[:,0])
    minjl=np.min(ilkf[:,0])
    maxil=np.max(ilkf[:,1])
    minil=np.min(ilkf[:,1])

    zlkf=ilkf

#---- find search direction and width (MV to function...)--------------------
#  
# The search direction is found by defining 3 vectors (v,v1 and v2). 
# v is in the diection of the LKF from n=0...n=nb. v1 and v2 indicate the 
# search direction (there is 180 deg between v1 and v2). v1 is obtained by a 
# 'regle de la main droite'. The index (nail up) points in the direction of v 
# and the thumb indicates the direction of v1. 
# jl = zlkf[:,0], il=zlkf[:,1]
#----------------------------------------------------------------------------

    halfw1=np.zeros(nb)
    halfw2=np.zeros(nb)
 
    for n in range(nb):

        jl = int(zlkf[n,0])
        il = int(zlkf[n,1])
        j=jl+jshift-1
        i=il+ishift-1 

        if (n == 0):
            deltj=zlkf[1,0]-zlkf[0,0]
            delti=zlkf[1,1]-zlkf[0,1]
        elif (n == nb-1):
            deltj=zlkf[n,0]-zlkf[n-1,0]
            delti=zlkf[n,1]-zlkf[n-1,1]
        else:
            deltj=zlkf[n+1,0]-zlkf[n-1,0]
            delti=zlkf[n+1,1]-zlkf[n-1,1]
        
# v vector is written as av,bv = av x^ + bv y^ (x^,y^=unit vectors along i and j)
# v, v1 and v2 are either (1,0),(0,1),(-1,0),(0,-1),(1,1),(-1,1),(-1,-1),(1,-1)
# be careful vector as x^ component first!!!
    
        sdelt=abs(delti)+abs(deltj) # sum of |delti|+|deltj|
        sdelt=int(sdelt)
        deno=max(0.1,abs(delti))
        av=delti/deno
        deno=max(0.1,abs(deltj))
        bv=deltj/deno

# width is analysed if sdelt=1, 2 or 4. These corresponds to vectors v along 
# x^, y^ or at 45 deg.

        if (sdelt == 1 or sdelt == 2 or sdelt == 4 and dist[j,i] > mindist):
            av1=int(-bv) # v1 and v2 are found by rotating v by +-90 deg
            bv1=int(av)
            av2=int(bv)
            bv2=int(-av)
            LKFepsmax=eps_tot[j,i]
            target=frac*LKFepsmax

# set initial values (overwritten if found while s < dsearch)

            idelta=int(dsearch*av1)
            jdelta=int(dsearch*bv1)
            halfw1[n]=np.sqrt(idelta**2 + jdelta**2)
            idelta=int(dsearch*av2)
            jdelta=int(dsearch*bv2)
            halfw2[n]=np.sqrt(idelta**2 + jdelta**2)

# along v1
            s=0
            while s < dsearch:
                idelta=int((s+1)*av1)
                jdelta=int((s+1)*bv1)
                jj=j+jdelta
                ii=i+idelta
                if eps_tot[jj,ii] < target:
                    halfw1[n]=np.sqrt(idelta**2 + jdelta**2)
                    break
                s=s+1

# along v2
            s=0
            while s < dsearch: 
                idelta=int((s+1)*av2)
                jdelta=int((s+1)*bv2)
                jj=j+jdelta
                ii=i+idelta
                if eps_tot[jj,ii] < target:
                    halfw2[n]=np.sqrt(idelta**2 + jdelta**2)
                    break
                s=s+1        

#        tpcalc= np.sqrt(zlkf[n,4]**2 + zlkf[n,5]**2)
#        print(tpcalc) # tpcalc should be equal tp LKFepsmax...it is        
        
        else:
            halfw1[n]=np.nan
            halfw2[n]=np.nan

    zlkf=np.c_[zlkf, halfw1] # add halfw vectors to zlkf
    zlkf=np.c_[zlkf, halfw2]
    lkfs[l,]=zlkf        

    l=l+1

#----- save output file -----

np.save(path_fileout,lkfs,allow_pickle=True)
