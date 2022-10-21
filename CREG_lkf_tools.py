import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path
from lkf_tools.dataset import *

#----  CREG_lkf_detect --------------------------------------
#
# Prepares CREG netcdf outputs to be used by Nils' LKF
# detection algorithm. Also offers possibility to plot 
# detected LKFs on pan-Arctic grid (produce_plot=True).
#
#------------------------------------------------------------

def CREG_lkf_detect(date, creggrid, grid_path, data_path, store_path, fileout, kvalue, produce_plot):

    print(fileout)

    cregflag=1 # used to customize code for CREG applications.

#----- open netcdf file -----

    creg_nc = xr.open_dataset(data_path)

# WATCHOUT: in code U,V,A, shr,div and vor are collocated. The lat,lon at these 
#           points are ULAT and ULON...In our case all these variables are at the
#           T points...this is why TLAT,TLON are renamed ULAT, ULON. 

    creg_nc = creg_nc.rename({'ULON':'LONTP', 'ULAT':'LATTP'})
    creg_nc = creg_nc.rename({'divu':'div', 'shear':'shr', 'aice':'A', 
                              'uvel':'U', 'vvel':'V', 'TLON':'ULON', 'TLAT':'ULAT'})

#----- open grid coordinate file -----

    grid_nc = xr.open_dataset(grid_path)
    grid_nc = grid_nc.rename({'e1t':'DXU', 'e2t':'DYV'})

    creg_nc = xr.Dataset.merge(creg_nc, grid_nc)

#creg_nc = creg_nc.rename({'ni':'x', 'nj':'y','divu':'div', 'shear':'shr', 'aice':'A', 
#                          'uvel':'U', 'vvel':'V', 'TLON':'ULON', 'TLAT':'ULAT'})

#---- process data and detect LKFs ---

# il y presentement un bug dans les sorties des i,j. i est lkf[:,0] et j lkf[:,1].
# Voir courriel de Nils du 4 oct 2022. Pour corriger les i,j je dois faire:
# i = lkf[:,0] + lkf_data.index_x[0][0]
# j = lkf[:,1] + lkf_data.index_y[0][0]
#
# pour creg025:
# lkf_data.index_x[0][0]=93
# lkf_data.index_y[0][0]=329
# 
# pour creg12:
# lkf_data.index_x[0][0]=278
# lkf_data.index_y[0][0]=985

    print('call process_dataset')

    lkf_data = process_dataset(fileout,creg=cregflag,  output_path=store_path,
                               xarray=creg_nc, skeleton_kernel=kvalue, t_red=1)

    lkf_data.detect_lkfs(indexes=[0])

#---- plot LKFs ---------------------

    if (produce_plot):
        
        fig = plt.figure(figsize=[10, 5])

        ax = plt.subplot(1, 1, 1, projection=ccrs.Orthographic(0, 90))

        ax.coastlines(zorder=3)

        pcm = ax.pcolormesh(lkf_data.lon[max([0,lkf_data.index_y[0][0]-1]):lkf_data.index_y[0][-1]+2:lkf_data.red_fac,
                                         max([0,lkf_data.index_x[0][0]-1]):lkf_data.index_x[0][-1]+2:lkf_data.red_fac],
                            lkf_data.lat[max([0,lkf_data.index_y[0][0]-1]):lkf_data.index_y[0][-1]+2:lkf_data.red_fac,
                                         max([0,lkf_data.index_x[0][0]-1]):lkf_data.index_x[0][-1]+2:lkf_data.red_fac],
                            np.sum(lkf_data.eps_tot_list,axis=0),transform=ccrs.PlateCarree(),vmin=0,vmax=1e-1,cmap='Greys_r')

        it = lkf_data.indexes[-1]

        lkfs = np.load(lkf_data.lkfpath.joinpath('lkf_%s_%03i.npy' %(lkf_data.netcdf_file.split('/')[-1].split('.')[0],(it+1))),allow_pickle=True)

        i=0
        for ilkf in lkfs:
            if i == 60:
                print(ilkf.shape)

            i=i+1
            if np.min(ilkf[:,2])<-150 and np.max(ilkf[:,2]>150):
                ilkf[ilkf[:,2]<0,2]+=360
            ax.plot(ilkf[:,2],ilkf[:,3],transform=ccrs.PlateCarree())

        print('saving image')
        filefig='mapLKF_'+creggrid+'_'+date+.'png'
        plt.colorbar(pcm,label='total deformation')
        plt.savefig(filefig)

#----  CREG_lkf_analysis --------------------------------------
#
# Analyses detected LKFs in order to calculate half widths
# of LKFs. Detected LKF points have maximum values of eps_tot. 
# Given a certain direction of an LKF, the half widths are 
# then calculated in perpendicular directions referred to as 
# search directions. The search is done for a max number of 
# cells (dsearch) in a given direction. The half widths are 
# obtained when eps_tot < frac * (eps_tot_max). Half widths 
# are only calculated when the distance of a LKF point 
# (grid cell) from land is larger than mindist.
#
#------------------------------------------------------------

def CREG_lkf_analysis(date,creggrid,path_filedist,path_filein,path_fileout,data_path,dsearch,frac,mindist):
    
    print('working on date:')
    print(date)

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

#---- analyse  lkfs -----

    l=0
    for ilkf in lkfs:
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

#----  CREG_lkf_concatenate_width ---------------------------
#
# Concatenate half widths of a given LKF file in two vectors.
# 
#------------------------------------------------------------

def CREG_lkf_concatenate_width (date,path_filein, hwidth):
    
# ilkf[n,7] = hwidth1
# ilkf[n,8] = hwidth2

    print(date)
    tpvect=[]
    l=0
    lt=0
    lr=0
    lkfs = np.load(path_filein,allow_pickle=True)
    print(lkfs.shape)
    
    for ilkf in lkfs:
        nb=ilkf.shape[0] # nb of points in LKF i
        for n in range(nb):
            lt=lt+1
            if hwidth == 1:
                if ilkf[n,7] > 0:
                    tpvect.append(ilkf[n,7])
                    l=l+1
                else:
                    lr=lr+1
            elif hwidth == 2:
                if ilkf[n,8] > 0:
                    tpvect.append(ilkf[n,8])
                    l=l+1
                else:
                    lr=lr+1
                    
            else:
                print('wrong hwidth value') 
            
    print('final number of LKF points:')
    print(l)
    print('number of rejected LKF points:')
    print(lr)
    
    return tpvect
