#load_ext autoreload
#autoreload 2

import xarray as xr
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path
from lkf_tools.dataset import *

#----- input -----

cregflag=1 # used to customize code for CREG applications.
creggrid='creg025' # creg025 or creg12
EXP='run6f'
ddate='2005032903_000'
main_dir='/home/jfl001/data/runsLemieux_et_al_2022/'
main_dir_grid='/home/socn000/data/eccc-ppp5/env_rhel-8-icelake-64/datafiles/constants/oce/repository/master/CONCEPTS/'
store_main_dir='/home/jfl001/data/Lemieux2022/LKF_diag'
kvalue=7
produce_plot=True

#----- define paths and file name --------

store_path=os.path.join(store_main_dir+'/'+creggrid+'/'+EXP+'/')
fileout=ddate + '_' + EXP

data_path=os.path.join(main_dir+creggrid+'/'+EXP+'/netcdf/'+ddate+'.nc')
if (creggrid == 'creg025'):
    grid_path=os.path.join(main_dir_grid+'/creg025pe/grid/coordinates_CREG025_LIM.nc')
elif (creggrid == 'creg12'):
    grid_path=os.path.join(main_dir_grid+'/creg012pe/grid/coordinates_CREG12_ext.nc')
else:
    print ("Wrong choice of grid")

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
#            print(ilkf)
        i=i+1
        if np.min(ilkf[:,2])<-150 and np.max(ilkf[:,2]>150):
            ilkf[ilkf[:,2]<0,2]+=360
        ax.plot(ilkf[:,2],ilkf[:,3],transform=ccrs.PlateCarree())

    plt.colorbar(pcm,label='total deformation')
    plt.savefig('testing12.png')

