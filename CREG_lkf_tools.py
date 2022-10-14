import xarray as xr
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path
from lkf_tools.dataset import *

def CREG_lkf_detect(date, grid_path, data_path, store_path, fileout, kvalue, produce_plot):

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
        plt.colorbar(pcm,label='total deformation')
        plt.savefig('testing12.png')
