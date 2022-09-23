#load_ext autoreload
#autoreload 2

import xarray as xr
import os
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

lkf_data = process_dataset(fileout,creg=cregflag,
                           output_path=store_path, xarray=creg_nc, t_red=1)

lkf_data.detect_lkfs(indexes=[0])
