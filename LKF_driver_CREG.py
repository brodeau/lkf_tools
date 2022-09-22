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

#-----------------

data_path=os.path.join(main_dir+creggrid+'/'+EXP+'/netcdf/'+ddate+'.nc')
#print(data_path)

#----- open netcdf file -----

creg_nc = xr.open_dataset(data_path)

# WATCHOUT: in code U,V,A, shr,div and vor are collocated. The lat,lon at these 
#           points are ULAT and ULON...In our case all these variables are at the
#           T points...this is why TLAT,TLON are renamed ULAT, ULON. 

creg_nc = creg_nc.rename({'ULON':'LONTP', 'ULAT':'LATTP'})
creg_nc = creg_nc.rename({'divu':'div', 'shear':'shr', 'aice':'A', 
                          'uvel':'U', 'vvel':'V', 'TLON':'ULON', 'TLAT':'ULAT'})

#creg_nc = creg_nc.rename({'ni':'x', 'nj':'y','divu':'div', 'shear':'shr', 'aice':'A', 
#                          'uvel':'U', 'vvel':'V', 'TLON':'ULON', 'TLAT':'ULAT'})

#creg_nc = creg_nc.assign(A=np.isfinite(creg_nc.div).astype('float'))
#creg_nc.A[0,:,:].plot()

lkf_data = process_dataset(data_path,creg=cregflag,
                           output_path='/home/jfl001/data/LKF_diag', xarray=creg_nc)

lkf_data.detect_lkfs(indexes=[0])
