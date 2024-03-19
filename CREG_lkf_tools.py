import numpy as np
import xarray as xr
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path
#from scipy.optimize import fsolve
from shapely.geometry import LineString
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

    creg_nc = creg_nc.rename({'divu':'div', 'shear':'shr', 'vort':'vor', 'aice':'A', 
                              'uvel':'U', 'vvel':'V'})

#----- open grid coordinate file -----

# WATCHOUT: in the code here U,V,A, shr,div and vor are collocated. 
#           All are ok for my creg12 runs except uvel and vvel (U point) but they
#           are not needed for detect. I would have to check this if I use lkf_tracking
#           ULAT and ULON here are in fact at the T point. This is why I use nav_lat and 
#           nav_lon to define these. 

    grid_nc = xr.open_dataset(grid_path)
    grid_nc = grid_nc.rename({'e1t':'DXU', 'e2t':'DYV', 'nav_lon':'ULON', 'nav_lat':'ULAT'})

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
        filefig='mapLKF_'+creggrid+'_'+date+'.png'
        plt.colorbar(pcm,label='total deformation')
        plt.savefig(filefig)

#----  CREG_lkf_calc_width ----------------------------------
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

def CREG_lkf_calc_width(date,creggrid,path_filedist,path_filein,path_fileout,data_path,dsearch,frac,mindist):
    
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

#----  CREG_lkf_density -------------------------------------
#
# Called by CREG_driver_LKF_density to calculate contribution 
# to density from a single LKF file (containing many detected
# LKFs). There is no criterion applied here for the distance 
# to land. This could be done later when plotting the density. 
#
#------------------------------------------------------------

def CREG_lkf_density(date,creggrid,path_filein):
    
    print('working on date:')
    print(date)

#----- open npy file -----

    lkfs = np.load(path_filein,allow_pickle=True)
    print(lkfs.shape)

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
        nx=528
        ny=735
        jshift=329
        ishift=93
    elif (creggrid == 'creg12'):
        nx=1580
        ny=2198
        jshift=985
        ishift=278
    else:
        print ("Wrong choice of grid")

#---- calc density contribution from group of LKFs -----

    TPdens= np.zeros((ny,nx))

    l=0
    for ilkf in lkfs:
        nb=ilkf.shape[0] # nb of points in LKF i
        zlkf=ilkf

        for n in range(nb):

            jl = int(zlkf[n,0])
            il = int(zlkf[n,1])
            j=jl+jshift-1
            i=il+ishift-1 
            TPdens[j,i]=TPdens[j,i]+1.0
        l=l+1

    return TPdens


#------------------------------------------------------------
#  Functions for CREG_lkf_pairs_and_angles
#------------------------------------------------------------

#---- identify type of intersection -------------------------

def identify_int(index1,nb1,index2,nb2):
# index1 and index 2 identify element closest to intersection pt
    if index1<=1 or index2<=1 or index1>=nb1-1 or index2>=nb2-1:
        int_type=2 # T or Y intersecting LKFs
    else:
        int_type=1 # X intersecting LKFs...possible conjugate faults    

    return int_type

#---- calculate intersection angle --------------------------

def calc_int_angle(ptype1,coeff1,ptype2,coeff2):
    m1=coeff1[0]
    m2=coeff2[0]
    minval=1e-12    

    if ptype1==1 and ptype2==2:
        denom=max(minval,abs(m2))
        m2=np.sign(m2)*1.0/denom # convert dx/dy to -dy/dx 
    elif ptype1==2 and ptype2==1:
        denom=max(minval,abs(m1))
        m1=np.sign(m1)*1.0/denom # convert dx/dy to -dy/dx 
 
    deno=max(minval,abs(1+m1*m2)) # avoid div by zero
    int_angle=np.arctan(abs((m1-m2)/deno))
    int_angle = int_angle*180/np.pi # convert to deg

    return int_angle

#---- check if two rectangles overlap -----------------------

def overlap(ibl1, jbl1, itr1, jtr1, ibl2, jbl2, itr2, jtr2):
    if itr1 < ibl2 or ibl1 > itr2 or jtr1 < jbl2 or jbl1 > jtr2:
        ovlflag=False
    else:
        ovlflag=True

    return ovlflag

#---- get i,j at intersection point -----------------------

def get_ij_intersection(intersec):
    if intersec.geom_type == 'Point':
        iint=intersec.x
        jint=intersec.y
    elif intersec.geom_type == 'LineString': #take 1st pt of string
        iint=intersec.xy[0][0]
        jint=intersec.xy[1][0]
    elif intersec.geom_type == 'MultiPoint': #take 1st pt of multipt
        iint=intersec.geoms[0].x
        jint=intersec.geoms[0].y
    else:
        iint=intersec.geoms[0].xy[0][0]
        jint=intersec.geoms[0].xy[1][0]
        #print('wowowo',intersec.geoms[0].xy[0][0],intersec.geoms[0].xy[1][0])

    return iint,jint

#---- polyfit over intersection zone ------------------------

def get_polyfit(vari, varj, xf, yf, pdeg, nbsub):
    if vari >= varj: # y=p(x)
        coeff = np.polyfit(xf,yf,pdeg)
        yfunc = np.poly1d(coeff)
        xpf = np.linspace(xf[0],xf[nbsub],50)
        ypf=yfunc(xpf)
        ptype=1
    else:
        coeff = np.polyfit(yf,xf,pdeg) # x=p(y)
        xfunc = np.poly1d(coeff)
        ypf = np.linspace(yf[0],yf[nbsub],50)
        xpf=xfunc(ypf)
        ptype=2
        
    return xpf,ypf,ptype,coeff

#---- add extra pt at start --------------------------------

def extra_pt_start(i0, i1, j0, j1):
    delti=i1-i0
    deltj=j1-j0
    ist=i0-delti
    jst=j0-deltj
        
    return ist,jst

def extra_pt_end(im1, im2, jm1, jm2):
    delti=im1-im2
    deltj=jm1-jm2
    iend=im1+delti
    jend=jm1+deltj
        
    return iend,jend

#----  CREG_lkf_pairs_and_angles ----------------------------
#
# Identifies pairs of LKFs that intersect and calc the 
# intersection angle. 
#
#------------------------------------------------------------

def CREG_lkf_pairs_and_angles(date,creggrid,path_filein,data_pathnc):
    
    print('working on date:')
    print(date)

    pdeg=1 # degree of polynomial for fit
    vort_val=0.0 # value used to extrapolate vort at both ends of LKF

#----- open npy file -----

    lkfs = np.load(path_filein,allow_pickle=True)
    print(lkfs.shape)

#----- open netcdf file --

    creg_nc = xr.open_dataset(data_pathnc)
    vort = creg_nc.vort[0,:,:]/100.0

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
        nx=528
        ny=735
        jshift=329
        ishift=93
    elif (creggrid == 'creg12'):
        nx=1580
        ny=2198
        jshift=985
        ishift=278
    else:
        print ("Wrong choice of grid")

#---- identify pairs of intersecting LKFs -----

    npot=0
    nt=0

#        j=jl+jshift-1
#        i=il+ishift-1

    for ind1, lkf1 in enumerate(lkfs):
        ntp1=lkf1.shape[0] # nb of points in LKF
        maxj1=np.max(lkf1[:,0])
        minj1=np.min(lkf1[:,0])
        maxi1=np.max(lkf1[:,1])
        mini1=np.min(lkf1[:,1])
        j1=lkf1[:,0]
        i1=lkf1[:,1]

        #--- vorticity along LKF1 ---
        vort1=np.zeros(ntp1)
        for n1 in range(ntp1):
            jj=int(j1[n1])+jshift-1
            ii=int(i1[n1])+ishift-1
            vort1[n1]=vort[jj,ii]

        #--- extrapolate vort1 by one pt on both sides to match shape of j1ext,i1ext ---
        vort1=np.append(vort_val,vort1) # start
        vort1=np.append(vort1,vort_val) # end

        #line1=LineString(np.column_stack((i1,j1)))
        
        #--- extrapolate LKF1 by one pt on both sides ---
        ist,jst=extra_pt_start(i1[0], i1[1], j1[0], j1[1])
        j1ext=np.append(jst,j1)
        i1ext=np.append(ist,i1)

        iend,jend=extra_pt_end(i1[-1], i1[-2], j1[-1], j1[-2])
        j1ext=np.append(j1ext,jend)
        i1ext=np.append(i1ext,iend)
        nb1=i1ext.shape[0] # nb of points in LKF
        line1ext=LineString(np.column_stack((i1ext,j1ext)))

        for ind2, lkf2 in enumerate(lkfs):
            if ind2 > ind1:
                ntp2=lkf2.shape[0] # nb of points in LKF
                maxj2=np.max(lkf2[:,0])
                minj2=np.min(lkf2[:,0])
                maxi2=np.max(lkf2[:,1])
                mini2=np.min(lkf2[:,1])
#                ovlflag=overlap(mini1, minj1, maxi1, maxj1, mini2, minj2, maxi2, maxj2)
                ovlflag=overlap(mini1-1, minj1-1, maxi1+1, maxj1+1, mini2-1, minj2-1, maxi2+1, maxj2+1)
                if ovlflag:
                    npot=npot+1
                    j2=lkf2[:,0]
                    i2=lkf2[:,1]
                    #line2=LineString(np.column_stack((i2,j2)))
                    
                    #--- extrapolate LKF1 by one pt on both sides ---
                    ist,jst=extra_pt_start(i2[0], i2[1], j2[0], j2[1])
                    j2ext=np.append(jst,j2)
                    i2ext=np.append(ist,i2)

                    iend,jend=extra_pt_end(i2[-1], i2[-2], j2[-1], j2[-2])
                    j2ext=np.append(j2ext,jend)
                    i2ext=np.append(i2ext,iend)
                    nb2=i2ext.shape[0]
                    line2ext=LineString(np.column_stack((i2ext,j2ext)))

                    #intersec=line1.intersection(line2) # inters. pt between line1,2
                    intersec=line1ext.intersection(line2ext) # inters. pt between line1,2
                    dlt=5
                    if not intersec.is_empty:

                        iint,jint=get_ij_intersection(intersec) # get i,j at intersection point
                        nt=nt+1
                        deltai=abs(i1ext-iint) # delta i between i1 and intersec i
                        deltaj=abs(j1ext-jint) # delta i between j1 and intersec j
                        index1=np.argmin(deltai+deltaj)
                        deltai=abs(i2ext-iint) # delta i between i2 and intersec i
                        deltaj=abs(j2ext-jint) # delta i between j2 and intersec j
                        index2=np.argmin(deltai+deltaj)
                        
                        #--- define part of array close to intersec for polyfit
                        min_ind1=max(0, index1-dlt)
                        max_ind1=min(index1+dlt,nb1)
                        xf1=i1ext[min_ind1:max_ind1+1]
                        nbsub1=xf1.shape[0]-1
                        yf1=j1ext[min_ind1:max_ind1+1]
                        vari1=max(xf1)-min(xf1) # variation of i1 in pts used for polyfit                    
                        varj1=max(yf1)-min(yf1)

                        xpf1,ypf1,ptype1,coeff1=get_polyfit(vari1,varj1,xf1,yf1,pdeg,nbsub1) # polyfit LKF1
      
                        min_ind2=max(0, index2-dlt)
                        max_ind2=min(index2+dlt,nb2)
                        xf2=i2ext[min_ind2:max_ind2+1]
                        nbsub2=xf2.shape[0]-1
                        yf2=j2ext[min_ind2:max_ind2+1]
                        vari2=max(xf2)-min(xf2) # variation of i2 in pts used for polyfit
                        varj2=max(yf2)-min(yf2)

                        xpf2,ypf2,ptype2,coeff2=get_polyfit(vari2,varj2,xf2,yf2,pdeg,nbsub2) # polyfit LKF2
                        
                        #--- identify if intersection is X, T or Y
                        int_type=identify_int(index1,nb1,index2,nb2)

                        #--- calc intersection angle (returns acute angle)
                        int_angle=calc_int_angle(ptype1,coeff1,ptype2,coeff2)
                        #### NEED TO CHECK IF THEY INTERSECT or should I???###
                        # ELIMINATE MULTIPT???
                        
                        if int_type==1: # possible conjugate fault lines
                            #--- vorticity along LKF2 ---
                            vort2=np.zeros(ntp2)
                            for n2 in range(ntp2):
                                jj=int(j2[n2])+jshift-1
                                ii=int(i2[n2])+ishift-1
                                vort2[n2]=vort[jj,ii]

                                #--- extrapolate vort1 by one pt on both sides to match shape of j1ext,i1ext ---
                                vort2=np.append(vort_val,vort2) # start
                                vort2=np.append(vort2,vort_val) # end

                        print(ind1,ind2)
                        cc=0
                        figg=1
                        if ind1 == 140 and ind2 == 144:
                            print('ptype',ptype1,ptype2,int_type)
                            print(index1,nb1,index2,nb2)
                            print(coeff1[0], coeff2[0])
                            print('angle',int_angle)
                            print('ind1',min_ind1,index1,max_ind1,nb1)
                            print('ind2',min_ind2,index2,max_ind2,nb2)
                            if figg==1:
                                plt.plot(i1ext,j1ext,'.m')
                                plt.plot(line1ext.xy[0],line1ext.xy[1], '-m')
                                plt.plot( xpf1,ypf1,'g')
                                plt.plot(i2ext,j2ext,'*b')
                                plt.plot(line2ext.xy[0],line2ext.xy[1], '-b')
                                plt.plot( xpf2,ypf2,'orange')
                                plt.plot( intersec.x,intersec.y, 'sr')
                                if cc==1:
                                    plt.xlim(iint-10, iint+10)
                                    plt.ylim(jint-10, jint+10)
                                plt.show()
                            elif figg==2:
                                maxjtp=int(max(maxj1,maxj2))
                                minjtp=int(min(minj1,minj2))
                                maxitp=int(max(maxi1,maxi2))
                                minitp=int(min(mini1,mini2))
                                # shift indices from local grid to creg grid
                                maxj=maxjtp+jshift-1
                                minj=minjtp+jshift-1
                                maxi=maxitp+ishift-1
                                mini=minitp+ishift-1
                                print(minj,maxj,mini,maxi)
                                vortpc=vort[minj:maxj,mini:maxi]                                
                                plt.pcolor(vortpc)
                                plt.colorbar()
                                #plt.savefig('testing12.png')
                                plt.show()
