import matplotlib.pyplot as plt
import numpy as np
import os
import sys

import fitsio
from astropy.table import join,Table
import healpy as hp

from LSS.tabulated_cosmo import TabulatedDESI
cosmo = TabulatedDESI()
dis_dc = cosmo.comoving_radial_distance

outdir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/plots/'
qt = 'COMP_TILE'

nside = 256
nest = True
zcol = 'Z_not4clus'

def mknz(zin,wl,fcr,bs=0.01,zmin=0.01,zmax=1.6):
    #cd = distance(om,1-om)
    ranf = fitsio.read_header(fcr,ext=1) #should have originally had 2500/deg2 density, so can convert to area
    area = ranf['NAXIS2']/2500.
    print('area is '+str(area))
    
    
    nbin = int((zmax-zmin)/bs)
    zhist = np.histogram(zin,bins=nbin,range=(zmin,zmax),weights=wl)
    zs = []
    nz = []
    for i in range(0,nbin):
        zl = zhist[1][i]
        zh = zhist[1][i+1]
        zm = (zh+zl)/2.
        zs.append(zm)
        voli = area/(360.*360./np.pi)*4.*np.pi/3.*(dis_dc(zh)**3.-dis_dc(zl)**3.)
        nbarz =  zhist[0][i]/voli
        nz.append(nbarz)
    return zs,nz


tps = ['QSO','LRG','BGS_ANY','BGS_BRIGHT','ELG','ELG_LOP','ELG_LOPnotqso']
for tp in tps:
    rf = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+tp+'_0_full.ran.fits'
    rt = fitsio.read_header(rf,ext=1)
    area = rt['NAXIS2']/2500
    dt = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+tp+'_full.dat.fits')

    wz = dt['ZWARN']*0 == 0
    wz &= dt['ZWARN'] != 1.e20
    wz &= dt['ZWARN'] != 999999
    wz &= dt['LOCATION_ASSIGNED'] == 1

    if tp == 'QSO':
        #good redshifts are currently just the ones that should have been defined in the QSO file when merged in full
        wg = dt[zcol]*0 == 0
        wg &= dt[zcol] != 999999
        wg &= dt[zcol] != 1.e20
    
    if tp[:3] == 'ELG':
        wg = dt['o2c'] > 0.9

    if tp == 'LRG':
        # Custom DELTACHI2 vs z cut from Rongpu
        #wg = dt['ZWARN'] == 0
        #drz = (10**(3 - 3.5*dt[zcol]))
        #mask_bad = (drz>30) & (dt['DELTACHI2']<30)
        #mask_bad |= (drz<30) & (dt['DELTACHI2']<drz)
        #mask_bad |= (dt['DELTACHI2']<10)
        #wg &= dt[zcol]<1.4
        #wg &= (~mask_bad)
        wg = dt['DELTACHI2'] > 15
        wg &= dt['ZWARN'] == 0
        wg &= dt[zcol]<1.5


    if tp[:3] == 'BGS':
        wg = dt['DELTACHI2'] > 40
    print(tp+':')
    print('area: '+str(area))
    print('# of good obs: '+str(len(dt[wz])))
    print('# of good z: '+str(len(dt[wz&wg])))
    print('completeness: '+str(round(len(dt[wz])/len(dt),3)))

    zl = dt[wg&wz][zcol]
    wl = 1./dt[wg&wz]['FRACZ_TILELOCID']
    bs=0.02
    zmin=0.01
    zmax=1.61
    if tp == 'QSO':
        zmin = 0.6
        zmax = 4.5
        dz = 0.05
    if tp[:3] == 'BGS':
        zmin = 0.01
        zmax = 1
        dz = 0.02
        plt.ylim(0,.05)
    
    svdir = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/fuji/LSScats/3.1/'
    svfs = svdir + tp+'_S_nz.txt'
    svfn = svdir + tp+'_N_nz.txt'
    if tp == 'LRG':
        svfs = svdir + 'LRG_main_S_nz.txt'
        svfn = svdir + 'LRG_main_N_nz.txt'
    if tp == 'ELG_LOP':
        svfs = svdir + 'ELG_HIP_S_nz.txt'
        svfn = svdir + 'ELG_HIP_N_nz.txt'
    if tp ==  'ELG_LOPnotqso': 
        svfs = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/LSScats/test/ELG_HIPnotqso_S_nz.dat'
        svfn = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/LSScats/test/ELG_HIPnotqso_N_nz.dat'
    svzs = np.loadtxt(svfs).transpose()  
    svzn = np.loadtxt(svfn).transpose()   
    plt.plot(svzs[0],(svzs[3]+svzn[3])/2.,label='SV3 fuji, v3.1')

    zm,nz = mknz(zl,wl,rf,bs=bs,zmin=zmin,zmax=zmax)
    if tp[:3] == 'ELG':
        plt.ylim(0,0.0015)
    plt.plot(zm,nz,label='main survey daily')
    plt.xlabel('Z')
    plt.ylabel(r'$n(z)~((h/$Mpc)$^3)$',labelpad=0)
    plt.legend()
    plt.title(tp)
    plt.grid(True)
    
    del dt

    plt.savefig(outdir+tp+'_nz.png')
    plt.clf()

    print(tp+' done')