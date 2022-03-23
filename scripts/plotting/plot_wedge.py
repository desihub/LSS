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

nside = 256
zcol = 'Z_not4clus'
ram = 130
rax = 150
ra0 = (ram+rax)/2.
decm = -1
decx = 1
zmin = 0
zmax = 3.5


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


tps = ['QSO','LRG','BGS_ANY','ELG']
cl = ['y','r','g','b']
for tp,c in zip(tps,cl):
    cols = ['RA','DEC',zcol,'ZWARN','DELTACHI2','LOCATION_ASSIGNED']
    if tp == 'ELG':
        cols.append('o2c')
    dt = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+tp+'zdone_full.dat.fits',columns=cols)
    sel = dt['RA'] > ram
    sel &= dt['RA'] < rax
    sel &= dt['DEC'] > decm
    sel &= dt['DEC'] < decx
    #sel &= dt[zcol] < zmax
    #sel &= dt[zcol] > zmin

    dt = dt[sel]
    
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
        wg = dt['ZWARN'] == 0
        drz = (10**(3 - 3.5*dt[zcol]))
        mask_bad = (drz>30) & (dt['DELTACHI2']<30)
        mask_bad |= (drz<30) & (dt['DELTACHI2']<drz)
        mask_bad |= (dt['DELTACHI2']<10)
        wg &= dt[zcol]<1.4
        wg &= (~mask_bad)


    if tp[:3] == 'BGS':
        wg = dt['DELTACHI2'] > 40
    print(tp+':')
    
    print('# of good obs: '+str(len(dt[wz])))
    print('# of good z: '+str(len(dt[wz&wg])))
    print('completeness: '+str(round(len(dt[wz])/len(dt),3)))

    dt = dt[wg&wz]
    sel = dt[zcol] < zmax
    sel &= dt[zcol] > zmin
    dt = dt[sel]

    
    r = dis_dc(dt[zcol])
    th = (90-dt['DEC'])*np.pi/180.
    phi = (dt['RA']-ra0)*np.pi/180
    x = r*np.cos(phi)*np.sin(th)
    y = r*np.sin(phi)*np.sin(th)
    z = r*np.cos(th)
    plt.plot(x,y,',',color=c)
    #plt.show()
    
    
    del dt


    print(tp+' done')
plt.show()