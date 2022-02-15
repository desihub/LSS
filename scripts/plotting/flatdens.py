import matplotlib.pyplot as plt
import numpy as np
import os
import sys

import fitsio
from astropy.table import join,Table


outdir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/plots/'
qt = 'COMP_TILE'

nside = 256
nest = True

tps = ['QSO','LRG','BGS_ANY','BGS_BRIGHT','ELG','ELG_LOP','ELG_LOPnotqso']
for tp in tps:
    rf = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+tp+'zdone_1_full.ran.fits'
    dt = Table.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+tp+'zdone_full.dat.fits')
    wz = dt['ZWARN']*0 == 0
    wz &= dt['ZWARN'] != 1.e20
    wz &= dt['ZWARN'] != 999999
    wz &= dt['LOCATION_ASSIGNED'] == 1
    dt = dt[wz]
    if tp == 'QSO':
        #good redshifts are currently just the ones that should have been defined in the QSO file when merged in full
        wg = dt[zcol]*0 == 0
        wg &= dt[zcol] != 999999
        wg &= dt[zcol] != 1.e20
        wg &= dt[zcol] > 0.8
        wg &= dt[zcol] < 2.1
    
    if tp[:3] == 'ELG':
        wg = dt['o2c'] > 0.9
        wg &= dt[zcol] > 0.8
        wg &= dt[zcol] < 1.6

    if tp == 'LRG':
        # Custom DELTACHI2 vs z cut from Rongpu
        wg = dt['ZWARN'] == 0
        drz = (10**(3 - 3.5*dt[zcol]))
        mask_bad = (drz>30) & (dt['DELTACHI2']<30)
        mask_bad |= (drz<30) & (dt['DELTACHI2']<drz)
        mask_bad |= (dt['DELTACHI2']<10)
        wg &= dt[zcol]<1.4
        wg &= (~mask_bad)
        wg &= dt[zcol] > 0.4
        wg &= dt[zcol] < 1.1


    if tp[:3] == 'BGS':
        wg = dt['DELTACHI2'] > 40
        wg &= dt[zcol] > 0.1
        wg &= dt[zcol] < .5

    dt = dt[wg]
    dt.keep_columns(['RA','DEC',zcol,'FRACZ_TILELOCID'])
    dt['WEIGHT'] = 1/dt['FRACZ_TILELOCID']
    rt = fitsio.read(rf,columns=['RA','DEC'])
    print(tp)
    wp,od = densvar.get_hpdens(rt,dt,datweights='WEIGHT',sz=.2,vm=.8,vx=1.2)

    pixls = np.arange(12*nside*nside,dtype=int)
    th,phi = hp.pix2ang(nside,pixls[wp],nest=nest)
    ra,dec = densvar.thphi2radec(th,phi)
    
    wr = ra > 300
    ra[wr] -=360
    vs = 1.2
    vm = 0.8

    plt.scatter(ra,np.sin(dec*np.pi/180),c=od,s=.1,edgecolor='none',vmax=vx,vmin=vm)
    plt.xlabel('RA')
    plt.ylabel('sin(DEC)')
    plt.colorbar()
    plt.title(tp +' completenss weighted over-density')


    plt.savefig(outdir+tp+'_compweighteddens.png')
    plt.clf()
    del dt
    del rt
    print(tp+' done')