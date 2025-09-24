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

zcol = 'Z_not4clus'
#ram = 130
#rax = 220
ram = 0
rax = 360
ra0 = (ram+rax)/2.
decm = -0.5
decx = .5
zmin = 0
zmax = 3.5


#plt.figure()
fig, ax = plt.subplots(dpi=1000)
ax.set_aspect('equal')
ax.patch.set_facecolor('black')
#ax.patch.set_alpha(1)

msdic = {'QSO':.24,'ELG':.21,'LRG':.21,'BGS_ANY':.1}

tps = ['QSO','LRG','BGS_ANY','ELG']
cl = ['y','r','lime','b']
zordl = [2,5,3,1]
for tp,c,zo in zip(tps,cl,zordl):
    cols = ['RA','DEC',zcol,'ZWARN','DELTACHI2','LOCATION_ASSIGNED']
    if tp == 'ELG':
        cols.append('o2c')
        zmin = 0.6
    dt = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+tp+'_full.dat.fits',columns=cols)
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
    ax.plot(x,y,'s',color=c,zorder=zo,ms=msdic[tp],lw=0,mew=0,)
    if tp == 'QSO':
        sel = dt[zcol] > 2.1
        ax.plot(x[sel],y[sel],'s',color='white',zorder=zo,ms=msdic[tp],lw=0,mew=0)
    #plt.show()
    
    
    del dt


    print(tp+' done')
#plt.axis('off') 
for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(bottom=False, labelbottom=False,
               left=False, labelleft=False)   
plt.savefig('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/plots/wedge_all.png')
plt.show()
