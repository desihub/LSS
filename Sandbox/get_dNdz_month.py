import numpy as np
import fitsio
from astropy.table import Table, join, unique,vstack
import os
import sys
from matplotlib import pyplot as plt

outdir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/dNdzmonth/'
mtld = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
sel = mtld['SURVEY'] == 'main'
yms = np.unique(mtld[sel]['LASTNIGHT']//100)
print('months to go through are:')
print(yms)

def dndz_monthall(yearmonths,tp,zcol='Z_not4clus'):
    
    if tp != 'ELGnotqso' and tp != 'ELGandQSO':
        dt = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+tp+'zdone_full.dat.fits')
    else:
        dt = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/ELGzdone_full.dat.fits')

    wg = dt['ZWARN'] != 999999
    wg &= dt['GOODHARDLOC'] == 1
    zmin = 0
    zmax = 2

    if tp == 'LRG':
        drz = (10**(3 - 3.5*dt[zcol]))
        mask_bad = (drz>30) & (dt['DELTACHI2']<30)
        mask_bad |= (drz<30) & (dt['DELTACHI2']<drz)
        mask_bad |= (dt['DELTACHI2']<10)
        wz = dt[zcol]<1.4
        wz &= (~mask_bad)    
    if tp[:3] == 'ELG':
        wg &= dt['o2c'] != 1e20
        wz = dt['o2c'] > 0.9
        if tp == 'ELGnotqso':
            wg &= (dt['DESI_TARGET'] & 4) == 0
        if tp == 'ELGandQSO':
            wg &= (dt['DESI_TARGET'] & 4) > 0
            wz |= dt['SPECTYPE'] == 'QSO'
            zmin = 0
            zmax = 4.5
    if tp == 'QSO':
        wz = dt['SPECTYPE'] == 'QSO'
        zmin = 0
        zmax = 4.5
    if tp[:3] == 'BGS':
        wz = dt['DELTACHI2'] > 40
        zmin = 0
        zmax = 1
    zl = dt[wg&wz][zcol]
    wl = 1./dt[wg&wz]['FRACZ_TILELOCID']
    fractot = len(dt[wg&wz])/len(dt[wg])
    for yearmonth in yearmonths:
        sel = mtld['LASTNIGHT']//100 == yearmonth
        #print(len(mtld[sel]))
        tids = np.unique(mtld[sel]['TILEID'])
        sd = np.isin(dt['TILEID'],tids)
        ntls = len(np.unique(dt[sd]['TILEID']))
        zlm = dt[wg&wz&sd][zcol]
        wlm = 1./dt[wg&wz&sd]['FRACZ_TILELOCID']
        fracm = len(dt[wg&wz&sd])/len(dt[wg&sd])
        plt.hist(zl,bins=50,density=True,weights=wl,histtype='step',label='all; ssr '+str(round(fractot,3)),range=(zmin,zmax))
        plt.hist(zlm,bins=50,density=True,weights=wlm,histtype='step',label=str(yearmonth)+', '+str(ntls)+' tiles; ssr '+str(round(fracm,3)),range=(zmin,zmax))
        plt.title(tp)
        plt.xlabel('Z')
        plt.ylabel('dN/dz')
        plt.legend()
        #if tp == 'LRG' or tp[:3] == 'ELG':
        #    plt.xlim(0,2)
        if tp == 'ELG_LOPnotqso' or tp == 'ELGnotqso': 
            plt.ylim(0,1.7)
        if tp == 'ELGandQSO': 
            plt.ylim(0,0.8)
        if tp == 'BGS_ANY':
            plt.ylim(0,3.8)
        if tp == 'BGS_BRIGHT':
            plt.ylim(0,4.1)
        if tp == 'LRG':
            plt.ylim(0,2.2)
        if tp == 'QSO':
            plt.ylim(0,0.7)
        plt.savefig(outdir+tp+str(yearmonth)+'.png')
        del zlm
        del wlm
        plt.clf()
    del dt

tps = ['LRG','QSO','ELGnotqso','ELG_LOPnotqso','ELGandQSO','BGS_ANY','BGS_BRIGHT']
for tp in tps:
	dndz_monthall(yms,tp)