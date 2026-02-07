import numpy as np
import fitsio
from astropy.table import Table, join, unique,vstack
import os
import sys
from matplotlib import pyplot as plt
import argparse

import LSS.common_tools as common
from LSS.globals import main


outdir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/dNdzmonth/'
mtld = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
sel = mtld['SURVEY'] == 'main'
mtld = mtld[sel]
yms = np.unique(mtld['LASTNIGHT']//100)
print('months to go through are:')
print(yms)

def dndz_monthall(yearmonths,tp,zcol='Z_not4clus'):

    mainp = main(tp,'daily')
    if tp == 'QSO':
        dt = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/datcomb_QSO_tarspecwdup_zdone.fits')
        #print(dt.dtype.names)
        dt = common.cut_specdat(dt,mainp.badfib)
        sel = dt['PRIORITY'] == 3400 #select QSO on 1st obs
        dt = dt[sel]
        azf ='/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/QSO_catalog.fits'
        arz = Table(fitsio.read(azf))
        arz.keep_columns(['TARGETID','LOCATION','TILEID','Z','Z_QN'])
        arz['TILEID'] = arz['TILEID'].astype(int)

        dt = join(dt,arz,keys=['TARGETID','TILEID','LOCATION'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
        #print(dt.dtype.names)
        zcol = 'Z'
        
    elif tp != 'ELGnotqso' and tp != 'ELGandQSO':
        dt = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/'+tp+'_full.dat.fits')
        #
    else:
        dt = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/ELG_full.dat.fits')
        #dt = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/datcomb_'+type+'_tarspecwdup_zdone.fits')

    wg = np.ones(len(dt),dtype='bool')
    if tp != 'QSO':
        wg &= dt['ZWARN'] != 999999
        wg &= dt['ZWARN'] != 1e20
        wg &= dt['ZWARN']*0 == 0
        wg &= dt['GOODHARDLOC'] == 1
    
    zmin = 0
    zmax = 2

    if tp == 'LRG' or tp == 'LGE':
        #drz = (10**(3 - 3.5*dt[zcol]))
        #mask_bad = (drz>30) & (dt['DELTACHI2']<30)
        #mask_bad |= (drz<30) & (dt['DELTACHI2']<drz)
        #mask_bad |= (dt['DELTACHI2']<10)
        wz = dt['DELTACHI2'] > 15
        wz &= dt['ZWARN'] == 0
        wz &= dt[zcol]<1.5
        #wz &= (~mask_bad)    
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
        #wz = dt['SPECTYPE'] == 'QSO'
        wz = dt[zcol]*0 == 0
        wz &= dt[zcol] != 999999
        zmin = 0
        zmax = 4.5
    if tp[:3] == 'BGS':
        wz = dt['DELTACHI2'] > 40
        zmin = 0
        zmax = 1
    zl = dt[wg&wz][zcol]
    #wl = 1./dt[wg&wz]['FRACZ_TILELOCID']
    fractot = len(dt[wg&wz])/len(dt[wg])
    for yearmonth in yearmonths:
        sel = mtld['LASTNIGHT']//100 == yearmonth
        #print(len(mtld[sel]))
        tids = np.unique(mtld[sel]['TILEID'])
        sd = np.isin(dt['TILEID'],tids)
        ntls = len(np.unique(dt[sd]['TILEID']))
        zlm = dt[wg&wz&sd][zcol]
        #wlm = 1./dt[wg&wz&sd]['FRACZ_TILELOCID']
        if len(dt[wg&sd]) > 0:
            fracm = len(dt[wg&wz&sd])/len(dt[wg&sd])
            print(tp,yearmonth,'success rate '+str(fracm))
            #plt.hist(zl,bins=50,density=True,weights=wl,histtype='step',label='all; ssr '+str(round(fractot,3)),range=(zmin,zmax))
            #plt.hist(zlm,bins=50,density=True,weights=wlm,histtype='step',label=str(yearmonth)+', '+str(ntls)+' tiles; ssr '+str(round(fracm,3)),range=(zmin,zmax))
            plt.hist(zl,bins=50,density=True,histtype='step',label='all; ssr '+str(round(fractot,3)),range=(zmin,zmax))
            plt.hist(zlm,bins=50,density=True,histtype='step',label=str(yearmonth)+', '+str(ntls)+' tiles; ssr '+str(round(fracm,3)),range=(zmin,zmax))
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
            #del wlm
            plt.clf()
    del dt

parser = argparse.ArgumentParser()
parser.add_argument("--tracer",help="tracer type to make plots for ",default='all')
args = parser.parse_args()


if args.tracer == 'all':
    tps = ['LGE','ELGnotqso','ELG_LOPnotqso','ELGandQSO','BGS_ANY','BGS_BRIGHT','LRG','QSO']
else:
    tps = [args.tracer]
for tp in tps:
	dndz_monthall(yms,tp)
