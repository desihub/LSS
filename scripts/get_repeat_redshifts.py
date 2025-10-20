import sys
import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import fitsio
import glob
import argparse
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt

from LSS import common_tools as common
from LSS.main import cattools as ct
from LSS.globals import main

parser = argparse.ArgumentParser()
parser.add_argument("--survey",help="set of tiles, e.g., Y1 for DR1, DA2 for DR2, main for all of main survey",default='DA2')
parser.add_argument("--verspec",help="version for redshifts",default='loa-v1')
parser.add_argument("--outdir",help="directory for output",default=os.getenv('SCRATCH'))
args = parser.parse_args()

def get_repeats(indat,zcol='Z'):
    tids,cnts = np.unique(indat['TARGETID'],return_counts=True)
    rtids = tids[cnts>1]
    sel_r = np.isin(indat['TARGETID'],rtids)
    specflr = indat[sel_r]
    print(len(specflr))
    specflr.sort('TARGETID')
    dzl = []
    zl1 = []
    zl2 = []
    tidl = []
    ind = 0
    while ind < len(specflr):
        tid = specflr[ind]['TARGETID']
        z1 = specflr[ind][zcol]
        ind2 = 1
        #while (ind+ind2) < len(specf_LRGr):
        while specflr[ind+ind2]['TARGETID'] == tid:
            zn = specflr[ind+ind2][zcol]
            dzl.append((z1-zn)/(1+z1))
            zl1.append(z1)
            zl2.append(zn)
            tidl.append(tid)
            ind2 += 1
            if ind+ind2 >= len(specflr):
                break
        ind += ind2
        if ind%1e4 == 0:
            print(ind)    
    res = Table()
    res['TARGETID'] = tidl
    res['Z1'] = zl1
    res['Z2'] = zl2
    return res
    
#do BGS

mainp = main('BGS',args.verspec)#get settings for dark time

mt = mainp.mtld
tiles = mainp.tiles
imbits = mainp.imbits #mask bits applied to targeting
ebits = mainp.ebits #extra mask bits we think should be applied


tsnrcut = mainp.tsnrcut
dchi2 = mainp.dchi2
tnsrcol = mainp.tsnrcol        
zmin = mainp.zmin
zmax = mainp.zmax
badfib = mainp.badfib

#get set of tiles
wd = mt['SURVEY'] == 'main'
wd &= mt['ZDONE'] == 'true'
wd &= mt['FAPRGRM'] == 'bright'
if args.survey == 'Y1':
    wd &= mt['ZDATE'] < 20220900
if args.survey == 'DA2':
    wd &= mt['ZDATE'] < 20240410
mtld = mt[wd]

ldirspec = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/'+args.verspec+'/'
specfo = ldirspec+'datcomb_bright_spec_zdone.fits'
specf = Table(fitsio.read(specfo.replace('global','dvs_ro')))
sel = np.isin(specf['TILEID'],mtld['TILEID'])
specf = specf[sel]
specf = Table(specf)
specf.keep_columns(['TARGETID','Z','ZWARN','DELTACHI2','LOCATION','BGS_TARGET','TILEID','TSNR2_LRG','TSNR2_ELG','ZWARN_MTL','COADD_FIBERSTATUS','FIBER','LASTNIGHT'])
specf = common.cut_specdat(specf,badfib=mainp.badfib_td,tsnr_min=tsnrcut,tsnr_col=tnsrcol,fibstatusbits=mainp.badfib_status,remove_badfiber_spike_nz=True,mask_petal_nights=True)

#do BGS
sel_BGS = specf['BGS_TARGET'] > 0
sel_gz = common.goodz_infull('BGS',specf,zcol='Z')
specfl = specf[sel_BGS&sel_gz]

bgsr = get_repeats(specfl)
bgsr.write(args.outdir+'/BGSrepeats.fits',overwrite=True)

sel = abs((bgsr['Z1']-bgsr['Z2'])/(1+bgsr['Z1'])) > 0.003
print('fraction of repeat BGS measurements with (Z1-Z2)/(1+Z1) > 0.003:')
print(np.sum(sel)/len(bgsr))
    
#do dark time tracers

mainp = main('LRG',args.verspec)#get settings for dark time

mt = mainp.mtld
tiles = mainp.tiles
imbits = mainp.imbits #mask bits applied to targeting
ebits = mainp.ebits #extra mask bits we think should be applied


tsnrcut = mainp.tsnrcut
dchi2 = mainp.dchi2
tnsrcol = mainp.tsnrcol        
zmin = mainp.zmin
zmax = mainp.zmax
badfib = mainp.badfib

#get set of tiles
wd = mt['SURVEY'] == 'main'
wd &= mt['ZDONE'] == 'true'
wd &= mt['FAPRGRM'] == 'dark'
if args.survey == 'Y1':
    wd &= mt['ZDATE'] < 20220900
if args.survey == 'DA2':
    wd &= mt['ZDATE'] < 20240410
mtld = mt[wd]

ldirspec = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/'+args.specver+'/'
specfo = ldirspec+'datcomb_dark_spec_zdone.fits'
specf = Table(fitsio.read(specfo.replace('global','dvs_ro')))
sel = np.isin(specf['TILEID'],mtld['TILEID'])
specf = specf[sel]
specf = Table(specf)
specf.keep_columns(['TARGETID','Z','ZWARN','DELTACHI2','LOCATION','DESI_TARGET','TILEID','TSNR2_LRG','TSNR2_ELG','ZWARN_MTL','COADD_FIBERSTATUS','FIBER','LASTNIGHT'])
specf = common.cut_specdat(specf,badfib=mainp.badfib_td,tsnr_min=tsnrcut,tsnr_col=tnsrcol,fibstatusbits=mainp.badfib_status,remove_badfiber_spike_nz=True,mask_petal_nights=True)

#do LRGs
sel_LRG = (specf['DESI_TARGET']) & 1 > 0
sel_gz = common.goodz_infull('LRG',specf,zcol='Z')
specfl = specf[sel_LRG&sel_gz]

lrgr = get_repeats(specfl)
lrgr.write(args.outdir+'/LRGrepeats.fits',overwrite=True)

sel = abs((lrgr['Z1']-lrgr['Z2'])/(1+lrgr['Z1'])) > 0.003
print('fraction of repeat LRG measurements with (Z1-Z2)/(1+Z1) > 0.003:')
print(np.sum(sel)/len(lrgr))

#do QSO
qsof = fitsio.read(mainp.qsozf,columns=['TARGETID','LOCATION','TILEID','Z'])
specf = join(specf,qsof,keys=['TARGETID','TILEID','LOCATION'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
selqso = specf['Z_QF']!=999999
selqso &= (specf['DESI_TARGET'] & 4) > 0
specfq = specf[selqso]
qsor = get_repeats(specfq,zcol='Z_QF')
qsor.write(args.outdir+'/QSOrepeats.fits',overwrite=True)
sel = abs((qsor['Z1']-qsor['Z2'])/(1+qsor['Z1'])) > 0.01
print('fraction of repeat QSO measurements with (Z1-Z2)/(1+Z1) > 0.01:')
print(np.sum(sel)/len(qsor))

#do elg
elgf = fitsio.read(mainp.elgzf,columns=['TARGETID','LOCATION','TILEID','OII_FLUX','OII_FLUX_IVAR'])
specf = join(specf,elgf,keys=['TARGETID','TILEID','LOCATION'],join_type='left')
o2c = np.log10(specf['OII_FLUX'] * np.sqrt(specf['OII_FLUX_IVAR']))+0.2*np.log10(specf['DELTACHI2'])
w = (o2c*0) != 0
w |= specf['OII_FLUX'] < 0
o2c[w] = -20
specf['o2c'] = o2c
sel_ELG = (specf['DESI_TARGET'] & 2) > 0
sel_gz = common.goodz_infull('ELG',specf,zcol='Z')
specfe = specf[sel_ELG&sel_gz]

elgr = get_repeats(specfe)
elgr.write(args.outdir+'/QSOrepeats.fits',overwrite=True)
sel = abs((elgr['Z1']-elgr['Z2'])/(1+elgr['Z1'])) > 0.001
print('fraction of repeat ELG measurements with (Z1-Z2)/(1+Z1) > 0.001:')
print(np.sum(sel)/len(elgr))


