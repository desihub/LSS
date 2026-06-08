#standard python
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
from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desitarget import targetmask
from desitarget.internal import sharedmem
from desimodel.footprint import is_point_in_desi
from desitarget import targetmask

import LSS.main.cattools as ct
import LSS.common_tools as common
import LSS.mocktools as mocktools
from LSS.globals import main

#import LSS.mkCat_singletile.fa4lsscat as fa
#from LSS.globals import main

if os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
#parser.add_argument("--tracer", help="tracer type to be selected")
parser.add_argument("--realization",type=int)
parser.add_argument("--prog", default="DARK")
parser.add_argument("--veto",default='_noveto')
#parser.add_argument("--mockdir", help="directory when pota mock data is",default='/global/cfs/cdirs/desi/users/acarnero/y1mock/SecondGen/clustering/')
parser.add_argument("--base_dir", help="base directory for input/output",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/')
parser.add_argument("--random_dir",help="where to find the data randoms",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v0.6/')
parser.add_argument("--specdata_dir",help="where to find the spec data ",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/')
parser.add_argument("--minr", help="minimum number for random files",default=0,type=int)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 40 are available (use parallel script for all)",default=1,type=int) 
parser.add_argument("--mockver", default=None, help = "which mocks to use. use abacus2ffa for Abacus 2nd gen fast fiber assignment")
parser.add_argument("--alt_indir", default = None)
parser.add_argument("--alt_outdir", default = None)
parser.add_argument("--tracer", default = None)
parser.add_argument("--remove_unassigned", default = 'y', help = 'set to n if dont want to include unassigned targets in catalog')


args = parser.parse_args()
print(args)

rm = int(args.minr)
rx = int(args.maxr)

notqso = ''
#if args.notqso == 'y':
#    notqso = 'notqso'

#tracer = args.tracer

# if tracer == 'LRG':
#     zmin = 0.4
#     zmax = 1.1
# 
# elif tracer == 'ELG_LOP':
#     zmin = 0.8
#     zmax = 1.6
# 
# elif tracer == 'QSO':
#     zmin = 0.8
#     zmax = 2.1
# 
# else:
#     sys.exit('tracer type '+args.tracer+' not supported (yet)')

if args.mockver != 'abacus2ffa':
    mockdir = args.base_dir+'mock'+str(args.realization)+'/'
    in_data_fn = mockdir+'pota-'+args.prog+'.fits'
    print(in_data_fn)
    cols = ['LOCATION',
    'FIBER',
    'TARGETID',
    'RA',
    'DEC','RSDZ',
    'PRIORITY_INIT',
    'PRIORITY',
    'DESI_TARGET','BRICKID','NOBS_G',
    'NOBS_R',
    'NOBS_Z',
    'MASKBITS','ZWARN',
    'COLLISION',
    'TILEID']
    mock_data = fitsio.read(in_data_fn,columns=cols)
    selcoll = mock_data['COLLISION'] == False
    mock_data = mock_data[selcoll]

elif args.mockver == 'abacus2ffa':
    mockdir = args.alt_indir
    in_data_fn = mockdir + args.tracer + '_joined_Ab2_m' + str(args.realization) + '.fits'
    mock_data = fitsio.read(in_data_fn)
    if args.remove_unassigned == 'y':
        mock_data = mock_data[mock_data["WEIGHT_IIP"] != 1e+20]

if args.prog == 'DARK':
    #bit = targetmask.desi_mask[args.tracer]
    bittest = targetmask.desi_mask
    desitarg='DESI_TARGET'
    tracers = ['LRG','QSO','ELG_LOP']
    if args.mockver == 'abacus2ffa':
        tracers = [args.tracer]

ndattot = len(mock_data)

lmockdat_noveto = len(mock_data)
mainp = main('LRG','iron','Y1') #needed for bad fiber list

if 'gtl' in args.veto:
    tilelocid = 10000*mock_data['TILEID']+mock_data['LOCATION']
    specfo = args.specdata_dir+'datcomb_'+args.prog.lower()+'_spec_zdone.fits'
    print('loading specf file '+specfo)
    specf = Table(fitsio.read(specfo))
    print(len(np.unique(specf['TILEID'])))
    specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    print('loaded specf file '+specfo)
    specfc = common.cut_specdat(specf,badfib=mainp.badfib)
    gtl = np.unique(specfc['TILELOCID'])
    goodtl = np.isin(tilelocid,gtl)
    mock_data = mock_data[goodtl]
    print(lmockdat_noveto,len(mock_data))

def ran_col_assign(randoms,data,sample_columns):
    inds = np.random.choice(len(data),len(randoms))
    dshuf = data[inds]
    for col in sample_columns:
        randoms[col] =  dshuf[col]
    return randoms

def apply_imaging_veto(ff,reccircmasks,ebits):
    if reccircmasks is not None:
        for maskfn in reccircmasks:
            mask = common.maskcircandrec(ff,maskfn)
            ff = ff[~mask]

    if ebits is not None:
        print('number before imaging mask '+str(len(ff)))
        if ebits == 'lrg_mask':
            sel = ff['lrg_mask'] == 0
            ff = ff[sel]
        else:
            ff = common.cutphotmask(ff,ebits)
        print('number after imaging mask '+str(len(ff)))
    return ff

    
for tracer in tracers:
    out_data_fn = mockdir+tracer+'_complete'+args.veto+'_clustering.dat.fits'
    out_data_froot = mockdir+tracer+'_complete'+args.veto+'_'

    mainp = main(tracer,'iron','Y1')
    bit = bittest[tracer]#targetmask.desi_mask[tracer]
    seltar = mock_data[desitarg] & bit > 0
    mock_data_tr = mock_data[seltar]
    lmockdat_noveto = len(mock_data_tr)
    print('length before/after cut to target type '+tracer)
    print(ndattot,len(mock_data))
    if tracer == 'LRG':
        zmin = 0.4
        zmax = 1.1

    elif (tracer == 'ELG_LOP') or (tracer == 'ELG'):
        zmin = 0.8
        zmax = 1.6

    elif tracer == 'QSO':
        zmin = 0.8
        zmax = 2.1


    
    mock_data_tr = Table(mock_data_tr)
    mock_data_tr = unique(mock_data_tr,keys=['TARGETID'])
    print('length after cutting to unique targetid',len(mock_data_tr))
    mock_data_tr.rename_column('RSDZ', 'Z')
    if args.mockver != 'abacus2ffa':
        mock_data_tr['WEIGHT'] = 1
    elif args.mockver == 'abacus2ffa':
        mock_data_tr['WEIGHT'] = mock_data_tr['WEIGHT_IIP']
    if 'imaging' in args.veto:
        if tracer == 'LRG':
            lrgmask = fitsio.read(args.base_dir+'forFA'+str(args.realization)+'_matched_input_full_lrg_imask.fits')
            mock_data_tr = join(mock_data_tr,lrgmask,keys=['TARGETID'])
            print(len(mock_data_tr))
        ebits = mainp.ebits
        reccircmasks = mainp.reccircmasks
        mock_data_tr = apply_imaging_veto(mock_data_tr,reccircmasks,ebits)
    selz = mock_data_tr['Z'] > zmin
    selz &= mock_data_tr['Z'] < zmax
    mock_data_tr = mock_data_tr[selz]
    print('length after cutting to redshift range',len(mock_data_tr))
    common.write_LSS(mock_data_tr,out_data_fn)

    def splitGC(flroot,datran='.dat',rann=0):
        import LSS.common_tools as common
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        app = 'clustering'+datran+'.fits'
        if datran == '.ran':
            app = str(rann)+'_clustering'+datran+'.fits'

        fn = Table(fitsio.read(flroot+app))
        c = SkyCoord(fn['RA']* u.deg,fn['DEC']* u.deg,frame='icrs')
        gc = c.transform_to('galactic')
        sel_ngc = gc.b > 0
        outf_ngc = flroot+'NGC_'+app
        common.write_LSS(fn[sel_ngc],outf_ngc)
        outf_sgc = flroot+'SGC_'+app
        common.write_LSS(fn[~sel_ngc],outf_sgc)

    splitGC(out_data_froot,'.dat')

    ran_samp_cols = ['Z','WEIGHT']


    for rann in range(rm,rx):
        tracerr = tracer
        if tracer == 'ELG_LOP':
            tracerr += 'notqso'
        if tracer == 'ELG':
            tracerr = 'ELG_LOPnotqso'
        in_ran_fn = args.random_dir+tracerr+'_'+str(rann)+'_full_noveto.ran.fits' #all noveto have same ra,dec, tracer becomes important for LRG imaging veto
        out_ran_fn = out_data_froot+str(rann)+'_clustering.ran.fits'
        rcols = ['RA','DEC','TILELOCID','NOBS_G','NOBS_R','NOBS_Z','MASKBITS']
        if tracerr == 'LRG':
            rcols.append('lrg_mask')
        ran = Table(fitsio.read(in_ran_fn,columns=rcols))
        if 'gtl' in args.veto:
            goodtl = np.isin(ran['TILELOCID'],gtl)
            ran = ran[goodtl]
        if 'imaging' in args.veto:
            ebits = mainp.ebits
            reccircmasks = mainp.reccircmasks
            ran = apply_imaging_veto(ran,reccircmasks,ebits)

        ran = ran_col_assign(ran,mock_data_tr,ran_samp_cols)
        common.write_LSS(ran,out_ran_fn)
        splitGC(out_data_froot,'.ran',rann)




