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
parser.add_argument("--veto",default='_imaging')
#parser.add_argument("--mockdir", help="directory when pota mock data is",default='/global/cfs/cdirs/desi/users/acarnero/y1mock/SecondGen/clustering/')
parser.add_argument("--base_dir", help="base directory for input/output",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/')
parser.add_argument("--random_dir",help="where to find the data randoms",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v0.6/')
parser.add_argument("--specdata_dir",help="where to find the spec data ",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/')
parser.add_argument("--minr", help="minimum number for random files",default=0,type=int)
parser.add_argument("--maxr", help="maximum for random files, default is all 18)",default=18,type=int) 
parser.add_argument("--tracer", default = 'LRG')
parser.add_argument("--par", default = 'n',help='whether to run random steps in parallel or not')
parser.add_argument("--apply_HPmapcut", default = 'y')
#parser.add_argument("--remove_unassigned", default = 'y', help = 'set to n if dont want to include unassigned targets in catalog')


args = parser.parse_args()
print(args)

rm = int(args.minr)
rx = int(args.maxr)

notqso = ''



if args.tracer == 'all':
    tracers = ['LRG','QSO','ELG_LOP']
else:
    tracers = [args.tracer]
    #if args.mockver == 'abacus2ffa':
    #    tracers = [args.tracer]

#ndattot = len(mock_data)

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

nproc = 9
    
for tracer in tracers:

    out_data_froot = args.base_dir+tracer+'_ffa'+args.veto+'_'
    if args.apply_HPmapcut == 'y':
        out_data_froot += 'HPmapcut'

    mainp = main(tracer,'iron','Y1')

    nran = rx-rm
    #if args.mkran == 'y':
    def _mkran(rann):
        
        tracerr = tracer
        if tracer == 'ELG_LOP':
            tracerr += 'notqso'
        if tracer == 'ELG':
            tracerr = 'ELG_LOPnotqso'
        in_ran_fn = args.random_dir+tracerr+'_'+str(rann)+'_full_noveto.ran.fits' #all noveto have same ra,dec, tracer becomes important for LRG imaging veto
        out_ran_fn = out_data_froot+str(rann)+'_full.ran.fits'
        rcols = ['RA','DEC','TILELOCID','NOBS_G','NOBS_R','NOBS_Z','MASKBITS','TARGETID','NTILE','GOODHARDLOC','PHOTSYS']
        if tracerr == 'LRG':
            rcols.append('lrg_mask')
        ran = Table(fitsio.read(in_ran_fn,columns=rcols))
        ran['FRAC_TLOBS_TILES'] = 1.
    
        print(len(ran),' random length before good loc cut')
        goodtl = ran['GOODHARDLOC']#np.isin(ran['TILELOCID'],gtl)
        ran = ran[goodtl]
        print(len(ran),' random length after good loc cut')
        if 'imaging' in args.veto:
            ebits = mainp.ebits
            reccircmasks = mainp.reccircmasks
            ran = apply_imaging_veto(ran,reccircmasks,ebits)
        if args.apply_HPmapcut == 'y':
            nside = 256
            lssmapdirout = args.random_dir+'/hpmaps/'
            if tracer[:3] == 'BGS':
                mapn = fitsio.read(lssmapdirout+'BGS_BRIGHT_mapprops_healpix_nested_nside'+str(nside)+'_N.fits')
                maps = fitsio.read(lssmapdirout+'BGS_BRIGHT_mapprops_healpix_nested_nside'+str(nside)+'_S.fits')
            
            else:
                mapn = fitsio.read(lssmapdirout+'QSO_mapprops_healpix_nested_nside'+str(nside)+'_N.fits')
                maps = fitsio.read(lssmapdirout+'QSO_mapprops_healpix_nested_nside'+str(nside)+'_S.fits')
            mapcuts = mainp.mapcuts
            ran = common.apply_map_veto_arrays(ran,mapn,maps,mapcuts,nside)
            print('random veto '+str(rann)+' done')

        common.write_LSS(ran,out_ran_fn)

    inds = np.arange(rm,rx)
    if args.par == 'y':
        from multiprocessing import Pool
        with Pool(processes=nproc) as pool:
            res = pool.map(_mkran, inds)
    else:
        for rn in inds:#range(rm,rx):
             _mkran(rn)
    
    


    


