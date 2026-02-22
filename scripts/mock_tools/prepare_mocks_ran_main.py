from astropy.io import fits # Access to FITS (Flexible Image Transport System) files.
from astropy.table import Table, hstack, vstack, Column # A class to represent tables of heterogeneous data.
import fitsio
import numpy as np
import glob
import os
import h5py
import argparse
import sys

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--mockver", help="type of mock to use",default='ab_firstgen')
parser.add_argument("--ranmin", help="number for the realization",default=1,type=int)
parser.add_argument("--ranmax", help="number for the realization",default=2,type=int)
parser.add_argument("--prog", help="dark or bright",default='dark')
parser.add_argument("--footprint", help="e.g., DA02 or SV3",default='DA02')
parser.add_argument("--base_output", help="base directory for output",default='/global/cfs/cdirs/desi/survey/catalogs/main/mocks/')
parser.add_argument("--nproc", help="run different random number in parallel?",default=1)
parser.add_argument("--prep", help="prepare file for fiberassign?",default='y')
parser.add_argument("--runfa", help="run fiberassign",default='y')
parser.add_argument("--stepfa", help="only run fiberassign in runfa",default='n')


args = parser.parse_args()


def prep(rannum):

    if args.mockver == 'ab_firstgen':
        mockpath = '/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/'
    
        snum = str(100*rannum)
        file_name = mockpath+'ELG/z1.100/cutsky_ELG_random_S'+snum+'_1X.fits'
        out_file_name = args.base_output+'/FirstGenMocks/AbacusSummit/ran_forFA'+str(rannum)+'.fits'
        print('will write to '+out_file_name)
        if not os.path.exists(args.base_output+'/FirstGenMocks'):
            os.mkdir(args.base_output+'/FirstGenMocks')
            print('made '+args.base_output+'/FirstGenMocks')
        if not os.path.exists(args.base_output+'/FirstGenMocks/AbacusSummit'):
            os.mkdir(args.base_output+'/FirstGenMocks/AbacusSummit')
            print('made '+args.base_output+'/FirstGenMocks/AbacusSummit')
        mockdir = args.base_output+'/FirstGenMocks/AbacusSummit'
 

        #def mask(main=0, nz=0, Y5=0, sv3=0):
        #   return main * (2**3) + sv3 * (2**2) + Y5 * (2**1) + nz * (2**0)
        
        if args.prep == 'y':
            targets = fitsio.read(file_name,columns=['RA','DEC','STATUS'])
            sel = targets['STATUS'] & 2 > 0 #cut to Y5 footprint
            data = targets[sel]
            targets = Table(targets)
            targets.remove_columns(['STATUS'])

    elif args.mockver == 'EZ_3gpc1year':
        mockdir = args.base_output+'/FA_EZ_1year/fiberassign_EZ_3gpc'
    else:
        sys.exit(args.mockver+' not supported')

    if args.prep == 'y':
        n=len(targets)
        targets['DESI_TARGET'] = np.ones(n, dtype='i8')
        targets['SUBPRIORITY'] = np.random.uniform(0, 1, n)
        targets['OBSCONDITIONS'] = np.zeros(n, dtype='i8')+int(3) 
        targets['NUMOBS_MORE'] = np.zeros(n, dtype='i8')+int(1) 
        targets['NUMOBS_INIT'] = np.zeros(n, dtype='i8')+int(1)
        targets['ZWARN'] = np.zeros(n, dtype='i8')+int(0)
        targets['TARGETID'] = np.arange(1,n+1)+10*n*rannum #each random file has approximately the same number, so this should keep the targetid distinct

        targets.write(out_file_name, overwrite = True)

        fits.setval(out_file_name, 'EXTNAME', value='TARGETS', ext=1)
        fits.setval(out_file_name, 'OBSCON', value=args.prog.upper(), ext=1)
        
    if args.runfa == 'y':
        #from LSS.mocktools import get_fba_mock_ran
        #get_fba_mock_ran(mockdir,rannum,survey=args.survey,prog=args.prog)
        targfn = out_file_name
        tile_fn = '/global/cfs/cdirs/desi/survey/catalogs/'+args.footprint+'/LSS/tiles-'+args.prog.upper()+'.fits'
        tiles = fitsio.read(tile_fn)
        ts = str(tiles['TILEID'][0]).zfill(6)
        #get info from origin fiberassign file
        fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
        rundate= fht['RUNDATE']
        rootdir = mockdir+'/'+args.footprint+'/ran'+str(rannum)+'_'+args.prog+'/'
        if not os.path.exists(rootdir):
            os.mkdir(rootdir)
            print('made '+rootdir)
        runstr = 'python fa_ran.py --infn '+targfn+' --outdir '+rootdir+' --program '+args.prog+' --rundate '+rundate +' --tilesfn '+tile_fn +' --numproc '+str(args.nproc)
        if args.stepfa == 'y':
            runstr += ' --steps fa'
        os.system(runstr)



if __name__ == '__main__':
    rx = args.ranmax
    rm = args.ranmin
#     if args.par == 'y':
#         from multiprocessing import Pool
#         from desitarget.internal import sharedmem
#         
#         N = rx-rm+1
#         inds = []
#         for i in range(rm,rx):
#             inds.append(i)
#         pool = sharedmem.MapReduce(np=N)
#         with pool:
#         
#             def reduce( r):
#                 print('chunk done')
#                 return r
#             pool.map(prep,inds,reduce=reduce)
# 
#         #p.map(doran,inds)
#     else:
    for i in range(rm,rx):
        prep(i)

