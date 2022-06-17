import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import argparse

import fitsio
from astropy.table import join,Table
import healpy as hp

from LSS.imaging import densvar

parser = argparse.ArgumentParser()
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='SV3')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='fuji')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
args = parser.parse_args()


indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
outdir = indir+'plots/'

if args.data == 'LSS':
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print('made '+outdir)


zcol = 'Z'
nran = 18

tps = [args.tracers]
if args.tracers == 'all':
    tps = ['QSO','LRG','BGS_BRIGHT','ELG_LOPnotqso']

zdw = ''#'zdone'

regl = ['_N','_S']

if args.survey == 'SV3' and args.tracers == 'all':
    tps = ['QSO','LRG','BGS_ANY','BGS_BRIGHT','ELG','ELG_HIP','ELG_HIPnotqso','ELGnotqso']
    zdw = ''
    if args.data != 'LSS':
        tps = ['QSO','LRG','ELG']
for tp in tps:
    

    dtf = fitsio.read(indir+tp+zdw+'_full.dat.fits')


    dcl = []
    for reg in regl:
        dtc = fitsio.read(indir+tp+zdw+reg+'_clustering.dat.fits')
        dcl.append(dtc)
    dc = np.concatenate(dcl)
    df = join(dtf,dc,keys=['TARGETID'],join_type='left')
    selgz = df['Z'].mask == False
    selo = df['ZWARN'] != 999999
    mean_gz = sum(df[selgz]['WEIGHT_ZFAIL'])/len(df[selo])
    print('number with good z, sum of weight_zfail,  number with good obs')
    print(len(df[selgz]),sum(df[selgz]['WEIGHT_ZFAIL']),len(df[selo]))
    for pt in range(0,10):
        fmin = pt*500
        fmax = fmin+500
        sel_fib = df['FIBER'] >= fmin
        sel_fib &= df['FIBER'] < fmax
        fib_obs,n_obs = np.unique(df[sel_fib&selo]['FIBER'],return_counts=True)
        nw_goodz = []
        n_goodz = []
        dat_gz = df[selgz&sel_fib]
        for fib in fib_obs:
            sel_fn = dat_gz['FIBER'] == fib
            nw_goodz.append(sum(dat_gz[sel_fn]['WEIGHT_ZFAIL']))
            n_goodz.append(len(dat_gz[sel_fn]))
        nw_goodz = np.array(nw_goodz)
        n_goodz = np.array(n_goodz)
        n_fail = n_obs-n_goodz
        plt.errorbar(fib_obs,nw_goodz/n_obs,np.sqrt(n_fail)/n_obs,fmt='ko')
        plt.plot(fib_obs,np.ones(len(fib_obs))*mean_gz,'k--')
        plt.xlabel('FIBER')
        plt.ylabel('fraction with good z in clus cat')
        plt.title(tp+' on petal '+str(pt))
        plt.show()
        plt.errorbar(fib_obs,nw_goodz/n_obs,np.sqrt(n_fail)/n_obs,fmt='ko')
        plt.xlabel('FIBER')
        plt.ylabel('fraction with good z in clus cat')
        plt.title(tp+' on petal '+str(pt))
        plt.ylim(mean_gz-0.1,mean_gz+0.1)
        plt.show()
                    


