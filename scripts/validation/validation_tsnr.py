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
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
args = parser.parse_args()


indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
outdir = indir+'plots/tsnr/'

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
    

    df = fitsio.read(indir+tp+zdw+'_full.dat.fits')


    #dcl = []
    #for reg in regl:
    #    dtc = fitsio.read(indir+tp+zdw+reg+'_clustering.dat.fits')
    #    dcl.append(dtc)
    #dc = np.concatenate(dcl)
    #df = join(dtf,dc,keys=['TARGETID'],join_type='left')
    selgz = common.goodz_infull(tp[:3],df)#df['Z'].mask == False
    selo = df['ZWARN'] != 999999
    mean_gz = sum(df[selgz]['WEIGHT_ZFAIL'])/len(df[selo])
    print('number with good z, sum of weight_zfail,  number with good obs')
    print(len(df[selgz]),sum(df[selgz]['WEIGHT_ZFAIL']),len(df[selo]))
    tsnr = df['TSNR2_'+tp[:3]]
    for pt in range(0,10):
        fmin = pt*500
        fmax = fmin+500
        sel_fib = df['FIBER'] >= fmin
        sel_fib &= df['FIBER'] < fmax
        
        
        n_obs,rbins = np.histogram(tsnr[sel_fib&selo],bins=10,range=(np.percentile(tsnr[selo],1),np.percentile(tsnr[selo],99)))
        n_goodw,_ = np.histogram(tsnr[sel_fib&selgz],bins=rbins,weights=df[sel_fib&selgz]['WEIGHT_ZFAIL'])
        n_good,_ = np.histogram(tsnr[sel_fib&selgz],bins=rbins)
        bs = rbins[1]-rbins[0]
        rl = rbins[:-1]+bs/2.
        err = np.sqrt(n_goodw*(1.-n_good/n_obs))/n_obs
        plt.errorbar(rl,n_goodw/n_obs,err,fmt='ok',label='with WEIGHT_ZFAIL')
        plt.plot(rl,n_good/n_obs,'--k',label='no weight')
        plt.xlabel('TSNR2_'+tp[:3])
        plt.ylabel('fraction with good z in clus cat')
        plt.title(tp+' on petal '+str(pt))
        plt.grid()
        plt.legend()
        plt.savefig(outdir+tp+'_'+str(pt)+'_relsuccess_tnsr.png')
        plt.clf()
    n_obs,rbins = np.histogram(tsnr[selo],bins=20,range=(np.percentile(tsnr[selo],1),np.percentile(tsnr[selo],99)))
    n_goodw,_ = np.histogram(tsnr[selgz],bins=rbins,weights=df[selgz]['WEIGHT_ZFAIL'])
    n_good,_ = np.histogram(tsnr[selgz],bins=rbins)
    bs = rbins[1]-rbins[0]
    rl = rbins[:-1]+bs/2.
    err = np.sqrt(n_goodw*(1.-n_good/n_obs))/n_obs
    plt.errorbar(rl,n_goodw/n_obs,err,fmt='ok',label='with WEIGHT_ZFAIL')
    plt.plot(rl,n_good/n_obs,'--k',label='no weight')
    plt.xlabel('TSNR2_'+tp[:3])
    plt.ylabel('fraction with good z in clus cat')
    plt.title(tp+' on all petals ')
    plt.grid()
    plt.legend()
    plt.savefig(outdir+tp+'_allpetals_relsuccess_tnsr.png')
    plt.clf()



