import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import sys
import argparse

import fitsio
from astropy.table import join,Table
import healpy as hp

from LSS.imaging import densvar
from LSS import common_tools as common

parser = argparse.ArgumentParser()
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--zmin",help="minimum redshift",default=None)
parser.add_argument("--zmax",help="maximum redshift",default=None)
args = parser.parse_args()


indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
outdir = indir+'plots/fiber/'

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
    
    figs = []
    dtf = fitsio.read(indir+tp+zdw+'_full.dat.fits')

    if args.survey == 'DA02' or args.survey == 'SV3':
        dcl = []
        for reg in regl:
            dtc = fitsio.read(indir+tp+zdw+reg+'_clustering.dat.fits')
            dcl.append(dtc)
        dc = np.concatenate(dcl)
        df = join(dtf,dc,keys=['TARGETID'],join_type='left')
        zcol = 'Z'
        selgz = df['Z'].mask == False
    else:
        df = dtf
        zcol = 'Z_not4clus'
        selgz = common.goodz_infull(tp[:3],df)
    zw = ''
    if args.zmin is not None:
        selgz &= df[zcol] > float(args.zmin)
        zw += str(args.zmin)
    if args.zmin is not None:
        selgz &= df[zcol] < float(args.zmax)
        zw += str(args.zmax)   
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
        err = np.sqrt(nw_goodz*n_fail/n_obs)/n_obs
        sel_4p = n_obs > 10
        fig = plt.figure()
        plt.errorbar(fib_obs[sel_4p],nw_goodz[sel_4p]/n_obs[sel_4p],err[sel_4p],fmt='.k')
        plt.plot(fib_obs,np.ones(len(fib_obs))*mean_gz,'r--')
        sel_3l = (mean_gz - (nw_goodz/n_obs))/err > 3
        for ii in range(0,len(fib_obs[sel_3l&sel_4p])):
            plt.text(fib_obs[sel_3l&sel_4p][ii],nw_goodz[sel_3l&sel_4p][ii]/n_obs[sel_3l&sel_4p][ii],str(fib_obs[sel_3l&sel_4p][ii]),color='red')
        plt.xlabel('FIBER')
        plt.ylabel('fraction with good z in clus cat')
        plt.title(tp+' on petal '+str(pt)+ ' '+zw)
        plt.grid()
        figs.append(fig)
        #plt.savefig(outdir+tp+'_'+str(pt)+zw+'_zfailweighted_relsuccess_fiber.png')
        #plt.clf()
        fig = plt.figure()
        plt.errorbar(fib_obs[sel_4p],nw_goodz[sel_4p]/n_obs[sel_4p],err[sel_4p],fmt='.k')
        plt.plot(fib_obs,np.ones(len(fib_obs))*mean_gz,'r--')
        plt.xlabel('FIBER')
        plt.ylabel('fraction with good z in clus cat')
        plt.title(tp+' on petal '+str(pt))
        plt.ylim(mean_gz-0.1,mean_gz+0.1)
        plt.grid()
        figs.append(fig)
        #plt.savefig(outdir+tp+'_'+str(pt)+zw+'_zfailweighted_relsuccess_fiber_zoom.png')
        #plt.clf()
    fibrel = df['FIBER'] -500*(df['FIBER']//500)
    fib_obs,n_obs = np.unique(fibrel[selo],return_counts=True)                
    sel_4p = n_obs > 10
    fibrel_gz = fibrel[selgz]
    dat_gz = df[selgz]
    nw_goodz = []
    n_goodz = []

    for fib in fib_obs:
        sel_fn = fibrel_gz == fib
        nw_goodz.append(sum(dat_gz[sel_fn]['WEIGHT_ZFAIL']))
        n_goodz.append(len(dat_gz[sel_fn]))
    nw_goodz = np.array(nw_goodz)
    n_goodz = np.array(n_goodz)
    n_fail = n_obs-n_goodz
    err = np.sqrt(nw_goodz*n_fail/n_obs)/n_obs
    fig = plt.figure()
    plt.errorbar(fib_obs[sel_4p],nw_goodz[sel_4p]/n_obs[sel_4p],err[sel_4p],fmt='.k')
    plt.plot(fib_obs,np.ones(len(fib_obs))*mean_gz,'r--')
    plt.xlabel('FIBER -500*FIBER//500')
    plt.ylabel('fraction with good z in clus cat '+zw)
    plt.grid()
    figs.append(fig)
    #plt.savefig(outdir+tp+'_stacked_zfailweighted'+zw+'_relsuccess_fiber_zoom.png')
    #plt.clf()

    with PdfPages(outdir+tp+'_zfailweighted_relsuccess_fiber.pdf') as pdf:
        for fig in figs:
            pdf.savefig(fig)
            plt.close()
    print('done with '+tp)



