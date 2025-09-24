#plot map of success/predicted success over focal plane
#uses "full" catalogs as inputs
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
        selgz = df['Z'].mask == False
    else:
        df = dtf
        selgz = common.goodz_infull(tp[:3],df)
    selo = df['ZWARN'] != 999999
    mean_gz = sum(df[selgz]['WEIGHT_ZFAIL'])/len(df[selo])
    print('number with good z, sum of weight_zfail,  number with good obs')
    print(len(df[selgz]),sum(df[selgz]['WEIGHT_ZFAIL']),len(df[selo]))
    focal_r = np.sqrt(df['FIBERASSIGN_X']**2.+df['FIBERASSIGN_Y']**2.)
    cnts_tot,bins = np.histogram(focal_r[selo],bins=20)
    cnts_wt,_= np.histogram(focal_r[selo&selgz],bins=bins,weights=df['WEIGHT_ZFAIL'][selo&selgz])
    cnts_good,_= np.histogram(focal_r[selo&selgz],bins=bins)
    nfail = cnts_tot-cnts_good
    err = np.sqrt(cnts_good*nfail/cnts_tot)/cnts_tot
    bs = bins[1]-bins[0]
    fig = plt.figure()
    plt.errorbar(bins[:-1]+bs/2,cnts_wt/cnts_tot,err,fmt='ko')
    plt.grid()
    plt.xlabel('focal plane radius (mm)')
    plt.ylabel('relative z success (with zfail weight)')
    plt.title(tp+' all petals')
    figs.append(fig)
    #plt.savefig(outdir+tp+'_allpetals_focalr_relsuccess.png')
    #plt.clf()

    #plt.plot(bins[:-1]+bs/2,np.ones(len(fib_obs))*mean_gz,'r--')

    for pt in range(0,10):
        fmin = pt*500
        fmax = fmin+500
        selfib = df['FIBER'] >= fmin
        selfib &= df['FIBER'] < fmax
        cnts_tot,bins = np.histogram(focal_r[selo&selfib],bins=20)
        cnts_wt,_= np.histogram(focal_r[selo&selgz&selfib],bins=bins,weights=df['WEIGHT_ZFAIL'][selo&selgz&selfib])
        cnts_good,_= np.histogram(focal_r[selo&selgz&selfib],bins=bins)
        nfail = cnts_tot-cnts_good
        err = np.sqrt(cnts_good*nfail/cnts_tot)/cnts_tot
        bs = bins[1]-bins[0]
        fig = plt.figure()
        plt.errorbar(bins[:-1]+bs/2,cnts_wt/cnts_tot,err,fmt='ko')
        plt.grid()
        plt.xlabel('focal plane radius (mm)')
        plt.ylabel('relative z success (with zfail weight)')
        plt.title(tp+' petal '+str(pt))
        figs.append(fig)
        #plt.savefig(outdir+tp+'_petal'+str(pt)+'_focalr_relsuccess.png')
        #plt.clf()

    with PdfPages(outdir+tp+'_focalr_relsuccess.pdf') as pdf:
        for fig in figs:
            pdf.savefig(fig)
            plt.close()
    print('done with '+tp)

