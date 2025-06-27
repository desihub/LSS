import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import sys
import argparse

import fitsio
from astropy.table import join,Table

import LSS.common_tools as common


parser = argparse.ArgumentParser()
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--hpmapcut",help="string to add for hpmapcut or not",default='_HPmapcut')
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


def normed_plot(data,seltot,sel,range=(80,200),col='TSNR2_ELG',wcol='WEIGHT_ZFAIL',cl='k',ps='o',lab=''):
    a,bins = np.histogram(data[seltot][col],range=range)
    b,_ = np.histogram(data[sel][col],bins=bins)
    c,_ = np.histogram(data[sel][col],bins=bins,weights=data[sel][wcol])
    sp = bins[1]-bins[0]
    err = np.sqrt(b*(1-b/a))/a  # binomial error
    normc = np.sum(c)/np.sum(a)
    cchi = np.sum((c/a/normc-1.)**2./(err/normc)**2.)
    plt.errorbar(bins[:-1]+sp/2.,c/a/normc,err/normc,fmt=cl+ps,label=lab+' weighted '+str(round(cchi,3)))
    normb = np.sum(b)/np.sum(a)
    bchi = np.sum((b/a/normb-1.)**2./(err/normb)**2.)
    plt.plot(bins[:-1]+sp/2.,b/a/normb,cl+'--',label=lab+' unweighted '+str(round(bchi,3)))



for tp in tps:
    
    if tp[:3] == 'QSO':
        zmin = 0.8
        zmax = 3.5
        dz = 0.3
    if tp[:3] == 'ELG':
        zmin = 0.6
        zmax = 1.6
        dz = 0.1
    if tp[:3] == 'LRG':
        zmin = 0.4
        zmax = 1.1
        dz = 0.1
    if tp[:3] == 'BGS':
        dz = 0.1
        zmin = 0.1
        zmax = 0.5
    df = fitsio.read(indir+tp+zdw+'_full'+args.hpmapcut+'.dat.fits')


    #dcl = []
    #for reg in regl:
    #    dtc = fitsio.read(indir+tp+zdw+reg+'_clustering.dat.fits')
    #    dcl.append(dtc)
    #dc = np.concatenate(dcl)
    #df = join(dtf,dc,keys=['TARGETID'],join_type='left')
    selgz = common.goodz_infull(tp[:3],df)#df['Z'].mask == False
    selo = df['ZWARN'] != 999999
    selo &= df['ZWARN']*0 == 0
    if tp[:3] == 'QSO':
        selo &= df['PRIORITY'] == 3400 #repeats throw things off
        dchi_cut =40 #was 30
        o2c_cut = 1.2 #was 0.9
        oiii_cut = 5
        selgal = (\
        ( ~selgz & (df['Z_RR'] > 0.01) ) &\ 
            ((df['DELTACHI2'] > dchi_cut) | (np.log10(df['OII_FLUX'] * df['OII_FLUX_IVAR']**0.5) > o2c_cut - 0.2 * np.log10(df['DELTACHI2']))\
            | (np.log10(df['OIII_FLUX'] * df['OIII_FLUX_IVAR'] > oiii_cut))\
            )  #...now 40,1.2,5...earlier 30,0.9,5
        selstar = ( ~selgz & (df['Z_RR'] < 0.01))     
        selo &= ~selgal 
        selo &= ~selstar

    mean_gz = sum(df[selgz]['WEIGHT_ZFAIL'])/len(df[selo])
    print('number with good z, sum of weight_zfail,  number with good obs')
    print(len(df[selgz]),sum(df[selgz]['WEIGHT_ZFAIL']),len(df[selo]))
    if tp[:3] == 'BGS':
        tsnrcol = 'TSNR2_BGS'
        rng=(880,2200)
    else:
        tsnrcol = 'TSNR2_ELG'#+tp[:3]
        rng=(80,200)
    
    zm = zmin
    figs = []
    while zm < zmax:
        selz = df['Z_not4clus'] > zm
        selz &= df['Z_not4clus'] < zm+dz
        seln = df['PHOTSYS'] == 'N'
        fig = plt.figure()
        normed_plot(df,selo&seln,selo&selgz&selz&seln,cl='r',ps='d',lab='N',col=tsnrcol,range=rng)
        normed_plot(df,selo&~seln,selo&selgz&selz&~seln,cl='b',ps='o',lab='S',col=tsnrcol,range=rng)
        plt.legend()
        plt.grid()
        plt.ylabel(tp+' relative z success')
        plt.xlabel(tsnrcol)
        plt.title(str(round(zm,3))+'<z<'+str(round(zm+dz,3)))
        #plt.savefig(outdir+tp+'_'+str(round(zm,3))+'ltzlt'+str(round(zm+dz,3))+'_relsuccess_tnsr.png')
        #plt.clf()
        figs.append(fig)
        zm += dz
    with PdfPages(outdir+tp+'_relsuccess_tnsr_zbins.pdf') as pdf:
        for fig in figs:
            pdf.savefig(fig)
            plt.close()
    print('done with '+tp)

