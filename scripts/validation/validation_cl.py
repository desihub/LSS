import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import argparse

import fitsio
from astropy.table import join,Table
import healpy as hp

from LSS.imaging import densvar
import LSS.common_tools as common


parser = argparse.ArgumentParser()
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='SV3')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='fuji')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--ps",help="point size for density map",default=1,type=float)
parser.add_argument("--dpi",help="resolution in saved density map in dots per inch",default=90,type=int)
args = parser.parse_args()


indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
outdir = indir+'plots/sky/'

randir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
ranall = fitsio.read(randir+'randoms-allsky-1-0.fits',columns=['RA','DEC'])
th,phi = densvar.radec2thphi(ranall['RA'],ranall['DEC'])
ranpix = hp.ang2pix(256,th,phi)
ranpall = np.zeros(12*256*256)
for pix in ranpix:
    ranpall[pix] += 1.

tps = [args.tracers]
if args.tracers == 'all':
    tps = ['QSO','LRG','BGS_BRIGHT','ELG_LOPnotqso']


def get_delta(dat,ran,racol='RA',decol='DEC',wts=None,wtspix=None,thresh=0,nest=False,appfrac=True,maskreg=None):#,ranpall=None
    th,phi = densvar.radec2thphi(dat[racol],dat[decol])
    datpix = hp.ang2pix(256,th,phi,nest=nest)
    datp = np.zeros(12*256*256)
    for i in range(0,len(datpix)):
        pix = datpix[i]
        if wts is not None:
            datp[pix] += wts[i]
        else:
            datp[pix] += 1.
    if wtspix is not None:
        datp *= wtspix
    th,phi = densvar.radec2thphi(ran[racol],ran[decol])
    ranpix = hp.ang2pix(256,th,phi,nest=nest)
    ranp = np.zeros(12*256*256)
    for pix in ranpix:
        ranp[pix] += 1.
    #ranp /= rannorm
    
    sel = ranp > thresh
    if maskreg is None:
        mnr = np.mean(datp[sel]/ranp[sel])
        print(mnr)
        delta = (datp/ranp/mnr -1)
    elif len(maskreg)==len(datp):
        sel &= maskreg
        mnr = np.mean(datp[sel]/ranp[sel])
        delta = (datp/ranp/mnr -1)
        
    else:
        regl = list(maskreg.keys())#['South','North','Des']
        delta = np.zeros(len(datp))
        for reg in regl:
            mr = maskreg[reg]
            mnr = np.mean(datp[sel&mr]/ranp[sel&mr])
            print(reg,mnr)
            delta[mr] = (datp[mr]/ranp[mr]/mnr -1)
    #if ranpall is not None:
    if appfrac:
        if nest:
            frac = ranp/ranpall_nest
        else:
            frac = ranp/ranpall
        delta *= frac
    delta[~sel] = hp.UNSEEN
    fsky = np.sum(ranp[sel])/np.sum(ranpall)
    return delta,fsky

zdw = ''

for tp in tps:
    
    dtf = fitsio.read(indir+tp+zdw+'_full.dat.fits')
    ran = fitsio.read(indir+tp+zdw+'_0_full.ran.fits')
    fnreg = indir+'/regressis_data/main_'+tp+'_256/RF/main_'+tp+'_imaging_weight_256.npy'
    rfw = np.load(fnreg,allow_pickle=True)
    maskreg = rfw.item()['mask_region']

    #seld = dtf['PHOTSYS'] == reg
    #dtf = dtf[seld]
    sel_gz = common.goodz_infull(tp[:3],dtf)
    sel_obs = dtf['ZWARN'] != 999999
    dtfoz = dtf[sel_obs&sel_gz]
    wt = 1./dtfoz['FRACZ_TILELOCID']*dtfoz['WEIGHT_ZFAIL']*dtfoz['WEIGHT_SYS']
    sel_nan = wt*0 != 0
    wt[sel_nan] = 1.
    if tp[:3] == 'LRG':
        zmin = 0.4
        zmax = 1.1
    if tp[:3] == 'BGS':
        zmin = 0.1
        zmax = 0.4
    if tp[:3] == 'ELG':
        zmin = 0.8
        zmax = 1.6
    if tp[:3] == 'QSO':
        zmin = 0.8
        zmax = 2.1

    sel_zr = dtfoz['Z_not4clus'] > zmin
    sel_zr &= dtfoz['Z_not4clus'] < zmax
    delta_raw,fsky = get_delta(dtf,ran,maskreg=maskreg)
    cl_raw = hp.anafast(delta_raw)
    ell = np.arange(len(cl_raw))
    delta_allz,_ = get_delta(dtfoz,ran,wts=wt,maskreg=maskreg)
    cl_allz = hp.anafast(delta_allz)
    delta_zr,_ = get_delta(dtfoz[sel_zr],ran,wts=wt[sel_zr],maskreg=maskreg)
    cl_zr = hp.anafast(delta_zr)
    plt.loglog(ell[1:],cl_raw[1:]/fsky-4.*np.pi/len(dtf),label='targets in Y1 area')
    plt.loglog(ell[1:],cl_allz[1:]/fsky-4.*np.pi/len(dtfoz),label='all z')
    plt.loglog(ell[1:],cl_zr[1:]/fsky--4.*np.pi/len(dtf[sel_zr]),label=str(zmin)+' < z < '+str(zmax))
    plt.title(tp)
    plt.legend()
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C_{\ell}$')
    plt.savefig(outdir+tp+'_cell.png')
    plt.clf()
    regl = list(maskreg.keys())
    for reg in regl:
        maskr = maskreg[reg]
        delta_reg,fsky_reg = get_delta(dtfoz[sel_zr],ran,wts=wt[sel_zr],maskreg=maskr)
        cl_reg = hp.anafast(delta_reg)
        plt.loglog(ell[1:],cl_reg[1:]/fsky_reg,label=reg)
    plt.title(tp+' '+str(zmin)+' < z < '+str(zmax))
    plt.legend()
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C_{\ell}$')
    plt.savefig(outdir+tp+'_cell_reg.png')
    plt.clf()

        
        
    