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
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
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
    tps = ['LRG','ELG_LOPnotqso','QSO','BGS_BRIGHT']

pix_list = np.arange(12*256*256)
th,phi = hp.pix2ang(256,pix_list)
ra,dec = densvar.thphi2radec(th,phi)

sindec = np.sin(dec*np.pi/180.)
cosdec = np.cos(dec*np.pi/180.)

sinra = np.sin(ra*np.pi/180.)
cosra = np.cos(ra*np.pi/180.)

def get_wtheta_auto(sindec,cosdec,sinra,cosra,odens,frac,thmin=0,thmax=10,bs=.1):
    '''
    sines and cosines of ra,dec coordinates, already cut by whatever masking
    overdensity (in same pixels)
    fractional area of same pixels
    '''
    odens /= frac #because it got multiplied by frac for cl
    fo = open('tempodenspczw.dat','w')
    for i in range(len(sindec)):
        fo.write(str(sinra[i])+' '+str(cosra[i])+' '+str(sindec[i])+' '+str(cosdec[i])+' '+str(odens[i])+' '+str(frac[i])+'\n ')
    fo.close()
    os.system('/global/homes/a/ajross/code/LSSanalysis/pix2p_linbin_test temp 1 '+str(bs)+' '+str(thmax))
    res = np.loadtxt('temp2ptPixclb.dat').transpose()
#     nbin = int((thmax-thmin)/bs)
#     odl = np.zeros(nbin)
#     fracl = np.zeros(nbin)
#     binedges = []#np.zeros(nbin+1)
#     th = thmax
#     while th > thmin:
#         be = np.cos(th*np.pi/180.)
#         binedges.append(be)
#         th -= 0.1
#     print(len(binedges),len(odl))
#     bin_angs = np.flip(np.arange(thmin+bs/2.,thmax,bs))
#     for ii in range(0,len(sindec)):
#         for jj in range(ii+1,len(cosdec)):
#             cosang = cosdec[ii]*cosdec[jj]*(cosra[ii]*cosra[jj] + sinra[ii]*sinra[jj]) + sindec[ii]*sindec[jj]
#             be = binedges[0]
#             ba = -1 #start at -1 because of condition below
#             while cosang > be:
#                 ba += 1
#                 be = binedges[ba+1]
#             
#             if ba > -1 and ba < nbin:
#                 odl[ba] += odens[ii]*odens[jj] #note, frac was already applied to odl
#                 fracl[ba] += frac[ii]*frac[jj] 
    return res[0],res[1]#bin_angs,odl/fracl

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
        mnr = np.sum(datp[sel])/np.sum(ranp[sel])
        print(mnr)
        delta = (datp/ranp/mnr -1)
    elif len(maskreg)==len(datp):
        if nest == False:
            maskreg = hp.reorder(maskreg,n2r=True)

        sel &= maskreg
        mnr = np.sum(datp[sel])/np.sum(ranp[sel])
        delta = (datp/ranp/mnr -1)
        
    else:
        regl = list(maskreg.keys())#['South','North','Des']
        delta = np.zeros(len(datp))
        for reg in regl:
            mr = maskreg[reg]
            if nest == False:
                mr = hp.reorder(mr,n2r=True)

            mnr = np.sum(datp[sel&mr])/np.sum(ranp[sel&mr]) 
            print(reg,mnr)
            delta[mr] = (datp[mr]/ranp[mr]/mnr -1)
    #if ranpall is not None:
    #if appfrac:
    if nest:
        frac = ranp/ranpall_nest
    else:
        frac = ranp/ranpall
    delta *= frac
    delta[~sel] = hp.UNSEEN
    fsky = np.sum(ranp[sel])/np.sum(ranpall)
    return delta,fsky,frac

zdw = ''

for tp in tps:
    print('doing '+tp)
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
    if 'FRAC_TLOBS_TILES' in list(dtfoz.dtype.names):
        print('using FRAC_TLOBS_TILES')
        wt *= 1/dtfoz['FRAC_TLOBS_TILES']
    
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
    delta_raw,fsky,frac = get_delta(dtf,ran,maskreg=maskreg)
    cl_raw = hp.anafast(delta_raw)
    ell = np.arange(len(cl_raw))
    delta_allz,_,_ = get_delta(dtfoz,ran,wts=wt,maskreg=maskreg)
    cl_allz = hp.anafast(delta_allz)
    delta_zr,_,_ = get_delta(dtfoz[sel_zr],ran,wts=wt[sel_zr],maskreg=maskreg)
    cl_zr = hp.anafast(delta_zr)
    print(len(dtf),np.sum(wt),np.sum(wt[sel_zr]))
    neff_oz = (np.sum(wt)+len(dtfoz))/2.
    neff_zr = (np.sum(wt[sel_zr])+len(dtfoz[sel_zr]))/2.
    lmax = -300
    plt.loglog(ell[1:lmax],cl_raw[1:lmax]/fsky-4.*np.pi*fsky/len(dtf),label='targets in Y1 area')
    plt.loglog(ell[1:lmax],cl_allz[1:lmax]/fsky-4.*np.pi*fsky/neff_oz,label='all z')
    plt.loglog(ell[1:lmax],cl_zr[1:lmax]/fsky-4.*np.pi*fsky/neff_zr,label=str(zmin)+' < z < '+str(zmax))
    plt.title(tp)
    plt.legend()
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C_{\ell}$')
    plt.savefig(outdir+tp+'_cell.png')
    plt.clf()
    print('doing w(theta)')
    sel = delta_raw != hp.UNSEEN
    angl,wth_raw = get_wtheta_auto(sindec[sel],cosdec[sel],sinra[sel],cosra[sel],delta_raw[sel],frac[sel])
    _,wth_allz = get_wtheta_auto(sindec[sel],cosdec[sel],sinra[sel],cosra[sel],delta_allz[sel],frac[sel])
    _,wth_zr = get_wtheta_auto(sindec[sel],cosdec[sel],sinra[sel],cosra[sel],delta_zr[sel],frac[sel])

    plt.plot(angl[:-1],1000*angl[:-1]*wth_raw[:-1],label='targets in Y1 area')
    plt.plot(angl[:-1],1000*angl[:-1]*wth_allz[:-1],label='all z')
    plt.plot(angl[:-1],1000*angl[:-1]*wth_zr[:-1],label=str(zmin)+' < z < '+str(zmax))
    plt.grid()
    plt.title(tp)
    plt.legend()
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$\theta\times w(\theta)\times 10^3$')
    plt.savefig(outdir+tp+'_wth.png')
    plt.clf()
    
    
    regl = list(maskreg.keys())
    cls = []
    cls_raw = []
    wths = []
    wths_raw = []
    fskys = []
    for reg in regl:
        maskr = maskreg[reg]
        delta_reg,fsky_reg,frac = get_delta(dtfoz[sel_zr],ran,wts=wt[sel_zr],maskreg=maskr)
        delta_reg_raw,_,_ = get_delta(dtf,ran,maskreg=maskr)
        cl_reg = hp.anafast(delta_reg)
        cls.append(cl_reg)
        cl_reg_raw = hp.anafast(delta_reg_raw)
        cls_raw.append(cl_reg_raw)

        fskys.append(fsky_reg)
        sel = delta_reg != hp.UNSEEN
        _,wth_reg = get_wtheta_auto(sindec[sel],cosdec[sel],sinra[sel],cosra[sel],delta_reg[sel],frac[sel])
        wths.append(wth_reg)
        _,wth_reg_raw = get_wtheta_auto(sindec[sel],cosdec[sel],sinra[sel],cosra[sel],delta_reg_raw[sel],frac[sel])
        wths_raw.append(wth_reg_raw)
       
    for cl,reg,fsky in zip(cls,regl,fskys):
        plt.loglog(ell[1:],cl[1:]/fsky,label=reg)
    plt.title(tp+' '+str(zmin)+' < z < '+str(zmax))
    plt.legend()
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C_{\ell}$')
    plt.savefig(outdir+tp+'_cell_reg.png')
    plt.clf()

    for cl,reg,fsky in zip(cls_raw,regl,fskys):
        plt.loglog(ell[1:],cl[1:]/fsky,label=reg)
    plt.title(tp+' targets in Y1')
    plt.legend()
    plt.xlabel(r'$\ell$')
    plt.ylabel(r'$C_{\ell}$')
    plt.savefig(outdir+tp+'_cell_regtar.png')
    plt.clf()

    
    for wth,reg in zip(wths,regl):
        plt.plot(angl[:-1],1000*angl[:-1]*wth[:-1],label=reg)
    plt.title(tp+' '+str(zmin)+' < z < '+str(zmax))
    plt.legend()
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$\theta\times w(\theta)\times 10^3$')
    plt.savefig(outdir+tp+'_wth_reg.png')
    plt.clf()

    for wth,reg in zip(wths_raw,regl):
        plt.plot(angl[:-1],1000*angl[:-1]*wth[:-1],label=reg)
    plt.title(tp+' targets in Y1')
    plt.legend()
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$\theta\times w(\theta)\times 10^3$')
    plt.savefig(outdir+tp+'_wth_regtar.png')
    plt.clf()
    
        
        
    