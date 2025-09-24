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
parser.add_argument("--version", help="catalog version",default='v1/unblinded')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--ps",help="point size for density map",default=1,type=float)
parser.add_argument("--dpi",help="resolution in saved density map in dots per inch",default=90,type=int)
args = parser.parse_args()


indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
outdir = indir+'plots/angular_power/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

randir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
ranall = fitsio.read(randir+'randoms-allsky-1-0.fits',columns=['RA','DEC'])
th,phi = densvar.radec2thphi(ranall['RA'],ranall['DEC'])
ranpix = hp.ang2pix(256,th,phi,nest=False)
ranpall = np.zeros(12*256*256)
for pix in ranpix:
    ranpall[pix] += 1.

tps = [args.tracers]
if args.tracers == 'all':
    tps = ['LRG','ELG_LOPnotqso','QSO','BGS_BRIGHT']

pix_list = np.arange(12*256*256)
th,phi = hp.pix2ang(256,pix_list)
ra,dec = densvar.thphi2radec(th,phi)


def get_delta(tp,zmin,zmax,reg,racol='RA',decol='DEC',wts='default',thresh=0.1,nest=False,appfrac=True,maskreg=None):#,ranpall=None
    dat = fitsio.read(indir+tp+reg+'_clustering.dat.fits')
    sel = dat['Z'] > zmin
    sel &= dat['Z'] < zmax
    dat = dat[sel]
    ran = fitsio.read(indir+tp+reg+'_0_clustering.ran.fits')
    sel = ran['Z'] > zmin
    sel &= ran['Z'] < zmax
    ran = ran[sel]
    
    th,phi = densvar.radec2thphi(dat[racol],dat[decol])
    datpix = hp.ang2pix(256,th,phi,nest=nest)
    datp = np.zeros(12*256*256)
    for i in range(0,len(datpix)):
        pix = datpix[i]
        if wts == 'default':
            datp[pix] += dat[i]['WEIGHT']
        #else:
        #    datp[pix] += 1.
    th,phi = densvar.radec2thphi(ran[racol],ran[decol])
    ranpix = hp.ang2pix(256,th,phi,nest=nest)
    ranp = np.zeros(12*256*256)
    for i in range(0,len(ranpix)):
        pix = ranpix[i]
        if wts == 'default':
            ranp[pix] += ran[i]['WEIGHT']
    if nest:
        frac = ranp/ranpall_nest
    else:
        frac = ranp/ranpall
    seltest = (frac*0 == 0)
    print('the range for pixel fraction is '+str(min(frac[seltest])),str(max(frac[seltest])))
    sel_thresh = frac > thresh
    sel_thresh &= seltest

    mnr = np.sum(datp[sel_thresh])/np.sum(ranp[sel_thresh])
    print('mean data/random is '+str(mnr))
    delta = (datp/ranp/mnr -1)
    delta *= frac
    delta[~sel_thresh] = hp.UNSEEN
    fsky = np.sum(ranp[sel_thresh])/np.sum(ranpall)
    Ngal = np.sum(datp[sel_thresh])
    return delta,fsky,frac,Ngal

zdw = ''

regl = ['_NGC','_SGC']

tpl = ['LRG','LRG','LRG','ELG_LOPnotqso','ELG_LOPnotqso','QSO','BGS_BRIGHT-21.5']
zrl = [(0.4,0.6),(0.6,0.8),(0.8,1.1),(0.8,1.1),(1.1,1.6),(0.8,2.1),(0.1,0.4)]

for reg in regl:
    for i in range(0,len(tpl)):
        tpi = tpl[i]
        zri = zrl[i]
        delta_i,fsky_i,frac_i,Ngal_i = get_delta(tpi,zri[0],zri[1],reg)
        for j in range(i,len(tpl)):
            tpj = tpl[j]
            zrj = zrl[j]
            fname_out = outdir + tpi+'zr'+str(zri[0])+'-'+str(zri[1])+'_cross_'+tpj+'zr'+str(zrj[0])+'-'+str(zrj[1])+reg
            delta_j,fsky_j,frac_j,Ngal_j = get_delta(tpj,zrj[0],zrj[1],reg)
            cl_ij = hp.anafast(delta_i,delta_j)
            ell = np.arange(len(cl_ij))
            lmax = -300
            fsky_eff = np.sqrt(fsky_i*fsky_j) #I doubt this is actually correct...should somehow be cross-correlation of mask?
            Neff = np.sqrt(Ngal_i*Ngal_j) #I also doubt this is correct...
            plt.loglog(ell[1:lmax],cl_ij[1:lmax]/fsky_eff)#-4.*np.pi*fsky_eff/Neff)
            plt.title(reg.strip('_')+' ' +tpi+' '+str(zri[0])+'<z<'+str(zri[1])+' x '+tpj+' '+str(zrj[0])+'<z<'+str(zrj[1]))
            plt.legend()
            plt.xlabel(r'$\ell$')
            plt.ylabel(r'$C_{\ell}$')
            plt.savefig(fname_out+'_cell.png')
            plt.clf()
            
            

# for tp in tps:
#     print('doing '+tp)
#     dtf = fitsio.read(indir+tp+zdw+'_full.dat.fits')
#     ran = fitsio.read(indir+tp+zdw+'_0_full.ran.fits')
#     fnreg = indir+'/regressis_data/main_'+tp+'_256/RF/main_'+tp+'_imaging_weight_256.npy'
#     rfw = np.load(fnreg,allow_pickle=True)
#     maskreg = rfw.item()['mask_region']
# 
#     #seld = dtf['PHOTSYS'] == reg
#     #dtf = dtf[seld]
#     sel_gz = common.goodz_infull(tp[:3],dtf)
#     sel_obs = dtf['ZWARN'] != 999999
#     dtfoz = dtf[sel_obs&sel_gz]
#     wt = 1./dtfoz['FRACZ_TILELOCID']*dtfoz['WEIGHT_ZFAIL']*dtfoz['WEIGHT_SYS']
#     if 'FRAC_TLOBS_TILES' in list(dtfoz.dtype.names):
#         print('using FRAC_TLOBS_TILES')
#         wt *= 1/dtfoz['FRAC_TLOBS_TILES']
#     
#     sel_nan = wt*0 != 0
#     wt[sel_nan] = 1.
#     if tp[:3] == 'LRG':
#         zmin = 0.4
#         zmax = 1.1
#     if tp[:3] == 'BGS':
#         zmin = 0.1
#         zmax = 0.4
#     if tp[:3] == 'ELG':
#         zmin = 0.8
#         zmax = 1.6
#     if tp[:3] == 'QSO':
#         zmin = 0.8
#         zmax = 2.1
# 
#     sel_zr = dtfoz['Z_not4clus'] > zmin
#     sel_zr &= dtfoz['Z_not4clus'] < zmax
#     delta_raw,fsky,frac = get_delta(dtf,ran,maskreg=maskreg)
#     cl_raw = hp.anafast(delta_raw)
#     ell = np.arange(len(cl_raw))
#     delta_allz,_,_ = get_delta(dtfoz,ran,wts=wt,maskreg=maskreg)
#     cl_allz = hp.anafast(delta_allz)
#     delta_zr,_,_ = get_delta(dtfoz[sel_zr],ran,wts=wt[sel_zr],maskreg=maskreg)
#     cl_zr = hp.anafast(delta_zr)
#     print(len(dtf),np.sum(wt),np.sum(wt[sel_zr]))
#     neff_oz = (np.sum(wt)+len(dtfoz))/2.
#     neff_zr = (np.sum(wt[sel_zr])+len(dtfoz[sel_zr]))/2.
#     lmax = -300
#     plt.loglog(ell[1:lmax],cl_raw[1:lmax]/fsky-4.*np.pi*fsky/len(dtf),label='targets in Y1 area')
#     plt.loglog(ell[1:lmax],cl_allz[1:lmax]/fsky-4.*np.pi*fsky/neff_oz,label='all z')
#     plt.loglog(ell[1:lmax],cl_zr[1:lmax]/fsky-4.*np.pi*fsky/neff_zr,label=str(zmin)+' < z < '+str(zmax))
#     plt.title(tp)
#     plt.legend()
#     plt.xlabel(r'$\ell$')
#     plt.ylabel(r'$C_{\ell}$')
#     plt.savefig(outdir+tp+'_cell.png')
#     plt.clf()
#     print('doing w(theta)')
#     sel = delta_raw != hp.UNSEEN
#     angl,wth_raw = get_wtheta_auto(sindec[sel],cosdec[sel],sinra[sel],cosra[sel],delta_raw[sel],frac[sel])
#     _,wth_allz = get_wtheta_auto(sindec[sel],cosdec[sel],sinra[sel],cosra[sel],delta_allz[sel],frac[sel])
#     _,wth_zr = get_wtheta_auto(sindec[sel],cosdec[sel],sinra[sel],cosra[sel],delta_zr[sel],frac[sel])
# 
#     plt.plot(angl[:-1],1000*angl[:-1]*wth_raw[:-1],label='targets in Y1 area')
#     plt.plot(angl[:-1],1000*angl[:-1]*wth_allz[:-1],label='all z')
#     plt.plot(angl[:-1],1000*angl[:-1]*wth_zr[:-1],label=str(zmin)+' < z < '+str(zmax))
#     plt.grid()
#     plt.title(tp)
#     plt.legend()
#     plt.xlabel(r'$\theta$')
#     plt.ylabel(r'$\theta\times w(\theta)\times 10^3$')
#     plt.savefig(outdir+tp+'_wth.png')
#     plt.clf()
#     
#     
#     regl = list(maskreg.keys())
#     cls = []
#     cls_raw = []
#     wths = []
#     wths_raw = []
#     fskys = []
#     for reg in regl:
#         maskr = maskreg[reg]
#         delta_reg,fsky_reg,frac = get_delta(dtfoz[sel_zr],ran,wts=wt[sel_zr],maskreg=maskr)
#         delta_reg_raw,_,_ = get_delta(dtf,ran,maskreg=maskr)
#         cl_reg = hp.anafast(delta_reg)
#         cls.append(cl_reg)
#         cl_reg_raw = hp.anafast(delta_reg_raw)
#         cls_raw.append(cl_reg_raw)
# 
#         fskys.append(fsky_reg)
#         sel = delta_reg != hp.UNSEEN
#         _,wth_reg = get_wtheta_auto(sindec[sel],cosdec[sel],sinra[sel],cosra[sel],delta_reg[sel],frac[sel])
#         wths.append(wth_reg)
#         _,wth_reg_raw = get_wtheta_auto(sindec[sel],cosdec[sel],sinra[sel],cosra[sel],delta_reg_raw[sel],frac[sel])
#         wths_raw.append(wth_reg_raw)
#        
#     for cl,reg,fsky in zip(cls,regl,fskys):
#         plt.loglog(ell[1:],cl[1:]/fsky,label=reg)
#     plt.title(tp+' '+str(zmin)+' < z < '+str(zmax))
#     plt.legend()
#     plt.xlabel(r'$\ell$')
#     plt.ylabel(r'$C_{\ell}$')
#     plt.savefig(outdir+tp+'_cell_reg.png')
#     plt.clf()
# 
#     for cl,reg,fsky in zip(cls_raw,regl,fskys):
#         plt.loglog(ell[1:],cl[1:]/fsky,label=reg)
#     plt.title(tp+' targets in Y1')
#     plt.legend()
#     plt.xlabel(r'$\ell$')
#     plt.ylabel(r'$C_{\ell}$')
#     plt.savefig(outdir+tp+'_cell_regtar.png')
#     plt.clf()
# 
#     
#     for wth,reg in zip(wths,regl):
#         plt.plot(angl[:-1],1000*angl[:-1]*wth[:-1],label=reg)
#     plt.title(tp+' '+str(zmin)+' < z < '+str(zmax))
#     plt.legend()
#     plt.xlabel(r'$\theta$')
#     plt.ylabel(r'$\theta\times w(\theta)\times 10^3$')
#     plt.savefig(outdir+tp+'_wth_reg.png')
#     plt.clf()
# 
#     for wth,reg in zip(wths_raw,regl):
#         plt.plot(angl[:-1],1000*angl[:-1]*wth[:-1],label=reg)
#     plt.title(tp+' targets in Y1')
#     plt.legend()
#     plt.xlabel(r'$\theta$')
#     plt.ylabel(r'$\theta\times w(\theta)\times 10^3$')
#     plt.savefig(outdir+tp+'_wth_regtar.png')
#     plt.clf()
#     
#         
#         
#     