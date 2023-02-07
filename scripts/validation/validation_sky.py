import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import argparse

import fitsio
from astropy.table import join,Table
import healpy as hp

#from LSS.imaging import densvar
import LSS.common_tools as common


parser = argparse.ArgumentParser()
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='SV3')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='fuji')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--ps",help="point size for density map",default=1,type=float)
parser.add_argument("--nside",help="point size for density map",default=64,type=int)
parser.add_argument("--dpi",help="resolution in saved density map in dots per inch",default=90,type=int)
args = parser.parse_args()


indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
outdir = indir+'plots/sky/'

if args.data == 'LSS':
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print('made '+outdir)


qt = 'COMP_TILE'

nside = args.nside
nest = True
zcol = 'Z_not4clus'
nran = 18

tps = [args.tracers]
if args.tracers == 'all':
    tps = ['QSO','LRG','BGS_BRIGHT','ELG_LOPnotqso']

zdw = ''#'zdone'

#regl = ['_N','_S']
regl = ['S','N']

def gethpmap(dl,weights=None):
    rth,rphi = (-dl['DEC']+90.)*np.pi/180.,dl['RA']*np.pi/180. 
    rpix = hp.ang2pix(nside,rth,rphi,nest=nest)
    wts = np.ones(len(rth))
    if weights is not None:
        wts = dl[weights]
    pixlr = np.zeros(12*nside*nside)
    for pix,wt in zip(rpix,wts):
        
        pixlr[pix] += wt
    return pixlr


if args.survey == 'SV3' and args.tracers == 'all':
    tps = ['QSO','LRG','BGS_ANY','BGS_BRIGHT','ELG','ELG_HIP','ELG_HIPnotqso','ELGnotqso']
    zdw = ''
    if args.data != 'LSS':
        tps = ['QSO','LRG','ELG']
for tp in tps:
    
    dtfh = fitsio.read_header(indir+tp+zdw+'_full_noveto.dat.fits',ext=1)
    for nr in range(0,nran):
        rffh = fitsio.read_header(indir+tp+zdw+'_'+str(nr)+'_full_noveto.ran.fits',ext=1)   
        print(tp+' full no veto number density '+str(dtfh['NAXIS2']/rffh['NAXIS2']*2500)+' per deg2, using random '+str(nr))

    dtfh = fitsio.read_header(indir+tp+zdw+'_full.dat.fits',ext=1)
    for nr in range(0,nran):
        rffh = fitsio.read_header(indir+tp+zdw+'_'+str(nr)+'_full.ran.fits',ext=1)   
        print(tp+' full (with veto) number density '+str(dtfh['NAXIS2']/rffh['NAXIS2']*2500)+' per deg2, using random '+str(nr))

#     dtf = fitsio.read(indir+tp+zdw+'_clustering.dat.fits')
#     nc = np.sum(dtf['WEIGHT'])
#     for nr in range(0,nran):
#         rffh = fitsio.read_header(indir+tp+zdw+'_'+str(nr)+'_clustering.ran.fits',ext=1)   
#         print(tp+' weighted clustering number density '+str(nc/rffh['NAXIS2']*2500)+' per deg2, using random '+str(nr))


    if args.data == 'LSS':
        for reg in regl:
            #dtf = fitsio.read(indir+tp+zdw+reg+'_clustering.dat.fits')
            #rf = indir+tp+zdw+reg+'_0_clustering.ran.fits'
            dtf = Table(fitsio.read(indir+tp+zdw+'_full.dat.fits'))
            seld = dtf['PHOTSYS'] == reg
            dtf = dtf[seld]
            sel_gz = common.goodz_infull(tp[:3],dtf)
            sel_obs = dtf['ZWARN'] != 999999
            dtf = dtf[sel_obs&sel_gz]
            #dtf['WEIGHT'] = 1./dtf['FRACZ_TILELOCID']*dtf['WEIGHT_ZFAIL']*dtf['WEIGHT_SYS']
            dtf['WEIGHT'] = 1./dtf['COMP_TILE']*dtf['WEIGHT_ZFAIL']*dtf['WEIGHT_SYS']
            sel_nan = dtf['WEIGHT']*0 != 0
            if len(dtf[sel_nan]) != 0:
                print(str(len(dtf[sel_nan]))+ ' nan weights')
                #print(np.unique(dtf[sel_nan]['FRACZ_TILELOCID']))
                #print(np.unique(dtf[sel_nan]['WEIGHT_ZFAIL']))
                #print(np.unique(dtf[sel_nan]['WEIGHT_SYS']))
            rf = indir+tp+zdw+'_0_full.ran.fits'
            rt = fitsio.read(rf)
            selr = rt['PHOTSYS'] == reg
            rt = rt[selr]
            rad = dtf['RA']
            wr = rad > 300
            rad[wr] -=360
    
            rar = rt['RA']
            wr = rar > 300
            rar[wr] -=360
   
            plt.plot(rad,np.sin(dtf['DEC']*np.pi/180.),'k,',label='data')
            plt.plot(rar,np.sin(rt['DEC']*np.pi/180.),'r,',label='randoms')
            plt.legend(labelcolor='linecolor')
            plt.xlabel('RA')
            plt.ylabel('sin(DEC)')
            plt.title(tp+' randoms plotted over data')
            plt.savefig(outdir+tp+reg+'_ranodat.png')
            plt.clf()

            plt.plot(rar,np.sin(rt['DEC']*np.pi/180.),'r,',label='randoms')
            plt.plot(rad,np.sin(dtf['DEC']*np.pi/180.),'k,',label='data')    
            plt.legend(labelcolor='linecolor')
            plt.xlabel('RA')
            plt.ylabel('sin(DEC)')
            plt.title(tp+' data plotted over randoms')
            plt.savefig(outdir+tp+reg+'_datoran.png')
            plt.clf()
    
    
            titlb = ' weighted over-density for '
            if tp == 'QSO':
                #good redshifts are currently just the ones that should have been defined in the QSO file when merged in full
                wg = dtf[zcol] > 0.8
                wg &= dtf[zcol] < 2.1
                titl = tp +titlb+ '0.8<z<2.1'
    
            if tp[:3] == 'ELG':
                wg = dtf[zcol] > 0.8
                wg &= dtf[zcol] < 1.6
                titl = tp +titlb+ '0.8<z<1.6'

            if tp == 'LRG':
                wg = dtf[zcol] > 0.4
                wg &= dtf[zcol] < 1.1
                titl = tp +titlb+ '0.4<z<1.1'


            if tp[:3] == 'BGS':

                wg = dtf[zcol] > 0.1
                wg &= dtf[zcol] < .4
                titl = tp +titlb+ '0.1<z<0.4'

            dtf = dtf[wg]
            #print(reg,len(dtf))
            rpix = gethpmap(rt)
            dpix = gethpmap(dtf,weights='WEIGHT')
            wp = (rpix > 0) 
            od = dpix[wp]/rpix[wp]
            od = od/np.mean(od)

            pixls = np.arange(12*nside*nside,dtype=int)
            th,phi = hp.pix2ang(nside,pixls[wp],nest=nest)
            ra,dec = 180./np.pi*phi,-(180./np.pi*th-90)#densvar.thphi2radec(th,phi)
            print(np.min(ra),np.max(ra))
    
            if args.survey != 'DA02':
                wr = ra > 300
                ra[wr] -=360
            vx = 1.25
            vm = 0.75
            print(np.min(ra),np.max(ra),np.min(od),np.max(od))
            nside_fac = (256/nside)**2.
            size_fac = 2
            sin_dec = np.sin(dec*np.pi/180)
            yr = (np.max(sin_dec)-np.min(sin_dec))*1.05
            xr = (np.max(ra)-np.min(ra))*1.1/90
            xfac = 2.*size_fac
            yfac = 2.3*size_fac
            fig = plt.figure(figsize=(xr*xfac, yr*yfac))
            ax = fig.add_subplot(111)
            mp = plt.scatter(ra,sin_dec,c=od,edgecolor='none',vmax=vx,vmin=vm,s=args.ps*nside_fac*size_fac,marker='o')
            ax.set_aspect(90)
            plt.colorbar(mp, pad=0.01)
            
            plt.xlabel('RA')
            plt.ylabel('sin(DEC)')
            plt.title(titl)
            plt.grid()


            plt.savefig(outdir+tp+reg+'_weighteddens'+str(nside)+'.png',dpi=args.dpi)
            plt.clf()
            del dtf
            del rt
            print(tp+' done')