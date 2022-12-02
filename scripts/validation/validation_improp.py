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
parser.add_argument("--version", help="catalog version",default='EDAbeta')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='DA02')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='guadalupe')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--ps",help="point size for density map",default=1,type=float)
parser.add_argument("--dpi",help="resolution in saved density map in dots per inch",default=90,type=int)
args = parser.parse_args()

pixfn      = '/global/cfs/cdirs/desi/survey/catalogs/pixweight_maps_all/pixweight-1-dark.fits'#'/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/pixweight/sv3/resolve/dark/sv3pixweight-1-dark.fits'
hdr        = fitsio.read_header(pixfn,ext=1)
nside,nest = hdr['HPXNSIDE'],hdr['HPXNEST']
all_maps = fitsio.read(pixfn)

pixfn_ex      = '/global/cfs/cdirs/desi/survey/catalogs/pixweight_maps_all/pixweight_external.fits'
ex_maps = fitsio.read(pixfn_ex)


indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
outdir = indir+'plots/imaging/'

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
clrs = ['r','b']

if args.survey == 'SV3' and args.tracers == 'all':
    tps = ['QSO','LRG','BGS_ANY','BGS_BRIGHT','ELG','ELG_HIP','ELG_HIPnotqso','ELGnotqso']
    zdw = ''
    if args.data != 'LSS':
        tps = ['QSO','LRG','ELG']

maps = ['STARDENS','EBV','GALDEPTH_G', 'GALDEPTH_R','GALDEPTH_Z','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']
maps_ex = ['HALPHA','EBVreconMEANF15','CALIBG', 'CALIBR','CALIBZ']
#EBVdiff = False
EBVdiff = ex_maps['EBVreconMEANF15'] - all_maps['EBV']

nbin = 10

def get_pix(ra, dec):
    return hp.ang2pix(nside, np.radians(-dec+90), np.radians(ra), nest=nest)


for tp in tps:

    if EBVdiff is not None:
        parv = EBVdiff
        map = 'EBVnew-SFD'
        for reg,cl in zip(regl,clrs):
            if reg == '_S' or map[:5] != 'CALIB':
                dtf = fitsio.read(indir+tp+zdw+reg+'_clustering.dat.fits')
                dpix = get_pix(dtf['RA'],dtf['DEC'])
                rf = indir+tp+zdw+reg+'_0_clustering.ran.fits'
                rt = fitsio.read(rf)
                rpix = get_pix(rt['RA'],rt['DEC'])
                pixlg = np.zeros(nside*nside*12)
                pixlgw = np.zeros(nside*nside*12)
                for ii in range(0,len(dpix)):
                    pixlg[dpix[ii]] += dtf[ii]['WEIGHT_COMP']*dtf[ii]['WEIGHT_FKP']
                    pixlgw[dpix[ii]] += dtf[ii]['WEIGHT']*dtf[ii]['WEIGHT_FKP']
                pixlr = np.zeros(nside*nside*12)
                for ii in range(0,len(rpix)):
                    pixlr[rpix[ii]] += 1.
                wp = pixlr > 0
                print(len(parv[wp]))
                rh,bn = np.histogram(parv[wp],bins=nbin,weights=pixlr[wp],range=(np.percentile(parv[wp],1),np.percentile(parv[wp],99)))
                dh,_ = np.histogram(parv[wp],bins=bn,weights=pixlg[wp])
                dhw,_ = np.histogram(parv[wp],bins=bn,weights=pixlgw[wp])
                print(rh)
                print(dh)
                print(dhw)
                norm = sum(rh)/sum(dh)
                sv = dh/rh*norm
                normw = sum(rh)/sum(dhw)
                svw = dhw/rh*normw

                ep = np.sqrt(dh)/rh*norm
                bc = []
                for i in range(0,len(bn)-1):
                    bc.append((bn[i]+bn[i+1])/2.)
                lab = reg.strip('_')+',before weights'
                print(lab)    
                plt.plot(bc,sv,'--',label=lab,color=cl)
                print(bc,svw,ep)
                plt.errorbar(bc,svw,ep,fmt='o',label=reg.strip('_')+', with weights',color=cl)
        plt.legend()
        plt.xlabel(map)
        plt.ylabel('Ngal/<Ngal> ')
    
        plt.title(args.survey+' '+tp)
        plt.grid()
        plt.savefig(outdir+tp+'_densvs'+map+'.png')
        plt.clf()


    for map in maps_ex:
        parv = ex_maps[map]
        print(map)
        for reg,cl in zip(regl,clrs):
            if reg == '_S' or map[:5] != 'CALIB':
                dtf = fitsio.read(indir+tp+zdw+reg+'_clustering.dat.fits')
                dpix = get_pix(dtf['RA'],dtf['DEC'])
                rf = indir+tp+zdw+reg+'_0_clustering.ran.fits'
                rt = fitsio.read(rf)
                rpix = get_pix(rt['RA'],rt['DEC'])
                pixlg = np.zeros(nside*nside*12)
                pixlgw = np.zeros(nside*nside*12)
                for ii in range(0,len(dpix)):
                    pixlg[dpix[ii]] += dtf[ii]['WEIGHT_COMP']*dtf[ii]['WEIGHT_FKP']
                    pixlgw[dpix[ii]] += dtf[ii]['WEIGHT']*dtf[ii]['WEIGHT_FKP']
                pixlr = np.zeros(nside*nside*12)
                for ii in range(0,len(rpix)):
                    pixlr[rpix[ii]] += 1.
                wp = pixlr > 0
                print(len(parv[wp]))
                rh,bn = np.histogram(parv[wp],bins=nbin,weights=pixlr[wp],range=(np.percentile(parv[wp],1),np.percentile(parv[wp],99)))
                dh,_ = np.histogram(parv[wp],bins=bn,weights=pixlg[wp])
                dhw,_ = np.histogram(parv[wp],bins=bn,weights=pixlgw[wp])
                print(rh)
                print(dh)
                print(dhw)
                norm = sum(rh)/sum(dh)
                sv = dh/rh*norm
                normw = sum(rh)/sum(dhw)
                svw = dhw/rh*normw

                ep = np.sqrt(dh)/rh*norm
                bc = []
                for i in range(0,len(bn)-1):
                    bc.append((bn[i]+bn[i+1])/2.)
                lab = reg.strip('_')+',before weights'
                print(lab)    
                plt.plot(bc,sv,'--',label=lab,color=cl)
                print(bc,svw,ep)
                plt.errorbar(bc,svw,ep,fmt='o',label=reg.strip('_')+', with weights',color=cl)
        plt.legend()
        plt.xlabel(map)
        plt.ylabel('Ngal/<Ngal> ')
    
        plt.title(args.survey+' '+tp)
        plt.grid()
        plt.savefig(outdir+tp+'_densvs'+map+'.png')
        plt.clf()

    
    for map in maps:
        parv = all_maps[map]
        for reg,cl in zip(regl,clrs):
            dtf = fitsio.read(indir+tp+zdw+reg+'_clustering.dat.fits')
            dpix = get_pix(dtf['RA'],dtf['DEC'])
            rf = indir+tp+zdw+reg+'_0_clustering.ran.fits'
            rt = fitsio.read(rf)
            rpix = get_pix(rt['RA'],rt['DEC'])
            pixlg = np.zeros(nside*nside*12)
            pixlgw = np.zeros(nside*nside*12)
            for ii in range(0,len(dpix)):
                pixlg[dpix[ii]] += dtf[ii]['WEIGHT_COMP']*dtf[ii]['WEIGHT_FKP']
                pixlgw[dpix[ii]] += dtf[ii]['WEIGHT']*dtf[ii]['WEIGHT_FKP']
            pixlr = np.zeros(nside*nside*12)
            for ii in range(0,len(rpix)):
                pixlr[rpix[ii]] += 1.
            wp = pixlr > 0
            rh,bn = np.histogram(parv[wp],bins=nbin,weights=pixlr[wp],range=(np.percentile(parv[wp],1),np.percentile(parv[wp],99)))
            dh,_ = np.histogram(parv[wp],bins=bn,weights=pixlg[wp])
            dhw,_ = np.histogram(parv[wp],bins=bn,weights=pixlgw[wp])
            
            norm = sum(rh)/sum(dh)
            sv = dh/rh*norm
            normw = sum(rh)/sum(dhw)
            svw = dhw/rh*normw

            ep = np.sqrt(dh)/rh*norm
            bc = []
            for i in range(0,len(bn)-1):
                bc.append((bn[i]+bn[i+1])/2.)
            lab = reg.strip('_')+',before weights'
            print(lab)    
            plt.plot(bc,sv,'--',label=lab,color=cl)
            print(bc,svw,ep)
            plt.errorbar(bc,svw,ep,fmt='o',label=reg.strip('_')+', with weights',color=cl)
        plt.legend()
        plt.xlabel(map)
        plt.ylabel('Ngal/<Ngal> ')
    
        plt.title(args.survey+' '+tp)
        plt.grid()
        plt.savefig(outdir+tp+'_densvs'+map+'.png')
        plt.clf()

