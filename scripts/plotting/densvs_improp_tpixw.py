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

parser = argparse.ArgumentParser()
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--ps",help="point size for density map",default=1,type=float)
parser.add_argument("--test",help="if yes, just use one map from the list",default='n')
parser.add_argument("--dpi",help="resolution in saved density map in dots per inch",default=90,type=int)
args = parser.parse_args()

nside,nest = 256,True


indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
outdir = indir+'plots/imaging/'

if args.data == 'LSS':
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print('made '+outdir)


zcol = 'Z_not4clus'
nran = 18

tps = [args.tracers]
if args.tracers == 'all':
    tps = ['ELG_LOPnotqso','BGS_BRIGHT','QSO','LRG']

zdw = ''#'zdone'

regl = ['N','S']
clrs = ['r','b']

maps = ['STARDENS','CALIB_G',
 'CALIB_R',
 'CALIB_Z',
 'EBV_MPF_Mean_FW15',
 'EBV_SGF14',
 'BETA_ML',
 'HI',
 'KAPPA_PLANCK',
 'EBV',
 'PSFDEPTH_G',
 'PSFDEPTH_R',
 'PSFDEPTH_Z',
 'GALDEPTH_G',
 'GALDEPTH_R',
 'GALDEPTH_Z',
 'PSFDEPTH_W1',
 'PSFDEPTH_W2',
 'PSFSIZE_G',
 'PSFSIZE_R',
 'PSFSIZE_Z']

if args.test == 'y':
    maps = maps[0] 
 
dmaps = [('EBV','EBV_MPF_Mean_FW15'),('EBV','EBV_SGF14')]

sky_g = np.zeros(256*256*12)
f = fitsio.read('/global/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/sky_resid_map_256_north.fits')
pixr = f['HPXPIXEL']
pix_nest = hp.ring2nest(256,pixr)
for i in range(0,len(f)):
    pix = pix_nest[i]#f['HPXPIXEL'][i]
    sky_g[pix] = f['sky_median_g'][i]
f = fitsio.read('/global/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/sky_resid_map_256_south.fits')
pix = f['HPXPIXEL']
pix_nest = hp.ring2nest(256,pix)
for i in range(0,len(f)):
    pix = pix_nest[i]#f['HPXPIXEL'][i]
    sky_g[pix] = f['sky_median_g'][i]


nbin = 10

def get_pix(ra, dec):
    return hp.ang2pix(nside, np.radians(-dec+90), np.radians(ra), nest=nest)
    
def plot_reldens(parv,dt_reg,rt_reg,cl,reg):
    dpix = get_pix(dt_reg['RA'],dt_reg['DEC'])
    rpix = get_pix(rt_reg['RA'],rt_reg['DEC'])

    pixlg = np.zeros(nside*nside*12)
    pixlgw = np.zeros(nside*nside*12)
    for ii in range(0,len(dpix)):
        pixlg[dpix[ii]] += 1./dt_reg[ii]['FRACZ_TILELOCID']
        pixlgw[dpix[ii]] += dt_reg[ii]['WEIGHT_SYS']/dt_reg[ii]['FRACZ_TILELOCID']
    pixlr = np.zeros(nside*nside*12)
    for ii in range(0,len(rpix)):
        pixlr[rpix[ii]] += 1.
    wp = pixlr > 0
    wp &= pixlgw*0 == 0
    wp &= parv != hp.UNSEEN
    print(len(parv[wp]))
    rh,bn = np.histogram(parv[wp],bins=nbin,weights=pixlr[wp],range=(np.percentile(parv[wp],1),np.percentile(parv[wp],99)))
    dh,_ = np.histogram(parv[wp],bins=bn,weights=pixlg[wp])
    dhw,_ = np.histogram(parv[wp],bins=bn,weights=pixlgw[wp])
    print((np.percentile(parv[wp],1),np.percentile(parv[wp],99)))
    print(rh)
    print(dh)
    norm = sum(rh)/sum(dh)
    sv = dh/rh*norm
    normw = sum(rh)/sum(dhw)
    svw = dhw/rh*normw

    ep = np.sqrt(dh)/rh*norm
    bc = []
    for i in range(0,len(bn)-1):
        bc.append((bn[i]+bn[i+1])/2.)
    lab = reg+', full, no imsys weights'
    print(lab)    
    plt.errorbar(bc,sv,ep,fmt='o',label=lab,color=cl)
    plt.plot(bc,svw,'-',color=cl,label='with imsys weights')

    

for tp in tps:
    dtf = fitsio.read(indir+tp+zdw+'_full.dat.fits')
    seld = dtf['ZWARN'] != 999999
    seld &= dtf['ZWARN']*0 == 0

    cols = list(dtf.dtype.names)
    if 'Z' in cols:
        print(tp+' Z column already in full file')
        zcol = 'Z'
    else:
        zcol = 'Z_not4clus'

    figs = []

    yl = (0.8,1.1)
    if tp == 'LRG':
        z_suc= dtf['ZWARN']==0
        z_suc &= dtf['DELTACHI2']>15
        z_suc &= dtf[zcol]<1.1
        z_suc &= dtf[zcol] > 0.4
        zr = ' 0.4 < z < 1.1'

    if tp[:3] == 'ELG':
        z_suc = dtf['o2c'] > 0.9
        z_suc &= dtf[zcol]<1.6
        z_suc &= dtf[zcol]>0.8
        zr = ' 0.8 < z < 1.6'
        yl = (0.7,1.1)

    if tp == 'QSO':
        z_suc = dtf[zcol]*0 == 0
        z_suc &= dtf[zcol] != 999999
        z_suc &= dtf[zcol] != 1.e20
        z_suc &= dtf[zcol]<2.1
        z_suc &= dtf[zcol]>0.8
        zr = ' 0.8 < z < 2.1 '


    if tp[:3] == 'BGS':    
        z_suc = dtf['ZWARN']==0
        z_suc &= dtf['DELTACHI2']>40
        z_suc &= dtf[zcol]<0.4
        z_suc &= dtf[zcol]>0.1
        zr = ' 0.1 < z < 0.4 '

    seld &= z_suc

    dtf = dtf[seld]
    rf = indir+tp+zdw+'_0_full.ran.fits'
    rt = fitsio.read(rf)
    
    sag = np.load('/global/cfs/cdirs/desi/survey/catalogs/extra_regressis_maps/sagittarius_stream_256.npy')
    parv = sag
    map = 'sagstream'
    #for reg,cl in zip(regl,clrs):
    reg = 'S'
    cl = clrs[1]        
    sel_reg_d = dtf['PHOTSYS'] == reg
    sel_reg_r = rt['PHOTSYS'] == reg
    dt_reg = dtf[sel_reg_d]
    rt_reg = rt[sel_reg_r]
    plot_reldens(parv,dt_reg,rt_reg,cl,reg)
    plt.legend()
    plt.xlabel(map)
    plt.ylabel('Ngal/<Ngal> ')

    plt.title(args.survey+' '+tp+zr)
    plt.grid()
    plt.ylim(yl[0],yl[1])
    plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
    plt.clf()
    
    
    if sky_g is not None:
        fig = plt.figure()
        parv = sky_g
        map = 'g_sky_res'
        for reg,cl in zip(regl,clrs):
                
            sel_reg_d = dtf['PHOTSYS'] == reg
            sel_reg_r = rt['PHOTSYS'] == reg
            dt_reg = dtf[sel_reg_d]
            rt_reg = rt[sel_reg_r]
            plot_reldens(parv,dt_reg,rt_reg,cl,reg)

        plt.legend()
        plt.xlabel(map)
        plt.ylabel('Ngal/<Ngal> ')
    
        plt.title(args.survey+' '+tp+zr)
        plt.grid()
        plt.ylim(yl[0],yl[1])
        figs.append(fig)
        #plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
        #plt.clf()



    mf = fitsio.read(indir+'hpmaps/'+tp+zdw+'_mapprops_healpix_nested_nside256.fits')
    for map in maps:
        fig = plt.figure()
        parv = mf[map]
        print(map)
        for reg,cl in zip(regl,clrs):
            if reg == 'S' or map[:5] != 'CALIB':
                sel_reg_d = dtf['PHOTSYS'] == reg
                sel_reg_r = rt['PHOTSYS'] == reg
                dt_reg = dtf[sel_reg_d]
                rt_reg = rt[sel_reg_r]
                plot_reldens(parv,dt_reg,rt_reg,cl,reg)
        plt.legend()
        plt.xlabel(map)
        plt.ylabel('Ngal/<Ngal> ')
    
        plt.title(args.survey+' '+tp+zr)
        plt.grid()
        plt.ylim(yl[0],yl[1])
        figs.append(fig)
        #plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
        #plt.clf()
    
    for map_pair in dmaps:
        fig = plt.figure()
        m1 = mf[map_pair[0]]
        m2 = mf[map_pair[1]]
        sel = (m1 == hp.UNSEEN)
        sel |= (m2 == hp.UNSEEN)
        parv = m1-m2
        parv[sel] = hp.UNSEEN
        map = map_pair[0]+' - '+map_pair[1]
        for reg,cl in zip(regl,clrs):
            sel_reg_d = dtf['PHOTSYS'] == reg
            sel_reg_r = rt['PHOTSYS'] == reg
            dt_reg = dtf[sel_reg_d]
            rt_reg = rt[sel_reg_r]
            plot_reldens(parv,dt_reg,rt_reg,cl,reg)
        plt.legend()
        plt.xlabel(map)
        plt.ylabel('Ngal/<Ngal> ')
    
        plt.title(args.survey+' '+tp+zr)
        plt.grid()
        plt.ylim(yl[0],yl[1])
        figs.append(fig)
        #plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
        #plt.clf()
       
    with PdfPages(outdir+tp+'_densfullvsall.pdf') as pdf:
        for fig in figs:
            pdf.savefig(fig)
            plt.close()
