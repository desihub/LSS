import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import os
import sys
import argparse

import fitsio
from astropy.table import join, Table
import healpy as hp

parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="base directory for catalogs", default='/dvs_ro/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--outdir", help="output directory for the plots", default=None)
parser.add_argument("--version", help="catalog version", default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA", default='Y1')
parser.add_argument("--tracers", help="all runs all for given survey", default='all')
parser.add_argument("--verspec", help="version for redshifts", default='iron')
parser.add_argument("--data", help="LSS or mock directory", default='LSS')
parser.add_argument("--use_map_veto", help="string to add on the end of full file reflecting if hp maps were used to cut", default='_HPmapcut')
parser.add_argument("--weight_col", help="column name for the imaging weight; \"off\" if applying no weights", default='WEIGHT_SN')
args = parser.parse_args()

nside, nest = 256, True

indir = os.path.join(args.basedir, args.survey, args.data, args.verspec, 'LSScats', args.version)
if args.outdir is None:
    outdir = os.path.join(indir, 'plots/imaging/').replace('dvs_ro', 'global')
else:
    outdir = args.outdir

if args.data == 'LSS':
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print('made '+outdir)


zcol = 'Z_not4clus'

tps = [args.tracers]
# fkpfac_dict = {'ELG_LOPnotqso':.25, 'BGS_BRIGHT':0.1, 'QSO':1., 'LRG':0.25}
if args.tracers == 'all':
    tps = ['LRG', 'ELG_LOPnotqso', 'QSO', 'BGS_BRIGHT-21.5']


zdw = ''#'zdone'


def get_region(cat):
    if nest:
        des_pixels = np.where(Table.read('/global/cfs/cdirs/desi/users/rongpu/useful/in_des/hp_in_des_{}_nest.fits.gz'.format(nside))['in_des'])[0]
    else:
        des_pixels = np.where(Table.read('/global/cfs/cdirs/desi/users/rongpu/useful/in_des/hp_in_des_{}_ring.fits.gz'.format(nside))['in_des'])[0]

    mask_bm = cat['PHOTSYS']=='N'
    mask_des = np.in1d(cat['HPXPIXEL'], des_pixels)
    mask_decals =(~mask_bm) & (~mask_des)
    return mask_bm, mask_decals, mask_des


def plot_nzsplit(dt, rt, zmin, zmax, zbinsize=0.01):

    nzbin = int(((1.+zbinsize/10.)*(zmax-zmin))/zbinsize)
    for region in ['BASSMZLS', 'DES', 'DECALS']:
        seld = dt[region].copy()
        selr = rt[region].copy()
        area = np.sum(selr)/2500  # area per deg2 if 1 random file being used
        plt.hist(dt[seld]['Z_not4clus'], range=(zmin, zmax), bins=nzbin, weights=dt['dwt'][seld]/area, label=region, histtype='step')


def plot_nzsplit_ratio(dt, rt, zmin, zmax, zbinsize=0.01, ylim=None):
    nzbin = int(((1.+zbinsize/10.)*(zmax-zmin))/zbinsize)
    seld = dt['DECALS'].copy()
    selr = rt['DECALS'].copy()
    area = np.sum(selr)/2500  # area per deg2 if 1 random file being used
    print('area DECaLS, not DES is '+str(area))
    nzdecals, bin_edges = np.histogram(dt[seld]['Z_not4clus'], range=(zmin, zmax), bins=nzbin, weights=dt['dwt'][seld]/area)
    zl = (bin_edges[1:]+bin_edges[:-1])/2

    for region in ['BASSMZLS', 'DES']:
        seld = dt[region].copy()
        selr = rt[region].copy()
        area = np.sum(selr)/2500  # area per deg2 if 1 random file being used
        print('area '+region+' is '+str(area))
        nzp = np.histogram(dt[seld]['Z_not4clus'], range=(zmin, zmax), bins=nzbin, weights=dt['dwt'][seld]/area)[0]
        nzp_no_weight = np.histogram(dt[seld]['Z_not4clus'], range=(zmin, zmax), bins=nzbin)[0]
        # norm = np.sum(nzdecals)/np.sum(nzp)
        # plt.errorbar(zl, nzp/nzdecals*norm, np.sqrt(nzp)/nzdecals*norm, label=region+'/DECALS')
        plt.errorbar(zl, nzp/nzdecals, 1/np.sqrt(nzp_no_weight), label=region+'/DECALS')
    if ylim is not None:
        plt.ylim(ylim[0], ylim[1])


for tp in tps:
    outfn = os.path.join(outdir, tp+'_nzsplit_region.pdf')

    dt = Table(fitsio.read(os.path.join(indir, tp+zdw+'_full'+args.use_map_veto+'.dat.fits')))
    dt['HPXPIXEL'] = hp.ang2pix(nside, dt['RA'], dt['DEC'], nest=nest, lonlat=True)

    seld = dt['ZWARN'] != 999999
    seld &= dt['ZWARN']*0 == 0

    cols = list(dt.dtype.names)
    if 'Z' in cols:
        print(tp+' Z column already in full file')
        zcol = 'Z'
    else:
        zcol = 'Z_not4clus'

    zbinsize = 0.05
    if tp == 'LRG':
        z_suc= dt['ZWARN']==0
        z_suc &= dt['DELTACHI2']>15
        zmax = 1.1
        zmin = 0.4
        # zr = ' 0.4 < z < 1.1'

    if tp[:3] == 'ELG':
        z_suc = dt['o2c'] > 0.9
        zmax = 1.6
        zmin = 0.8

    if tp == 'QSO':
        z_suc = dt[zcol]*0 == 0
        z_suc &= dt[zcol] != 999999
        z_suc &= dt[zcol] != 1.e20
        zmax = 2.1
        zmin = 0.8
        zbinsize = 0.1
        #zr = ' 0.8 < z < 2.1 '

    if tp[:3] == 'BGS':
        z_suc = dt['ZWARN']==0
        z_suc &= dt['DELTACHI2']>40
        zmax = 0.4
        zmin = 0.1

        # zr = ' 0.1 < z < 0.4 '

    seld &= z_suc

    dt = dt[seld]
    tpr = tp
    if tp == 'BGS_BRIGHT-21.5':
        tpr = 'BGS_BRIGHT'
    rf = os.path.join(indir, tpr+zdw+'_0_full'+args.use_map_veto+'.ran.fits')

    rt = Table(fitsio.read(rf))
    rt['HPXPIXEL'] = hp.ang2pix(nside, rt['RA'], rt['DEC'], nest=nest, lonlat=True)

    dt['BASSMZLS'], dt['DECALS'], dt['DES'] = get_region(dt)
    rt['BASSMZLS'], rt['DECALS'], rt['DES'] = get_region(rt)

    dcomp = 1/dt['FRACZ_TILELOCID']
    if 'FRAC_TLOBS_TILES' in list(dt.dtype.names):
        # print('using FRAC_TLOBS_TILES')
        dcomp *= 1/dt['FRAC_TLOBS_TILES']
    dt['dwt'] = dcomp*dt['WEIGHT_ZFAIL']
    if args.weight_col!='off':
        dt['dwt'] *= dt[args.weight_col]

    nside, nest = 256, True
    figs = []

    fig = plt.figure()
    plot_nzsplit(dt, rt, zmin, zmax, zbinsize=zbinsize)
    plt.legend()
    plt.xlabel('redshift')
    plt.ylabel('counts per deg2 ')
    plt.title(args.survey+' '+tp+' ')
    plt.grid()
    figs.append(fig)

    fig = plt.figure()
    plot_nzsplit_ratio(dt, rt, zmin, zmax, zbinsize=zbinsize)
    plt.axhline(1.0, lw=3, alpha=0.5, color='k')
    plt.legend()
    plt.xlabel('redshift')
    plt.ylabel('dN/dz ratio to DECaLS')
    plt.title(args.survey+' '+tp)
    plt.grid()
    figs.append(fig)

    with PdfPages(outfn) as pdf:
        for fig in figs:
            pdf.savefig(fig)
            plt.close()
    print('done with '+tp+'; save to '+outfn)

