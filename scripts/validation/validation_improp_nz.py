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
parser.add_argument("--basedir", help="base directory for catalogs",default='/dvs_ro/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--use_map_veto",help="string to add on the end of full file reflecting if hp maps were used to cut",default='_HPmapcut')
parser.add_argument("--weight_col", help="column name for weight",default='WEIGHT_SN')
parser.add_argument("--nsplit",help="number of percentile bins to split into",default=2,type=int)
parser.add_argument("--mode",help="nz or ratio of nz to unsplit",default='ratio')
parser.add_argument("--ps",help="point size for density map",default=1,type=float)
parser.add_argument("--test",help="if yes, just use one map from the list",default='n')
parser.add_argument("--dpi",help="resolution in saved density map in dots per inch",default=90,type=int)
args = parser.parse_args()

nside,nest = 256,True


indir = args.basedir+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
outdir = indir+'plots/imaging/'
outdir = outdir.replace('dvs_ro','global')
if args.data == 'LSS':
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print('made '+outdir)


zcol = 'Z_not4clus'
nran = 18

tps = [args.tracers]
#fkpfac_dict = {'ELG_LOPnotqso':.25,'BGS_BRIGHT':0.1,'QSO':1.,'LRG':0.25}
if args.tracers == 'all':
    tps = ['LRG','ELG_LOPnotqso','QSO','BGS_BRIGHT-21.5']
    

zdw = ''#'zdone'

regl = ['N','S']
clrs = ['r','b']

maps = ['CALIB_G',
 'CALIB_R',
 'CALIB_Z',
 'STARDENS',
 'EBV',
 'EBV_CHIANG_SFDcorr',
 'EBV_MPF_Mean_FW15',
 'EBV_SGF14',
 'BETA_ML',
 'HI',
 'KAPPA_PLANCK',
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
    maps = [maps[0]] 
    #print(maps)
 
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
    
def plot_nzsplit(parv,dt_reg,rt_reg,zmin,zmax,reg,nsplit=2,zbinsize=0.01):
    dpix = get_pix(dt_reg['RA'],dt_reg['DEC'])
    rpix = get_pix(rt_reg['RA'],rt_reg['DEC'])
    dval = parv[dpix]
    rval = parv[rpix]
    perbs = 100/nsplit
    print(str(np.min(rval))+ ' make sure this minimum value is not null')
    minvl = [np.min(rval)]
    for i in range(1,nsplit):
        minvl.append(np.percentile(rval,i*perbs))
    minvl.append(np.max(rval))
    nzbin = int(((1.+zbinsize/10.)*(zmax-zmin))/zbinsize)
    #print(nzbin)
    dcomp = 1/dt_reg['FRACZ_TILELOCID']
    if 'FRAC_TLOBS_TILES' in list(dt_reg.dtype.names):
        #print('using FRAC_TLOBS_TILES')
        dcomp *= 1/dt_reg['FRAC_TLOBS_TILES']
    dwt = dcomp*dt_reg[args.weight_col]*dt_reg['WEIGHT_ZFAIL']
    for i in range(0,nsplit):
        seld = dval > minvl[i]
        seld &= dval < minvl[i+1]
        selr = rval > minvl[i]
        selr &= rval < minvl[i+1]
        area = len(rval[selr])/2500 #area per deg2 if 1 random file being used
        plt.hist(dt_reg[seld]['Z_not4clus'],range=(zmin,zmax),bins=nzbin,weights=dwt[seld]/area,label='percentile bin '+str(i+1),histtype='step')

def plot_nzsplit_ratio(parv,dt_reg,rt_reg,zmin,zmax,reg,nsplit=2,zbinsize=0.01,ylim=(0.95,1.05)):
    dpix = get_pix(dt_reg['RA'],dt_reg['DEC'])
    rpix = get_pix(rt_reg['RA'],rt_reg['DEC'])
    dval = parv[dpix]
    rval = parv[rpix]
    perbs = 100/nsplit
    print(str(np.min(rval))+ ' make sure this minimum value is not null')
    minvl = [np.min(rval)]
    for i in range(1,nsplit):
        minvl.append(np.percentile(rval,i*perbs))
    minvl.append(np.max(rval))
    nzbin = int(((1.+zbinsize/10.)*(zmax-zmin))/zbinsize)
    zl = np.arange(zmin+zbinsize/2,zmax,zbinsize)
    #print(nzbin)
    dcomp = 1/dt_reg['FRACZ_TILELOCID']
    if 'FRAC_TLOBS_TILES' in list(dt_reg.dtype.names):
        #print('using FRAC_TLOBS_TILES')
        dcomp *= 1/dt_reg['FRAC_TLOBS_TILES']
    dwt = dcomp*dt_reg[args.weight_col]*dt_reg['WEIGHT_ZFAIL']*dt_reg['WEIGHT_FKP']
    nzall = np.histogram(dt_reg['Z_not4clus'],range=(zmin,zmax),bins=nzbin,weights=dwt)
    nzs = []
    for i in range(0,int(nsplit/2)):
        seld = dval > minvl[i]
        seld &= dval < minvl[i+1]
        selr = rval > minvl[i]
        selr &= rval < minvl[i+1]
        area = len(rval[selr])/2500 #area per deg2 if 1 random file being used
        nzp = np.histogram(dt_reg[seld]['Z_not4clus'],range=(zmin,zmax),bins=nzbin,weights=dwt[seld])
        norm = np.sum(nzall[0])/np.sum(nzp[0])
        plt.errorbar(zl,nzp[0]/nzall[0]*norm,np.sqrt(nzp[0])/nzall[0]*norm,label=str(round((i+1)*perbs,2))+ ' percentile ')
    plt.ylim(ylim[0],ylim[1])

  

for tp in tps:
    tw = ''
    if args.test == 'y':
        tw = '_test'
    if args.mode == 'ratio':
        tw += '_ratio'
    outfn = outdir+tp+'_nzsplit'+str(args.nsplit)+tw+'.pdf'  

    dtf = fitsio.read(indir+tp+zdw+'_full'+args.use_map_veto+'.dat.fits')
    seld = dtf['ZWARN'] != 999999
    seld &= dtf['ZWARN']*0 == 0

    cols = list(dtf.dtype.names)
    if 'Z' in cols:
        print(tp+' Z column already in full file')
        zcol = 'Z'
    else:
        zcol = 'Z_not4clus'

    
    zbinsize = 0.05
    if tp == 'LRG':
        z_suc= dtf['ZWARN']==0
        z_suc &= dtf['DELTACHI2']>15
        zmax = 1.1
        zmin =  0.4
        #zr = ' 0.4 < z < 1.1'

    if tp[:3] == 'ELG':
        z_suc = dtf['o2c'] > 0.9
        zmax = 1.6
        zmin = 0.8

    if tp == 'QSO':
        z_suc = dtf[zcol]*0 == 0
        z_suc &= dtf[zcol] != 999999
        z_suc &= dtf[zcol] != 1.e20
        zmax = 2.1
        zmin = 0.8
        zbinsize = 0.1
        #zr = ' 0.8 < z < 2.1 '


    if tp[:3] == 'BGS':    
        z_suc = dtf['ZWARN']==0
        z_suc &= dtf['DELTACHI2']>40
        zmax = 0.4
        zmin = 0.1
        
        #zr = ' 0.1 < z < 0.4 '

    seld &= z_suc

    dtf = dtf[seld]
    tpr = tp
    if tp == 'BGS_BRIGHT-21.5':
        tpr = 'BGS_BRIGHT'
    rf = indir+tpr+zdw+'_0_full'+args.use_map_veto+'.ran.fits'
    rt = fitsio.read(rf)
    
        
    nside,nest = 256,True
    figs = []
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
    fig = plt.figure()
    
    if args.mode == 'nz':
        plot_nzsplit(parv,dt_reg,rt_reg,zmin,zmax,reg,nsplit=args.nsplit,zbinsize=zbinsize)
        plt.legend()
        plt.xlabel('redshift')
        plt.ylabel('counts per deg2 ')
    else:
        plot_nzsplit_ratio(parv,dt_reg,rt_reg,zmin,zmax,reg,nsplit=args.nsplit,zbinsize=zbinsize)
        plt.legend()
        plt.xlabel('redshift')
        plt.ylabel('dN/dz ratio to full sample ')
    
    plt.title(args.survey+' '+tp+' '+map)
    plt.grid()
    figs.append(fig)
    #plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
    #plt.clf()



    if sky_g is not None:
        
        parv = sky_g
        map = 'g_sky_res'
        for reg,cl in zip(regl,clrs):
            fig = plt.figure()
            sel_reg_d = dtf['PHOTSYS'] == reg
            sel_reg_r = rt['PHOTSYS'] == reg
            dt_reg = dtf[sel_reg_d]
            rt_reg = rt[sel_reg_r]
            if args.mode == 'nz':
                plot_nzsplit(parv,dt_reg,rt_reg,zmin,zmax,reg,nsplit=args.nsplit,zbinsize=zbinsize)

                plt.legend()
                plt.xlabel('redshift')
                plt.ylabel('counts per deg2 ')
            else:
                plot_nzsplit_ratio(parv,dt_reg,rt_reg,zmin,zmax,reg,nsplit=args.nsplit,zbinsize=zbinsize)
                plt.legend()
                plt.xlabel('redshift')
                plt.ylabel('dN/dz ratio to full sample ')

            plt.title(args.survey+' '+tp+' '+map+' '+reg)
            plt.grid()
            figs.append(fig)
        #plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
        #plt.clf()



    mfn = fitsio.read(indir+'hpmaps/'+tpr+zdw+'_mapprops_healpix_nested_nside256_N.fits')
    mfs = fitsio.read(indir+'hpmaps/'+tpr+zdw+'_mapprops_healpix_nested_nside256_S.fits')
    mfdic ={'N':mfn,'S':mfs}
    for map in maps:
        
        
        print(map)
        for reg,cl in zip(regl,clrs):
            mf = mfdic[reg]
            parv = mf[map]
            if reg == 'S' or map[:5] != 'CALIB':
                fig = plt.figure()
                sel_reg_d = dtf['PHOTSYS'] == reg
                sel_reg_r = rt['PHOTSYS'] == reg
                dt_reg = dtf[sel_reg_d]
                rt_reg = rt[sel_reg_r]
                if args.mode == 'nz':
                    plot_nzsplit(parv,dt_reg,rt_reg,zmin,zmax,reg,nsplit=args.nsplit,zbinsize=zbinsize)

                    plt.legend()
                    plt.xlabel('redshift')
                    plt.ylabel('counts per deg2 ')
                else:
                    plot_nzsplit_ratio(parv,dt_reg,rt_reg,zmin,zmax,reg,nsplit=args.nsplit,zbinsize=zbinsize)
                    plt.legend()
                    plt.xlabel('redshift')
                    plt.ylabel('dN/dz ratio to full sample ')

                plt.title(args.survey+' '+tp+' '+map+' '+reg)
                plt.grid()
                figs.append(fig)
        #plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
        #plt.clf()

    for map_pair in dmaps:
        
        m1 = mf[map_pair[0]]
        m2 = mf[map_pair[1]]
        sel = (m1 == hp.UNSEEN)
        sel |= (m2 == hp.UNSEEN)
        parv = m1-m2
        parv[sel] = hp.UNSEEN
        map = map_pair[0]+' - '+map_pair[1]
        for reg,cl in zip(regl,clrs):
            fig = plt.figure()
            sel_reg_d = dtf['PHOTSYS'] == reg
            sel_reg_r = rt['PHOTSYS'] == reg
            dt_reg = dtf[sel_reg_d]
            rt_reg = rt[sel_reg_r]
            if args.mode == 'nz':
                plot_nzsplit(parv,dt_reg,rt_reg,zmin,zmax,reg,nsplit=args.nsplit,zbinsize=zbinsize)

                plt.legend()
                plt.xlabel('redshift')
                plt.ylabel('counts per deg2 ')
            else:
                plot_nzsplit_ratio(parv,dt_reg,rt_reg,zmin,zmax,reg,nsplit=args.nsplit,zbinsize=zbinsize)
                plt.legend()
                plt.xlabel('redshift')
                plt.ylabel('dN/dz ratio to full sample ')

            plt.title(args.survey+' '+tp+' '+map+' '+reg)
            plt.grid()
            figs.append(fig)
        #plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
        #plt.clf()


    #        if do_ebvnew_diff == 'y':
    dirmap = '/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ebv/v0/kp3_maps/'
    nside = 256#64
    nest = False
    eclrs = ['gr','rz']
    for ec in eclrs:
        ebvn = fitsio.read(dirmap+'v0_desi_ebv_'+ec+'_'+str(nside)+'.fits')
        debv = ebvn['EBV_DESI_'+ec.upper()]-ebvn['EBV_SFD']
        parv = debv
        for reg,cl in zip(regl,clrs):
            fig = plt.figure()
            sel_reg_d = dtf['PHOTSYS'] == reg
            sel_reg_r = rt['PHOTSYS'] == reg
            dt_reg = dtf[sel_reg_d]
            rt_reg = rt[sel_reg_r]
            if args.mode == 'nz':
                plot_nzsplit(parv,dt_reg,rt_reg,zmin,zmax,reg,nsplit=args.nsplit,zbinsize=zbinsize)

                plt.legend()
                plt.xlabel('redshift')
                plt.ylabel('counts per deg2 ')
            else:
                plot_nzsplit_ratio(parv,dt_reg,rt_reg,zmin,zmax,reg,nsplit=args.nsplit,zbinsize=zbinsize)
                plt.legend()
                plt.xlabel('redshift')
                plt.ylabel('dN/dz ratio to full sample ')

            plt.title(args.survey+' '+tp+' EBV_DESI_'+ec.upper()+' - EBV_SFD '+reg)
            plt.legend()
            plt.grid()
            figs.append(fig)


#     ebvn = fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/test/initial_corrected_ebv_map_nside_64.fits')
#     nside = 64
#     nest = False
#     debv = np.zeros(nside*nside*12)
#     for i in range(0,len(ebvn)):
#         pix = ebvn[i]['HPXPIXEL']
#         sfdv = ebvn[i]['EBV_SFD']
#         nv = ebvn[i]['EBV_NEW'] 
#         debv[pix] = nv-sfdv
#     parv = debv
#     
#     for reg,cl in zip(regl,clrs):
#         fig = plt.figure()
#         sel_reg_d = dtf['PHOTSYS'] == reg
#         sel_reg_r = rt['PHOTSYS'] == reg
#         dt_reg = dtf[sel_reg_d]
#         rt_reg = rt[sel_reg_r]
#         if args.mode == 'nz':
#             plot_nzsplit(parv,dt_reg,rt_reg,zmin,zmax,reg,nsplit=args.nsplit,zbinsize=zbinsize)
# 
#             plt.legend()
#             plt.xlabel('redshift')
#             plt.ylabel('counts per deg2 ')
#         else:
#             plot_nzsplit_ratio(parv,dt_reg,rt_reg,zmin,zmax,reg,nsplit=args.nsplit,zbinsize=zbinsize)
#             plt.legend()
#             plt.xlabel('redshift')
#             plt.ylabel('dN/dz ratio to full sample ')
# 
#         plt.title(args.survey+' '+tp+' EBV_RZ - EBV_SFD '+reg)
#         plt.legend()
#         plt.grid()
#         figs.append(fig)

   
    with PdfPages(outfn) as pdf:
        for fig in figs:
            pdf.savefig(fig)
            plt.close()
    print('done with '+tp)