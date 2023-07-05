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
parser.add_argument("--basedir", help="base directory for catalogs",default='/global/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--mapmd", help="set of maps to use",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--ps",help="point size for density map",default=1,type=float)
parser.add_argument("--test",help="if yes, just use one map from the list",default='n')
parser.add_argument("--dpi",help="resolution in saved density map in dots per inch",default=90,type=int)
args = parser.parse_args()

nside,nest = 256,True


indir = args.basedir+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
outdir = indir+'plots/imaging/'

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

all_maps = ['CALIB_G',
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


 
all_dmaps = [('EBV','EBV_MPF_Mean_FW15'),('EBV','EBV_SGF14')]

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

mf = fitsio.read(indir+'hpmaps/'+tpr+zdw+'_mapprops_healpix_nested_nside256.fits')
sag = np.load('/global/cfs/cdirs/desi/survey/catalogs/extra_regressis_maps/sagittarius_stream_256.npy')

if args.mapmd == 'all':
    maps = all_maps
    dmaps = all_dmaps
    dosag = 'y'
    dosky_g = 'y'
    do_ebvnew_diff = 'y'
    

if args.test == 'y':
    maps = [maps[0]] 
    #print(maps)


nbin = 10

def get_pix(ra, dec):
    return hp.ang2pix(nside, np.radians(-dec+90), np.radians(ra), nest=nest)
    
def plot_reldens(parv,dt_reg,rt_reg,reg,titl='',cl='k',xlab='',yl = (0.8,1.1)):
    dpix = get_pix(dt_reg['RA'],dt_reg['DEC'])
    rpix = get_pix(rt_reg['RA'],rt_reg['DEC'])

    pixlg = np.zeros(nside*nside*12)
    pixlgw = np.zeros(nside*nside*12)
    dcomp = 1/dt_reg['FRACZ_TILELOCID']
    if 'FRAC_TLOBS_TILES' in list(dt_reg.dtype.names):
        print('using FRAC_TLOBS_TILES')
        dcomp *= 1/dt_reg['FRAC_TLOBS_TILES']
    for ii in range(0,len(dpix)):
        pixlg[dpix[ii]] += dt_reg[ii]['WEIGHT_FKP']*dcomp[ii]
        pixlgw[dpix[ii]] += dt_reg[ii]['WEIGHT_FKP']*dt_reg[ii]['WEIGHT_SYS']*dcomp[ii]
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

    meancomp = np.mean(1/dcomp)#np.mean(dt_reg['FRACZ_TILELOCID'])
    ep = np.sqrt(dh/meancomp)/rh*norm #put in mean completeness factor to account for completeness weighting
    
    chi2 = np.sum((svw-1)**2./ep**2.)
    chi2nw = np.sum((sv-1)**2./ep**2.)
    bc = []
    for i in range(0,len(bn)-1):
        bc.append((bn[i]+bn[i+1])/2.)
    labnw = r' no imsys weights, $\chi^2$='+str(round(chi2nw,3))
    labw = r'with imsys weights, $\chi^2$='+str(round(chi2,3))
    #print(lab)    
    plt.errorbar(bc,svw,ep,fmt='o',label=labw,color=cl)
    plt.plot(bc,sv,'-',color=cl,label=labnw)
    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel('Ngal/<Ngal> ')

    plt.title(titl)
    plt.grid()
    plt.ylim(yl[0],yl[1])
    return chi2
        

    

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

    

    yl = (0.8,1.1)    
    if tp == 'LRG':
        z_suc= dtf['ZWARN']==0
        z_suc &= dtf['DELTACHI2']>15
        z_suc &= dtf[zcol]<1.1
        z_suc &= dtf[zcol] > 0.4
        #zr = ' 0.4 < z < 1.1'

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
        #zr = ' 0.8 < z < 2.1 '


    if tp[:3] == 'BGS':    
        z_suc = dtf['ZWARN']==0
        z_suc &= dtf['DELTACHI2']>40
        z_suc &= dtf[zcol]<0.4
        z_suc &= dtf[zcol]>0.1
        #zr = ' 0.1 < z < 0.4 '

    seld &= z_suc

    dtf = dtf[seld]
    tpr = tp
    if tp == 'BGS_BRIGHT-21.5':
        tpr = 'BGS_BRIGHT'
    rf = indir+tpr+zdw+'_0_full.ran.fits'
    rt = fitsio.read(rf)
    
    zbins = [(0.4,0.6),(0.6,0.8),(0.8,1.1)]
    if tp[:3] == 'ELG':
        zbins = [(0.8,1.1),(1.1,1.6)]
    if tp == 'QSO':
        zbins = [(0.8,1.6),(1.6,2.1),(0.8,2.1)]
    if tp[:3] == 'BGS':
        zbins = [(0.1,0.4)]
    for zb in zbins:
        zmin = zb[0]
        zmax = zb[1]
        selz = dtf['Z_not4clus'] > zmin
        selz &= dtf['Z_not4clus'] < zmax
        zr = str(zmin)+'<z<'+str(zmax)       

        for reg,cl in zip(regl,clrs):
            sel_reg_d = dtf['PHOTSYS'] == reg
            sel_reg_r = rt['PHOTSYS'] == reg
            dt_reg = dtf[sel_reg_d&selz]
            rt_reg = rt[sel_reg_r]
            
            #reset for every loop through the maps        
            nside,nest = 256,True
            figs = []
            chi2tot = 0
            nmaptot = 0
            
            if dosag == 'y' and reg == 'S':
        
                parv = sag
                mp = 'sagstream'
                fig = plt.figure()
                chi2 = plot_reldens(parv,dt_reg,rt_reg,cl,reg,titl=args.survey+' '+tp+zr+' '+reg,xlab=mp,yl=yl)
                chi2tot += chi2
                nmaptot += 1
                figs.append(fig)
        #plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
        #plt.clf()
    


            if dosky_g == 'y':
                fig = plt.figure()
                parv = sky_g
                mp = 'g_sky_res'
                
                chi2 = plot_reldens(parv,dt_reg,rt_reg,cl,reg,xlab=mp,titl=args.survey+' '+tp+zr+' '+reg,yl=yl)
                figs.append(fig)
                chi2tot += chi2
                nmaptot += 1

                #plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
                #plt.clf()



        
            for mp in maps:
                fig = plt.figure()
                parv = mf[mp]
                print(mp)
                
                if reg == 'S' or mp[:5] != 'CALIB':
                    chi2 = plot_reldens(parv,dt_reg,rt_reg,cl,reg,yl=yl,xlab=mp,titl=args.survey+' '+tp+zr+' '+reg)
                    chi2tot += chi2
                    nmaptot += 1
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
                mp = map_pair[0]+' - '+map_pair[1]
                chi2 = plot_reldens(parv,dt_reg,rt_reg,cl,reg,yl=yl,xlab=mp,titl=args.survey+' '+tp+zr+' '+reg)
                chi2tot += chi2
                nmaptot += 1

                figs.append(fig)
                #plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
                #plt.clf()
    
            if do_ebvnew_diff == 'y':
                ebvn = fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/test/initial_corrected_ebv_map_nside_64.fits')
                nside = 64
                nest = False
                debv = np.zeros(nside*nside*12)
                for i in range(0,len(ebvn)):
                    pix = ebvn[i]['HPXPIXEL']
                    sfdv = ebvn[i]['EBV_SFD']
                    nv = ebvn[i]['EBV_NEW'] 
                    debv[pix] = nv-sfdv
                parv = debv
                fig = plt.figure()
                plot_reldens(parv,dt_reg,rt_reg,cl,reg,xlab='EBV_RZ - EBV_SFD',titl=args.survey+' '+tp+zr+' '+reg)
                figs.append(fig)
                chi2tot += chi2
                nmaptot += 1
    
       
            tw = ''
            if args.test == 'y':
                tw = '_test'
            with PdfPages(outdir+tp+zr+'_densfullvsall'+tw+' '+reg+'.pdf') as pdf:
                for fig in figs:
                    pdf.savefig(fig)
                    plt.close()
            
            print('results for '+tp+zr+' '+reg)
            print('total chi2 is '+str(chi2)+' for '+str(nmaptot)+ ' maps')
    print('done with '+tp)