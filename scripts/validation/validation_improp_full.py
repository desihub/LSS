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
from LSS import common_tools as common

from regressis import footprint
foot = footprint.DR9Footprint(256, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
north, south, des = foot.get_imaging_surveys()


parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="base directory for catalogs",default='/dvs_ro/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--use_map_veto",help="string to add on the end of full file reflecting if hp maps were used to cut",default='_HPmapcut')
parser.add_argument("--weight_col", help="column name for weight",default='WEIGHT_SYS')
parser.add_argument("--mapmd", help="set of maps to use",default='validate')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--ps",help="point size for density map",default=1,type=float)
parser.add_argument("--test",help="if yes, just use one map from the list",default='n')
parser.add_argument("--dpi",help="resolution in saved density map in dots per inch",default=90,type=int)
args = parser.parse_args()

nside,nest = 256,True


indir = args.basedir+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
outdir = indir+'plots/imaging/'
outdir = outdir.replace('dvs_ro','global')
outdir_txt = outdir+'ngalvsysfiles/'
print('writing to '+outdir)

if not os.path.exists(outdir_txt):
    os.makedirs(outdir_txt)


zcol = 'Z_not4clus'
nran = 18
GCnorm =False

tps = [args.tracers]
#fkpfac_dict = {'ELG_LOPnotqso':.25,'BGS_BRIGHT':0.1,'QSO':1.,'LRG':0.25}
if args.tracers == 'all':
    tps = ['LRG','ELG_LOPnotqso','QSO','BGS_BRIGHT-21.5']
    

zdw = ''#'zdone'

regl = ['S','N']
clrs = ['r','b']

all_maps = ['CALIB_G',
 'CALIB_R',
 'CALIB_Z',
 'STARDENS',
 'HALPHA',
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

lrg_mask_frac = np.zeros(256*256*12)
ranmap = np.zeros(256*256*12)
ranmap_lmask = np.zeros(256*256*12)
randir = '/dvs_ro/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
ran = fitsio.read(randir+'randoms-1-0.fits',columns=['RA','DEC'])
ran_lrgmask = fitsio.read('/dvs_ro/cfs/cdirs/desi/survey/catalogs/main/LSS/randoms-1-0lrgimask.fits')
th,phi = common.radec2thphi(ran['RA'],ran['DEC'])
ranpix = hp.ang2pix(256,th,phi,nest=True)
for pix,mvalue in zip(ranpix,ran_lrgmask['lrg_mask']):
    ranmap[pix] += 1
    if mvalue > 1:
        ranmap_lmask[pix] += 1
sel = ranmap > 0
lrg_mask_frac[sel] = ranmap_lmask[sel]/ranmap[sel]

sky_g = {'N':np.zeros(256*256*12),'S':np.zeros(256*256*12)}
sky_r = {'N':np.zeros(256*256*12),'S':np.zeros(256*256*12)}
sky_z = {'N':np.zeros(256*256*12),'S':np.zeros(256*256*12)}
f = fitsio.read('/dvs_ro/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/sky_resid_map_256_north.fits')
pixr = f['HPXPIXEL']
pix_nest = hp.ring2nest(256,pixr)
for i in range(0,len(f)):
    pix = pix_nest[i]#f['HPXPIXEL'][i]
    sky_g['N'][pix] = f['sky_median_g'][i]
    sky_r['N'][pix] = f['sky_median_r'][i]
    sky_z['N'][pix] = f['sky_median_z'][i]
f = fitsio.read('/dvs_ro/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/sky_resid_map_256_south.fits')
pix = f['HPXPIXEL']
pix_nest = hp.ring2nest(256,pix)
for i in range(0,len(f)):
    pix = pix_nest[i]#f['HPXPIXEL'][i]
    sky_g['S'][pix] = f['sky_median_g'][i]
    sky_r['S'][pix] = f['sky_median_r'][i]
    sky_z['S'][pix] = f['sky_median_z'][i]


sag = np.load('/dvs_ro/cfs/cdirs/desi/survey/catalogs/extra_regressis_maps/sagittarius_stream_256.npy')

if args.mapmd == 'all':
    maps = all_maps
    dmaps = all_dmaps
    dosag = 'y'
    dosky_g = 'y'
    dosky_r = 'y'
    dosky_z = 'y'
    do_ebvnew_diff = 'y'
    do_lrgmask = 'y'
    do_ebvnocib_diff = 'y'
    

if args.test == 'y':
    maps = [maps[0]] 
    #print(maps)


nbin = 10

def get_pix(ra, dec):
    return hp.ang2pix(nside, np.radians(-dec+90), np.radians(ra), nest=nest)
    
def plot_reldens(parv,pixlg,pixlgw,pixlr,titl='',cl='k',xlab='',yl = (0.8,1.1),desnorm=False):
#     from regressis import footprint
#     foot = footprint.DR9Footprint(256, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
#     north, south, des = foot.get_imaging_surveys()
#     dcomp = 1/dt_reg['FRACZ_TILELOCID']
#     dpix = get_pix(dt_reg['RA'],dt_reg['DEC'])
#     rpix = get_pix(rt_reg['RA'],rt_reg['DEC'])
#     seldesr = des[rpix]
#     seldesd = des[dpix]
#     norm_des = np.ones(len(dpix))
#     norm_desw = np.ones(len(dpix))
#     if sum(rpix[seldesr]) > 0 and desnorm:
#         des_ratio = np.sum(dt_reg['WEIGHT_FKP'][seldesd]*dcomp[seldesd])/len(rt_reg[seldesr])
#         notdes_ratio = np.sum(dt_reg['WEIGHT_FKP'][~seldesd]*dcomp[~seldesd])/len(rt_reg[~seldesr])
#         norm_desv = des_ratio/notdes_ratio
#         norm_des[~seldesd] = norm_desv
#         print(norm_desv)
#         des_ratiow = np.sum(dt_reg['WEIGHT_FKP'][seldesd]*dt_reg[args.weight_col][seldesd]*dcomp[seldesd])/len(rt_reg[seldesr])
#         notdes_ratiow = np.sum(dt_reg['WEIGHT_FKP'][~seldesd]*dt_reg[args.weight_col][~seldesd]*dcomp[~seldesd])/len(rt_reg[~seldesr])
#         norm_desvw = des_ratiow/notdes_ratiow
#         norm_desw[~seldesd] = norm_desvw
#         print(norm_desvw)
# 
#     pixlg = np.zeros(nside*nside*12)
#     pixlgw = np.zeros(nside*nside*12)
#     
#     #if 'FRAC_TLOBS_TILES' in list(dt_reg.dtype.names):
#     #    #print('using FRAC_TLOBS_TILES')
#     #    dcomp *= 1/dt_reg['FRAC_TLOBS_TILES']
#     for ii in range(0,len(dpix)):
#         pixlg[dpix[ii]] += dt_reg[ii]['WEIGHT_FKP']*dcomp[ii]*norm_des[ii]
#         pixlgw[dpix[ii]] += dt_reg[ii]['WEIGHT_FKP']*dt_reg[ii][args.weight_col]*dcomp[ii]*norm_desw[ii]
#     pixlr = np.zeros(nside*nside*12)
#     for ii in range(0,len(rpix)):
#         pixlr[rpix[ii]] += rt_reg[ii]['WEIGHT_FKP']*rt_reg[ii]['FRAC_TLOBS_TILES']
    wp = pixlr > 0
    wp &= pixlgw*0 == 0
    wp &= parv != hp.UNSEEN
    #print(len(parv[wp]))

    rh,bn = np.histogram(parv[wp],bins=nbin,weights=pixlr[wp],range=(np.percentile(parv[wp],.1),np.percentile(parv[wp],99.9)))
    dh,_ = np.histogram(parv[wp],bins=bn,weights=pixlg[wp])
    dhw,_ = np.histogram(parv[wp],bins=bn,weights=pixlgw[wp])
    #print((np.percentile(parv[wp],1),np.percentile(parv[wp],99)))
    #print(rh)
    #print(dh)
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

    plt.title(titl+' '+args.weight_col)
    plt.grid()
    plt.ylim(yl[0],yl[1])
    print(xlab,'weighted: '+str(chi2),'unweighted: '+str(chi2nw))
    fname = outdir_txt + 'ngalvs_'+xlab+titl.replace(' ','')+'.txt'
    fname = fname.replace(' - ','')
    fname = fname.replace('<','')
    fo = open(fname,'w')
    for i in range(0,len(bc)):
        fo.write(str(bc[i])+' '+str(svw[i])+' '+str(sv[i])+' '+str(ep[i])+'\n')
    fo.close()
    print('wrote to '+fname)
    return chi2,chi2nw
        

    

for tp in tps:
    if args.weight_col == 'WEIGHT_RF':
        GCnorm = True
    
    depthmd = 'GAL'
    if tp == 'QSO':
        depthmd = 'PSF'
    if args.mapmd == 'validate':
        maps = ['STARDENS','EBV_CHIANG_SFDcorr','HI',depthmd+'DEPTH_G',depthmd+'DEPTH_R',depthmd+'DEPTH_Z','PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']
        dmaps = []
        dosag = 'n'
        dosky_g = 'n'
        dosky_r = 'n'
        dosky_z = 'n'
        do_ebvnew_diff = 'y'
        do_lrgmask = 'n'
        do_ebvnocib_diff = 'n'
        print('doing validation for '+tp)
        
    if tp[:3] == 'ELG' or tp[:3] == 'BGS':
        if 'PSFDEPTH_W1' in maps:
            maps.remove('PSFDEPTH_W1')
    if tp[:3] == 'ELG' or tp[:3] == 'BGS' or tp[:3] == 'LRG':
        if 'PSFDEPTH_W2' in maps:
            maps.remove('PSFDEPTH_W2')

    fcd = indir+tp+zdw+'_full'+args.use_map_veto+'.dat.fits'
    tpr = tp
    if tp == 'BGS_BRIGHT-21.5':
        tpr = 'BGS_BRIGHT'

    rf = indir+tpr+zdw+'_0_full'+args.use_map_veto+'.ran.fits'
    dtf = fitsio.read(fcd)
    seld = dtf['ZWARN'] != 999999
    seld &= dtf['ZWARN']*0 == 0

    cols = list(dtf.dtype.names)
    if 'Z' in cols:
        print(tp+' Z column already in full file')
        zcol = 'Z'
    else:
        zcol = 'Z_not4clus'

    
    zmax = 1.6
    zmin = 0.01
    bs = 0.01

    yl = (0.8,1.1)    
    if tp == 'LRG':
        z_suc= dtf['ZWARN']==0
        z_suc &= dtf['DELTACHI2']>15
        z_suc &= dtf[zcol]<1.1
        z_suc &= dtf[zcol] > 0.4
        P0 = 10000
        nbar = 0.0004
        #zr = ' 0.4 < z < 1.1'

    if tp[:3] == 'ELG':
        z_suc = dtf['o2c'] > 0.9
        z_suc &= dtf[zcol]<1.6
        z_suc &= dtf[zcol]>0.8
        zr = ' 0.8 < z < 1.6'
        yl = (0.7,1.1)
        P0 = 4000
        nbar = 0.0005

    if tp == 'QSO':
        z_suc = dtf[zcol]*0 == 0
        z_suc &= dtf[zcol] != 999999
        z_suc &= dtf[zcol] != 1.e20
        z_suc &= dtf[zcol]<2.1
        z_suc &= dtf[zcol]>0.8
        #zr = ' 0.8 < z < 2.1 '
        zmax = 4
        bs = 0.02
        P0 = 6000
        nbar = 0.00002


    if tp[:3] == 'BGS':    
        z_suc = dtf['ZWARN']==0
        z_suc &= dtf['DELTACHI2']>40
        z_suc &= dtf[zcol]<0.4
        z_suc &= dtf[zcol]>0.1
        #zr = ' 0.1 < z < 0.4 '
        P0 = 7000
        nbar = 0.0005

    #nz = common.mknz_full(fcd,rf,tp,bs=bs,zmin=zmin,zmax=zmax)
    
    seld &= z_suc

    dtf = dtf[seld]
    #zl = dtf[zcol]
    #nl = np.zeros(len(zl))
    #for ii in range(0,len(zl)):
    #    z = zl[ii]
    #    zind = int((z-zmin)/bs)
    #    if z > zmin and z < zmax:
    #        nl[ii] = nz[zind]

    ntl = np.unique(dtf['NTILE'])
    comp_ntl = np.zeros(len(ntl))
    weight_ntl = np.zeros(len(ntl))
    fttl = np.zeros(len(ntl))
    for i in range(0,len(ntl)):
        sel = dtf['NTILE'] == ntl[i]
        mean_ntweight = np.mean(1/dtf[sel]['FRACZ_TILELOCID'])        
        weight_ntl[i] = mean_ntweight
        comp_ntl[i] = 1/mean_ntweight#*mean_fracobs_tiles
        mean_fracobs_tiles = np.mean(dtf[sel]['FRAC_TLOBS_TILES'])
        fttl[i] = mean_fracobs_tiles
    print(comp_ntl,fttl)
    comp_ntl = comp_ntl*fttl
    print('completeness per ntile:')
    print(comp_ntl)
    nx = nbar*comp_ntl[dtf['NTILE']-1]
    fkpl = 1/(1+nx*P0) #this is just the effect of the completeness varying on the fkp weight, no actual z dependence
    dtf = Table(dtf)
    #fd['WEIGHT'] = fd['WEIGHT_COMP']*fd['WEIGHT_SYS']*fd['WEIGHT_ZFAIL']/weight_ntl[fd['NTILE']-1]
    dtf['WEIGHT_FKP'] = 1/weight_ntl[dtf['NTILE']-1]*fkpl
    
    rt = fitsio.read(rf)
    nx = nbar*comp_ntl[rt['NTILE']-1]
    fkpl = 1/(1+nx*P0)
    rt = Table(rt)
    rt['WEIGHT_FKP'] = 1/weight_ntl[rt['NTILE']-1]*fkpl #randoms should now have weight that varies with completeness in same way as data
    
    mf = {'N':fitsio.read(indir+'hpmaps/'+tpr+zdw+'_mapprops_healpix_nested_nside256_N.fits'),\
    'S':fitsio.read(indir+'hpmaps/'+tpr+zdw+'_mapprops_healpix_nested_nside256_S.fits')}
    zbins = [(0.4,0.6),(0.6,0.8),(0.8,1.1)]
    desnorm = False
    if tp[:3] == 'ELG':
        zbins = [(0.8,1.1),(1.1,1.6)]
    if tp == 'QSO':
        zbins = [(0.8,1.6),(1.6,2.1),(0.8,2.1)]
        desnorm=True
        GCnorm = False
    if tp[:3] == 'BGS':
        zbins = [(0.1,0.4)]
    for zb in zbins:
        zmin = zb[0]
        zmax = zb[1]
        selz = dtf['Z_not4clus'] > zmin
        selz &= dtf['Z_not4clus'] < zmax
        zr = str(zmin)+'<z<'+str(zmax)       

        for reg,cl in zip(regl,clrs):
            if args.mapmd == 'validate':
                fo = open(outdir+tp+zr+'_densfullvsall'+'_'+reg+'_'+args.mapmd+args.weight_col+'_chi2.txt','w')
            sel_reg_d = dtf['PHOTSYS'] == reg
            sel_reg_r = rt['PHOTSYS'] == reg
            dt_reg = dtf[sel_reg_d&selz]
            rt_reg = rt[sel_reg_r]
            
            #reset for every loop through the maps        
            nside,nest = 256,True
            figs = []
            chi2tot = 0
            nmaptot = 0

            dcomp = 1/dt_reg['FRACZ_TILELOCID']
            dpix = get_pix(dt_reg['RA'],dt_reg['DEC'])
            rpix = get_pix(rt_reg['RA'],rt_reg['DEC'])

            norm_n = np.ones(len(dpix))
            norm_nw = np.ones(len(dpix))

            
            if reg == 'S' and GCnorm:
                seln = common.splitGC(dt_reg)
                seln_ran = common.splitGC(rt_reg)
                ransum_n = np.sum(rt_reg[seln_ran]['WEIGHT_FKP']*rt_reg[seln_ran]['FRAC_TLOBS_TILES'])
                ransum_s = np.sum(rt_reg[~seln_ran]['WEIGHT_FKP']*rt_reg[~seln_ran]['FRAC_TLOBS_TILES'])
                n_ratio = np.sum(dt_reg['WEIGHT_FKP'][seln]*dcomp[seln])/ransum_n                
                s_ratio = np.sum(dt_reg['WEIGHT_FKP'][~seln]*dcomp[~seln])/ransum_s
                norm_nv = n_ratio/s_ratio
                norm_n[~seln] = norm_nv

                n_ratiow = np.sum(dt_reg['WEIGHT_FKP'][seln]*dt_reg[args.weight_col][seln]*dcomp[seln])/ransum_n                
                s_ratiow = np.sum(dt_reg['WEIGHT_FKP'][~seln]*dt_reg[args.weight_col][~seln]*dcomp[~seln])/ransum_s
                norm_nvw = n_ratiow/s_ratiow
                norm_nw[~seln] = norm_nvw
                print(norm_nvw)

            seldesr = des[rpix]
            seldesd = des[dpix]
            norm_des = np.ones(len(dpix))
            norm_desw = np.ones(len(dpix))
            
            if sum(rpix[seldesr]) > 0 and desnorm:
                ransum_des = np.sum(rt_reg[seldesr]['WEIGHT_FKP']*rt_reg[seldesr]['FRAC_TLOBS_TILES'])
                ransum_notdes = np.sum(rt_reg[~seldesr]['WEIGHT_FKP']*rt_reg[~seldesr]['FRAC_TLOBS_TILES'])

                des_ratio = np.sum(dt_reg['WEIGHT_FKP'][seldesd]*dcomp[seldesd])/ransum_des
                notdes_ratio = np.sum(dt_reg['WEIGHT_FKP'][~seldesd]*dcomp[~seldesd])/ransum_notdes
                norm_desv = des_ratio/notdes_ratio
                norm_des[~seldesd] = norm_desv
                print(norm_desv)
                des_ratiow = np.sum(dt_reg['WEIGHT_FKP'][seldesd]*dt_reg[args.weight_col][seldesd]*dcomp[seldesd])/ransum_des
                
                notdes_ratiow = np.sum(dt_reg['WEIGHT_FKP'][~seldesd]*dt_reg[args.weight_col][~seldesd]*dcomp[~seldesd])/ransum_notdes
                norm_desvw = des_ratiow/notdes_ratiow
                norm_desw[~seldesd] = norm_desvw
                print(norm_desvw)
                #mult = dt_reg['WEIGHT_FKP'][seldesd]*dt_reg[args.weight_col][seldesd]*dcomp[seldesd]
                #print(np.sum(mult)/len)

            pixlg = np.zeros(nside*nside*12)
            pixlgw = np.zeros(nside*nside*12)
    
            #if 'FRAC_TLOBS_TILES' in list(dt_reg.dtype.names):
            #    #print('using FRAC_TLOBS_TILES')
            #    dcomp *= 1/dt_reg['FRAC_TLOBS_TILES']
            for ii in range(0,len(dpix)):
                pixlg[dpix[ii]] += dt_reg[ii]['WEIGHT_FKP']*dcomp[ii]*norm_des[ii]*norm_n[ii]
                pixlgw[dpix[ii]] += dt_reg[ii]['WEIGHT_FKP']*dt_reg[ii][args.weight_col]*dcomp[ii]*norm_desw[ii]*norm_nw[ii]
            pixlr = np.zeros(nside*nside*12)
            for ii in range(0,len(rpix)):
                pixlr[rpix[ii]] += rt_reg[ii]['WEIGHT_FKP']*rt_reg[ii]['FRAC_TLOBS_TILES']

            
            if dosag == 'y' and reg == 'S':
        
                parv = sag
                mp = 'sagstream'
                fig = plt.figure()
                chi2,chi2nw = plot_reldens(parv,pixlg,pixlgw,pixlr,cl=cl,titl=args.survey+' '+tp+zr+' '+reg,xlab=mp,yl=yl,desnorm=desnorm)
                chi2tot += chi2
                nmaptot += 1
                figs.append(fig)
        #plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
        #plt.clf()
    
            if do_lrgmask == 'y':
                fig = plt.figure()
                parv = lrg_mask_frac
                mp = 'fraction of area in LRG mask'
                
                chi2,chi2nw  = plot_reldens(parv,pixlg,pixlgw,pixlr,cl=cl,xlab=mp,titl=args.survey+' '+tp+zr+' '+reg,yl=yl,desnorm=desnorm)
                figs.append(fig)
                chi2tot += chi2
                nmaptot += 1


            if dosky_g == 'y':
                fig = plt.figure()
                parv = sky_g[reg]
                mp = 'g_sky_res'
                
                chi2,chi2nw  = plot_reldens(parv,pixlg,pixlgw,pixlr,cl=cl,xlab=mp,titl=args.survey+' '+tp+zr+' '+reg,yl=yl,desnorm=desnorm)
                figs.append(fig)
                chi2tot += chi2
                nmaptot += 1

            if dosky_r == 'y':
                fig = plt.figure()
                parv = sky_r[reg]
                mp = 'r_sky_res'
                
                chi2,chi2nw  = plot_reldens(parv,pixlg,pixlgw,pixlr,cl=cl,xlab=mp,titl=args.survey+' '+tp+zr+' '+reg,yl=yl,desnorm=desnorm)
                figs.append(fig)
                chi2tot += chi2
                nmaptot += 1

            if dosky_z == 'y':
                fig = plt.figure()
                parv = sky_z[reg]
                mp = 'z_sky_res'
                
                chi2,chi2nw  = plot_reldens(parv,pixlg,pixlgw,pixlr,cl=cl,xlab=mp,titl=args.survey+' '+tp+zr+' '+reg,yl=yl,desnorm=desnorm)
                figs.append(fig)
                chi2tot += chi2
                nmaptot += 1

                #plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
                #plt.clf()

            for map_pair in dmaps:
                fig = plt.figure()
                m1 = mf[reg][map_pair[0]]
                m2 = mf[reg][map_pair[1]]
                sel = (m1 == hp.UNSEEN)
                sel |= (m2 == hp.UNSEEN)
                parv = m1-m2
                parv[sel] = hp.UNSEEN
                mp = map_pair[0]+' - '+map_pair[1]
                chi2,chi2nw  = plot_reldens(parv,pixlg,pixlgw,pixlr,cl=cl,yl=yl,xlab=mp,titl=args.survey+' '+tp+zr+' '+reg,desnorm=desnorm)
                chi2tot += chi2
                nmaptot += 1

                figs.append(fig)
                #plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
                #plt.clf()


        
            for mp in maps:
                fig = plt.figure()
                parv = mf[reg][mp]
                #print(mp)
                #if mp == 'STARDENS':
                #    parv = np.log(parv)
                if reg == 'S' or mp[:5] != 'CALIB':
                    chi2,chi2nw = plot_reldens(parv,pixlg,pixlgw,pixlr,cl=cl,yl=yl,xlab=mp,titl=args.survey+' '+tp+zr+' '+reg,desnorm=desnorm)
                    chi2tot += chi2
                    nmaptot += 1
                    figs.append(fig)
                    if args.mapmd == 'validate':
                        fo.write(str(mp)+' '+str(chi2)+' '+str(chi2nw)+'\n')
                #plt.savefig(outdir+tp+'_densfullvs'+map+'.png')
                #plt.clf()
    
    
            if do_ebvnew_diff == 'y':
                #dirmap = '/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ebv/v0/kp3_maps/'
                dirmap = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/kp3_maps/'
                nside = 256#64
                nest = False
                eclrs = ['gr','rz']
                for ec in eclrs:
                    #ebvn = fitsio.read(dirmap+'v1_desi_ebv_'+ec+'_'+str(nside)+'.fits')
                    ebvn = fitsio.read(dirmap+'v1_desi_ebv_'+str(nside)+'.fits')
                    debv = ebvn['EBV_DESI_'+ec.upper()]-ebvn['EBV_SFD_'+ec.upper()]
                    parv = debv
                    fig = plt.figure()
                    chi2,chi2nw = plot_reldens(parv,hp.reorder(pixlg,n2r=True),hp.reorder(pixlgw,n2r=True),hp.reorder(pixlr,n2r=True),cl=cl,xlab='EBV_DESI_'+ec.upper()+' - EBV_SFD',titl=args.survey+' '+tp+zr+' '+reg,desnorm=desnorm)
                    figs.append(fig)
                    if args.mapmd == 'validate':
                        fo.write('EBV_DESI_'+ec.upper()+'-EBV_SFD'+' '+str(chi2)+' '+str(chi2nw)+'\n')

                    chi2tot += chi2
                    nmaptot += 1

            if do_ebvnocib_diff == 'y':
                #dirmap = '/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ebv/v0/kp3_maps/'
                dirmap = '/global/cfs/cdirs/desicollab/users/rongpu/data/ebv/desi_stars/kp3_maps/'
                nside = 256#64
                nest = False
                eclrs = ['gr','rz']
                ebvnocib = hp.reorder(mf[reg]['EBV_CHIANG_SFDcorr'],n2r=True)
                for ec in eclrs:
                    #ebvn = fitsio.read(dirmap+'v1_desi_ebv_'+ec+'_'+str(nside)+'.fits')
                    ebvn = fitsio.read(dirmap+'v1_desi_ebv_'+str(nside)+'.fits')
                    debv = ebvn['EBV_DESI_'+ec.upper()]-ebvnocib#-ebvn['EBV_SFD_'+ec.upper()]
                    parv = debv
                    fig = plt.figure()
                    chi2,chi2nw = plot_reldens(parv,hp.reorder(pixlg,n2r=True),hp.reorder(pixlgw,n2r=True),hp.reorder(pixlr,n2r=True),cl=cl,xlab='EBV_DESI_'+ec.upper()+' - EBV_SFDnoCIB',titl=args.survey+' '+tp+zr+' '+reg,desnorm=desnorm)
                    figs.append(fig)
                    if args.mapmd == 'validate':
                        fo.write('EBV_DESI_'+ec.upper()+'-EBV_SFDnoCIB'+' '+str(chi2)+' '+str(chi2nw)+'\n')

                    chi2tot += chi2
                    nmaptot += 1
    
       
            tw = ''
            if args.test == 'y':
                tw = '_test'
            with PdfPages(outdir+tp+zr+'_densfullvsall'+tw+'_'+reg+'_'+args.mapmd+args.weight_col+'.pdf') as pdf:
                for fig in figs:
                    pdf.savefig(fig)
                    plt.close()
            
            print('results for '+tp+zr+' '+reg +' using '+args.weight_col+' weights')
            print('total chi2 is '+str(chi2tot)+' for '+str(nmaptot)+ ' maps')
            if args.mapmd == 'validate':
                fo.write('total chi2 is '+str(chi2tot)+' for '+str(nmaptot)+ ' maps\n')
                fo.close()

    print('done with '+tp)