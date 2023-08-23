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
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--use_map_veto",help="string to add on the end of full file reflecting if hp maps were used to cut",default='_HPmapcut')
parser.add_argument("--compmd",help="extra completeness on data or random",default='ran')
parser.add_argument("--ps",help="point size for density map",default=.1,type=float)
parser.add_argument("--nside",help="point size for density map",default=64,type=int)
parser.add_argument("--dpi",help="resolution in saved density map in dots per inch",default=90,type=int)
args = parser.parse_args()


indir = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
outdir = indir+'plots/sky/'
outdir = outdir.replace('dvs_ro','global')
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
    tps = ['ELG_LOPnotqso','LRG','QSO','BGS_BRIGHT-21.5']

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

def plot_map_sindec(ra,sin_dec,od,vm,vx,titl,outf,size_fac=2):
	yr = (np.max(sin_dec)-np.min(sin_dec))*1.05
	xr = (np.max(ra)-np.min(ra))*1.1/90
	xfac = 2.*size_fac
	yfac = 2.3*size_fac
	fig = plt.figure(figsize=(xr*xfac, yr*yfac))
	ax = fig.add_subplot(111)
	mp = plt.scatter(ra,sin_dec,c=od,edgecolor='none',vmax=vx,vmin=vm,s=args.ps)#,marker=',')#s=args.ps*nside_fac*size_fac,marker='o')
	ax.set_aspect(90)
	plt.colorbar(mp, pad=0.01,shrink=2/2.3)
	
	plt.xlabel('RA')
	plt.ylabel('sin(DEC)')
	plt.title(titl)
	plt.grid()


	plt.savefig(outf)#,dpi=args.dpi)
	plt.clf()



if args.survey == 'SV3' and args.tracers == 'all':
    tps = ['QSO','LRG','BGS_ANY','BGS_BRIGHT','ELG','ELG_HIP','ELG_HIPnotqso','ELGnotqso']
    zdw = ''
    if args.data != 'LSS':
        tps = ['QSO','LRG','ELG']
for tp in tps:
    tpr = tp
    if tp == 'BGS_BRIGHT-21.5':
        tpr = 'BGS_BRIGHT'
    dtfh = fitsio.read_header(indir+tpr+zdw+'_full_noveto.dat.fits',ext=1)
    for nr in range(0,nran):
        rffh = fitsio.read_header(indir+tpr+zdw+'_'+str(nr)+'_full_noveto.ran.fits',ext=1)   
        print(tpr+' full no veto number density '+str(dtfh['NAXIS2']/rffh['NAXIS2']*2500)+' per deg2, using random '+str(nr))

    dtfh = fitsio.read_header(indir+tpr+zdw+'_full.dat.fits',ext=1)
    for nr in range(0,nran):
        rffh = fitsio.read_header(indir+tpr+zdw+'_'+str(nr)+'_full.ran.fits',ext=1)   
        print(tpr+' full (with veto) number density '+str(dtfh['NAXIS2']/rffh['NAXIS2']*2500)+' per deg2, using random '+str(nr))

#     dtf = fitsio.read(indir+tp+zdw+'_clustering.dat.fits')
#     nc = np.sum(dtf['WEIGHT'])
#     for nr in range(0,nran):
#         rffh = fitsio.read_header(indir+tp+zdw+'_'+str(nr)+'_clustering.ran.fits',ext=1)   
#         print(tp+' weighted clustering number density '+str(nc/rffh['NAXIS2']*2500)+' per deg2, using random '+str(nr))


    if args.data == 'LSS':
        ral = []
        sdecl = []
        odl = []
        odl_oc = []
        dt = Table(fitsio.read(indir+tp+zdw+'_full'+args.use_map_veto+'.dat.fits'))
        cols = list(dt.dtype.names)
        sel_gz = common.goodz_infull(tp[:3],dt)
        sel_obs = dt['ZWARN'] != 999999
        dt = dt[sel_obs&sel_gz]
        dt['WEIGHT_COMP'] = 1./dt['FRACZ_TILELOCID']
        if 'FRAC_TLOBS_TILES' in cols and args.compmd == 'dat':
            dt['WEIGHT_COMP'] *= 1/dt['FRAC_TLOBS_TILES']

        dt['WEIGHT'] = dt['WEIGHT_COMP']*dt['WEIGHT_ZFAIL']*dt['WEIGHT_SYS']
        sel_nan = dt['WEIGHT']*0 != 0
        if len(dt[sel_nan]) != 0:
            print(str(len(dt[sel_nan]))+ ' nan weights')
        rada = dt['RA']
        wr = rada > 300
        rada[wr] -=360

 
        rf = indir+tpr+zdw+'_0_full'+args.use_map_veto+'.ran.fits'
        rta = fitsio.read(rf)

        rara = rta['RA']
        wr = rara > 300
        rara[wr] -=360
        
        sindd = np.sin(dt['DEC']*np.pi/180.)
        sindr = np.sin(rta['DEC']*np.pi/180.)

        vm = np.min(dt['COMP_TILE'])
        vx = np.max(dt['COMP_TILE'])
        titl  = tp+' COMP_TILE'
        outf = outdir+tp+'_comptile.png'
        plot_map_sindec(rada,sindd,dt['COMP_TILE'],vm,vx,titl,outf)
        #vm = np.min(dt['WEIGHT_SYS'])
        #vx = np.max(dt['WEIGHT_SYS'])
        vm = 0.75
        vx = 1.25
        titl = tp+' WEIGHT_SYS'
        outf = outdir+tp+'_weightsys.png'
        plot_map_sindec(rada,sindd,dt['WEIGHT_SYS'],vm,vx,titl,outf)

        for reg in regl:
            
            seld = dt['PHOTSYS'] == reg
            dtf = dt[seld]
            rad = rada[seld]
            
            selr = rta['PHOTSYS'] == reg
            rt = rta[selr]
            rar = rara[selr]
    
   
            plt.plot(rad,sindd[seld],'k,',label='data')
            plt.plot(rar,sindr[selr],'r,',label='randoms')
            plt.legend(labelcolor='linecolor')
            plt.xlabel('RA')
            plt.ylabel('sin(DEC)')
            plt.title(tp+' randoms plotted over data')
            plt.savefig(outdir+tp+reg+'_ranodat.png')
            plt.clf()

            plt.plot(rar,sindr[selr],'r,',label='randoms')
            plt.plot(rad,sindd[seld],'k,',label='data')    
            plt.legend(labelcolor='linecolor')
            plt.xlabel('RA')
            plt.ylabel('sin(DEC)')
            plt.title(tp+' data plotted over randoms')
            plt.savefig(outdir+tp+reg+'_datoran.png')
            plt.clf()


    
    
            titlb = ' weighted relative density for '
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
            if args.compmd == 'dat':
                rpix = gethpmap(rt)
            else:
                rpix = gethpmap(rt,weights='FRAC_TLOBS_TILES')
            dpix = gethpmap(dtf,weights='WEIGHT')
            dpix_oc = gethpmap(dtf,weights='WEIGHT_COMP')
            wp = (rpix > 0) 
            od = dpix/rpix
            od = od/np.mean(od[wp])
            od_oc = dpix_oc/rpix
            od_oc = od_oc/np.mean(od_oc[wp])
            rth,rphi = (-dtf['DEC']+90.)*np.pi/180.,dtf['RA']*np.pi/180. 
            rpix = hp.ang2pix(nside,rth,rphi,nest=nest)
            odd = np.zeros(len(rpix))
            odd = od[rpix]
            odl.append(odd)
            odd_oc = np.zeros(len(rpix))
            odd_oc = od_oc[rpix]
            odl_oc.append(odd_oc)
            pixls = np.arange(12*nside*nside,dtype=int)
            th,phi = hp.pix2ang(nside,pixls[wp],nest=nest)
            ra,dec = 180./np.pi*phi,-(180./np.pi*th-90)#densvar.thphi2radec(th,phi)
            print(np.min(ra),np.max(ra))
    
            if args.survey != 'DA02':
                #wr = ra > 300
                #ra[wr] -=360
                rad = dtf['RA']
                wr = rad > 300
                rad[wr] -=360

            ral.append(rad)
            sin_dec = np.sin(dtf['DEC']*np.pi/180)#np.sin(dec*np.pi/180)
            sdecl.append(sin_dec)
            del dtf
            del rt
            
        ra = np.concatenate(ral)
        sin_dec = np.concatenate(sdecl)
        od = np.concatenate(odl)
        od_oc = np.concatenate(odl_oc)
        vx = 1.25
        vm = 0.75
        print(np.min(ra),np.max(ra),np.min(od),np.max(od))
        nside_fac = (256/nside)**2.
        
        outf = outdir+tp+'_weighteddens'+str(nside)+'.png'
        plot_map_sindec(ra,sin_dec,od,vm,vx,titl,outf)
        outf = outdir+tp+'_componly_weighteddens'+str(nside)+'.png'
        plot_map_sindec(ra,sin_dec,od_oc,vm,vx,titl+' only comp.',outf)

        print(tp+' done')