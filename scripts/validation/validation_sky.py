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
parser.add_argument("--weight_col", help="column name for weight",default='WEIGHT_SN')
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
        os.makedirs(outdir)
        print('made '+outdir)


qt = 'COMP_TILE'

nside = args.nside
nest = True
zcol = 'Z_not4clus'
nran = 18

tps = [args.tracers]
if args.tracers == 'all':
    tps = ['ELG_LOPnotqso','LRG','QSO','BGS_BRIGHT-21.35']

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

def plot_map_healpix(nside, v, pix=None, vmin=None, vmax=None, cmap='jet', title=None, save_path=None,
             xsize=None, dpi=None, show=True, timing=False, nest=False, coord=None, cbar_label=''):
    from astropy.coordinates import SkyCoord
    from astropy import units
    import healpy as hp
    from healpy.newvisufunc import projview, newprojplot
    # Font sizes for healpix maps
    fontsize_dict = {
        "xlabel": 9.5,
        "ylabel": 9.5,
        "title": 9.5,
        "xtick_label": 9.5,
        "ytick_label": 9.5,
        "cbar_label": 9.5,
        "cbar_tick_label": 9.5,
    }
    params = {'legend.fontsize': 'x-large',
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large',
          'figure.facecolor': 'w'}
    plt.rcParams.update(params)

    default_dpi = {32: 100, 64: 200, 128: 400, 256: 600, 512: 1200}
    default_xsize = {32: 1500, 64: 4000, 128: 4000, 256: 6000, 512: 12000}

    
    if xsize is None:
        xsize = default_xsize[nside]

    if dpi is None:
        dpi = default_dpi[nside]

    npix = hp.nside2npix(nside)

    v = np.array(v)

    # Density map
    hp_mask = np.zeros(npix, dtype=bool)
    if pix is None:
        map_values = v.copy()
        hp_mask[np.isfinite(map_values)] = True
    else:
        map_values = np.zeros(npix, dtype=v.dtype)
        map_values[pix] = v
        hp_mask[pix] = True
    mplot = hp.ma(map_values)
    mplot.mask = ~hp_mask

    # Galactic plane
    org = 120
    tmpn = 1000
    cs = SkyCoord(l=np.linspace(0, 360, tmpn) * units.deg, b=np.zeros(tmpn) * units.deg, frame="galactic")
    ras, decs = cs.icrs.ra.degree, cs.icrs.dec.degree
    ras = np.remainder(ras + 360 - org, 360)  # shift ra values
    ras[ras > 180] -= 360  # scale conversion to [-180, 180]
    ii = ras.argsort()
    ras, decs = ras[ii], decs[ii]

    if timing:
        time_start = time.time()

    projview(mplot, min=vmin, max=vmax,
             rot=(120, 0, 0), coord=coord, cmap=cmap, xsize=xsize,
             graticule=True, graticule_labels=True, projection_type="mollweide", nest=nest,
             title=title,
             xlabel='RA (deg)', ylabel='Dec (deg)',
             custom_xtick_labels=[r'$240\degree$', r'$180\degree$', r'$120\degree$', r'$60\degree$', r'$0\degree$'],
             fontsize=fontsize_dict, unit=cbar_label)
    newprojplot(theta=np.radians(90-decs), phi=np.radians(ras), color='k', lw=1)
    if save_path is not None:
        plt.savefig(save_path, bbox_inches="tight", dpi=dpi)
    if show:
        plt.show()
    else:
        plt.close()

    if timing:
        print('Done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
    print('use DESI < 22.2 environment for a white background')


if args.survey == 'SV3' and args.tracers == 'all':
    tps = ['QSO','LRG','BGS_ANY','BGS_BRIGHT','ELG','ELG_HIP','ELG_HIPnotqso','ELGnotqso']
    zdw = ''
    if args.data != 'LSS':
        tps = ['QSO','LRG','ELG']
for tp in tps:
    tpr = tp
    prog = 'dark'
    if 'BGS_BRIGHT-' in tp:
        tpr = 'BGS_BRIGHT'
    if 'BGS' in tp:    
        prog = 'bright'
    dtfh = fitsio.read_header(indir+tpr+zdw+'_full_noveto.dat.fits',ext=1)
    for nr in range(0,nran):
        rffh = fitsio.read_header(indir+prog+'_'+str(nr)+'_full_noveto.ran.fits',ext=1)   
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
        odlhp = []
        odl_ochp = []
        dt = Table(fitsio.read(indir+tp+zdw+'_full'+args.use_map_veto+'.dat.fits'))
        cols = list(dt.dtype.names)
        sel_gz = common.goodz_infull(tp[:3],dt)
        sel_obs = dt['ZWARN'] != 999999
        dt = dt[sel_obs&sel_gz]
        dt['WEIGHT_COMP'] = 1./dt['FRACZ_TILELOCID']
        if 'FRAC_TLOBS_TILES' in cols and args.compmd == 'dat':
            dt['WEIGHT_COMP'] *= 1/dt['FRAC_TLOBS_TILES']
        cols = list(dt.dtype.names)
        if args.weight_col not in cols:
            print('no '+args.weight_col+', getting set to 1 for plotting')
            dt[args.weight_col] = np.ones(len(dt))
        if 'WEIGHT_ZFAIL' not in cols:
            dt['WEIGHT_ZFAIL'] = np.ones(len(dt))
            print('no redshift failure weight, getting set to 1 for plotting')
        dt['WEIGHT'] = dt['WEIGHT_COMP']*dt['WEIGHT_ZFAIL']*dt[args.weight_col]
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

        vm = 0
        vx = 1
        titl  = tp+' FRAC_TLOBS_TILES data'
        outf = outdir+tp+'_datafractlobs.png'
        plot_map_sindec(rada,sindd,dt['FRAC_TLOBS_TILES'],vm,vx,titl,outf)

        titl  = tp+' FRAC_TLOBS_TILES random'
        outf = outdir+tp+'_randfractlobs.png'
        plot_map_sindec(rara,sindr,rta['FRAC_TLOBS_TILES'],vm,vx,titl,outf)


        #vm = np.min(dt['WEIGHT_SYS'])
        #vx = np.max(dt['WEIGHT_SYS'])
        vm = 0.75
        vx = 1.25
        titl = tp+' '+args.weight_col
        outf = outdir+tp+'_'+args.weight_col+'.png'
        plot_map_sindec(rada,sindd,dt[args.weight_col],vm,vx,titl,outf)

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
            
            # healpix map data
            oddhp = np.zeros(len(rpix))
            oddhp[wp] = od[wp]
            odlhp.append(oddhp)
            odd_ochp = np.zeros(len(rpix))
            odd_ochp[wp] = od_oc[wp]            
            odl_ochp.append(odd_ochp)
            
            # sindec map data
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
        
        # sindec map data    
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

        
        # healpix map data
        odhp  = np.full(len(wp),np.inf)
        od_ochp  = np.full(len(wp),np.inf)
        odhp[odlhp[0]>0]   = odlhp[0][odlhp[0]>0]/np.mean(odlhp[0][odlhp[0]>0])
        odhp[odlhp[1]>0]   = odlhp[1][odlhp[1]>0]/np.mean(odlhp[1][odlhp[1]>0])
        
        od_ochp[odl_ochp[0]>0]= odl_ochp[0][odl_ochp[0]>0]/np.mean(odl_ochp[0][odl_ochp[0]>0])
        od_ochp[odl_ochp[1]>0]= odl_ochp[1][odl_ochp[1]>0]/np.mean(odl_ochp[1][odl_ochp[1]>0])

        outf = outdir+tp+f'_weighted_all_healpix{nside}.png'
        plot_map_healpix(nside,odhp,   vmin=vm,vmax=vx,title=titl+' for Year-1, including all corrections',   save_path=outf,nest=nest,cmap='viridis')
        
        outf = outdir+tp+f'_weighted_noimaging_healpix{nside}.png'
        plot_map_healpix(nside,od_ochp,vmin=vm,vmax=vx,title=titl+' for Year-1, without imaging corrections',save_path=outf,nest=nest,cmap='viridis')

        
        
        print(tp+' done')
