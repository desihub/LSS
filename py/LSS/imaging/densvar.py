import fitsio
import astropy.io.fits as fits
from astropy.table import Table
import healpy as hp
import numpy as np
from matplotlib import pyplot as plt
from LSS import imsys_fitter as sf

pixfn      = '/global/cfs/cdirs/desi/survey/catalogs/pixweight_maps_all/pixweight-1-dark.fits'#'/global/cfs/cdirs/desi/target/catalogs/dr9/0.57.0/pixweight/sv3/resolve/dark/sv3pixweight-1-dark.fits'
pixfn_ext      = '/global/cfs/cdirs/desi/survey/catalogs/pixweight_maps_all/pixweight_external.fits'
syst_names = ['STARDENS','EBV', 'PSFDEPTH_G', 'PSFDEPTH_R',\
    'PSFDEPTH_Z','GALDEPTH_G', 'GALDEPTH_R','GALDEPTH_Z',\
    'PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']
syst_names_ext = ['HALPHA','EBVreconMEANF15', 'CALIBG', 'CALIBR',\
    'CALIBZ']

hdr        = fits.getheader(pixfn,1)
nside,nest = hdr['HPXNSIDE'],hdr['HPXNEST']
print(nside,nest)

R_G=3.214 # http://legacysurvey.org/dr8/catalogs/#galactic-extinction-coefficients
R_R=2.165
R_Z=1.211

dr = '9'

#fidf = 'targetDR9m44.fits'
ranf = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits'
#ranf = '/global/cscratch1/sd/ajross/tarcat/vtest/tv0.49.0/randomsDR9v0.49.0_0_masked.fits'

 #these get used to veto imaging area; combination of bits applied to ELGs and LRGs in DR8 targeting

def get_prop_map(name):
    if name in syst_names:
        return fitsio.read(pixfn)[name]
    if name in syst_names_ext:
        return fitsio.read(pixfn_ext)[name]
    parsp = name.split('-')
    if len(parsp) > 1: 
        if parsp[1] == 'EBV':
            ebv = fitsio.read(pixfn)['EBV']
            if '_R' in par:
                R_v = R_R
            if '_G' in par:
                R_v = R_g
            if '_Z' in par:
                R_v = R_Z
            return  10.**(-0.4*R_v*ebv*2.)*fitsio.read(pixfn)[parsp[0]]
        elif parsp[1] == 'X' or parsp[1] == 'DIV':
            name1 = parsp[0]
            if name1 in syst_names:
                par1 = fitsio.read(pixfn)[name1]
            if name1 in syst_names_ext:
                par1 = fitsio.read(pixfn_ext)[name1]
            name2 = parsp[2]
            if name2 in syst_names:
                par2 = fitsio.read(pixfn)[name2]
            if name2 in syst_names_ext:
                par2 = fitsio.read(pixfn_ext)[name2]
            if parsp[1] == 'X':
                return par1*par2
            if parsp[1] == 'DIV':
                return par1/par2

    if name == 'PSFTOT':
        parv = fitsio.read(pixfn)
        return  (parv['PSFSIZE_G'])*(parv['PSFSIZE_R'])*(parv['PSFSIZE_Z'])
    elif name == 'SN2TOT_FLAT':
        parv = fitsio.read(pixfn)
        ebv = parv['EBV']
        return  10.**(-0.4*R_G*ebv*2.)*parv['PSFDEPTH_G'] + 10.**(-0.4*R_R*ebv*2.)*parv['PSFDEPTH_R'] + 10.**(-0.4*R_Z*ebv*2.)*parv['PSFDEPTH_Z']
    elif name == 'SNTOT_FLAT':
        parv = fitsio.read(pixfn)
        ebv = parv['EBV']
        mp = 10.**(-0.4*R_G*ebv*2.)*parv['PSFDEPTH_G'] + 10.**(-0.4*R_R*ebv*2.)*parv['PSFDEPTH_R'] + 10.**(-0.4*R_Z*ebv*2.)*parv['PSFDEPTH_Z']
        return  np.sqrt(mp)
    elif name == 'SN2TOT_G':
        parv = fitsio.read(pixfn)
        ebv = parv['EBV']
        return  10.**(-0.4*R_G*ebv*2.)*parv['PSFDEPTH_G']
    sys.exit('invalid map name')
    


def mask(dd,mb=[1]):
    keep = (dd['NOBS_G']>0) & (dd['NOBS_R']>0) & (dd['NOBS_Z']>0)
    print(len(dd[keep]))
    
    keepelg = keep
    for bit in mb:
        keepelg &= ((dd['MASKBITS'] & 2**bit)==0)
    print(len(dd[keepelg]))
    dd = dd[keepelg] 
    return dd       

def masklc(dd,mb=[1]):
    #keep = (dd['input_nobs_g']>0) & (dd['input_nobs_r']>0) & (dd['input_nobs_z']>0)
    keep = (dd['nobs_g']>0) & (dd['nobs_r']>0) & (dd['nobs_z']>0)
    print(len(dd[keep]))
    
    keepelg = keep
    for bit in mb:
        keepelg &= ((dd['maskbits'] & 2**bit)==0)
    print(len(dd[keepelg]))
    dd = dd[keepelg] 
    return dd       


def sel_reg(ra,dec,reg):
    wra = (ra > 100-dec)
    wra &= (ra < 280 +dec)
    if reg == 'DN':
        w = dec < 32.375
        w &= wra
    if reg == 'DS':
        w = ~wra
        w &= dec > -25
    return w        


def get_pix(nside, ra, dec, nest=0):
    return hp.ang2pix(nside, np.radians(-dec+90), np.radians(ra), nest=nest)

def add_par(dd,par):
    th,phi =radec2thphi(dd['RA'],dd['DEC'])
    hpxr = hp.ang2pix(nside,th,phi,nest=nest)
    gvp = []
    parl = fitsio.read(pixfn)[par]
    for px in hpxr:
        gvp.append(parl[px])
    dd[par+'_pix'] = np.array(gvp)
    return dd        

def plot_relnz_pixpar(sample,par,reg,zmin=0.8,zmax=1.6,nbin=8,nper = 5,survey='main',specrel='daily',version='test',basedir='/global/cfs/cdirs/desi/survey/catalogs/'):
    indir = basedir+survey+'/LSS/'+specrel+'/LSScats/'+version+'/'
    rcols = ['RA','DEC','PHOTSYS']
    zd = ''
    if survey != 'SV3':
        zd = 'zdone'
    rcol = ['RA','DEC','PHOTSYS']
    rd = fitsio.read(indir+sample+zd+'_0_full.ran.fits',columns=rcol)
    if reg == 'DN' or reg == 'DS':
        sel = sel_reg(rd['RA'],rd['DEC'],reg)
        rd = rd[sel]
    else: 
        iss = rd['PHOTSYS'] == 'S'
        if reg == 'S':
            sel = iss
        if reg == 'N':
            sel = ~iss
        rd = rd[sel]
    
    rd = Table(rd)
    rd = add_par(rd,par)

    if sample[:3] == 'ELG':
        dcols = ['RA','DEC','Z_not4clus','ZWARN','PHOTSYS','o2c','FRACZ_TILELOCID']

    dd = fitsio.read(indir+sample+zd+'_full.dat.fits',columns=dcols)
    reglab = ''
    if reg == 'DN' or reg == 'DS':
        sel = sel_reg(dd['RA'],dd['DEC'],reg)
        dd = dd[sel]
        reglab = 'DECaLS SGC'
        if reg == 'DN':
            reglab = 'DECaLS NGC'
    else:
        iss = dd['PHOTSYS'] == 'S'
        if reg == 'S':
            sel = iss
            reglab = 'DECaLS'
        if reg == 'N':
            sel = ~iss
            reglab = 'BASS/MzLS'
        dd = dd[sel]
    sel = dd['ZWARN'] != 999999
    if sample[:3] == 'ELG':
        sel &= dd['o2c'] > 0.9
    print(len(dd),len(dd[sel]))
    dd = dd[sel]
    dd = Table(dd)
    dd = add_par(dd,par)
    perl = []
    div = 100/nper
    for i in range(nper+1):
        perl.append(div*i)
    print('percentiles are '+str(perl))
    gdp = []
    for per in perl:
        gdp.append(np.percentile(dd[par+'_pix'],per))
    print(gdp)
    #cl = ['b','r','k','purple','brown']
    dndz_ot,be = np.histogram(dd['Z_not4clus'],range=(zmin,zmax),bins=nbin,weights=1/dd['FRACZ_TILELOCID'])
    bs = (zmax-zmin)/nbin
    fac = len(rd)
    for i in range(0,nper):
        sd = dd[par+'_pix'] >gdp[i]
        sd &= dd[par+'_pix'] <gdp[i+1]
        dndz_ob,be = np.histogram(dd[sd]['Z_not4clus'],range=(zmin,zmax),bins=nbin,weights=1/dd[sd]['FRACZ_TILELOCID'])    
        #plt.plot(be[:-1]+0.05,dndz_ob/dndz_ot,'--',color=cl[i])
        sdi = rd[par+'_pix'] >gdp[i]
        sdi &= rd[par+'_pix'] <gdp[i+1]
        facb = len(rd[sdi])#/len(obi_sel[sd])
        print(facb)
        #print(fac,facb,len(obi_sel[sd]))
        #plt.plot(be[:-1]+0.05,dndz_db/dndz_dt,':',color=cl[i])
        #plt.plot(be[:-1]+0.05,dndz_ob/dndz_ib*facb,label=str(i))
        plt.plot(be[:-1]+bs/2.,(dndz_ob*fac/facb-dndz_ot)/dndz_ot,label=str(round(gdp[i],3))+'<'+par+'<'+str(round(gdp[i+1],3)))
    plt.legend()
    ol = np.zeros(len(be[:-1]))
    plt.plot(be[:-1]+bs/2.,ol,':')
    plt.xlabel('redshift')
    plt.ylabel('relative change in n(z)')
    plt.title(survey+' '+sample+' '+reglab)
    plt.grid(True)
    plt.show()

def plot_relnz_clus_pixpar(sample,par,reg,weightcol='WEIGHT',zmin=0.8,zmax=1.6,nbin=8,nper = 5,survey='DA02',specrel='guadalupe',version='test',basedir='/global/cfs/cdirs/desi/survey/catalogs/'):
    indir = basedir+survey+'/LSS/'+specrel+'/LSScats/'+version+'/'
    zd = ''
    if survey != 'SV3':
        zd = 'zdone'


    rd = fitsio.read(indir+sample+zd+reg+'_0_clustering.ran.fits')
    
    rd = Table(rd)
    rd = add_par(rd,par)


    dd = fitsio.read(indir+sample+zd+reg+'_clustering.dat.fits')
    reglab = ''
    if reg == '_DN' or reg == '_DS':
        reglab = 'DECaLS SGC'
        if reg == '_DN':
            reglab = 'DECaLS NGC'
    else:
        if reg == '_S':
            reglab = 'DECaLS'
        if reg == '_N':
            reglab = 'BASS/MzLS'
    dd = Table(dd)
    dd = add_par(dd,par)
    if 'RF' in weightcol:
        weights = np.ones(len(dd))
        weights *= dd['WEIGHT_RF']*dd['WEIGHT_COMP']
        dd[weightcol] = weights
    perl = []
    div = 100/nper
    for i in range(nper+1):
        perl.append(div*i)
    print('percentiles are '+str(perl))
    gdp = []
    for per in perl:
        gdp.append(np.percentile(dd[par+'_pix'],per))
    print(gdp)

    roundfac = int(np.log10(gdp[-1]))
    roundfac = 2-roundfac
    if roundfac < 0:
        roundfac = None
    print(par,roundfac)

    #cl = ['b','r','k','purple','brown']
    dndz_ot,be = np.histogram(dd['Z'],range=(zmin,zmax),bins=nbin,weights=dd[weightcol])
    dndz_nt,_ = np.histogram(dd['Z'],range=(zmin,zmax),bins=nbin)
    perr = np.sqrt(dndz_nt)/dndz_nt*np.sqrt(nper)
    
    bs = (zmax-zmin)/nbin
    plt.fill_between(be[:-1]+bs/2.,-perr,perr,color='gray')
    fac = len(rd)
    for i in range(0,nper):
        sd = dd[par+'_pix'] >gdp[i]
        sd &= dd[par+'_pix'] <gdp[i+1]
        dndz_ob,be = np.histogram(dd[sd]['Z'],range=(zmin,zmax),bins=nbin,weights=dd[sd][weightcol])    
        #plt.plot(be[:-1]+0.05,dndz_ob/dndz_ot,'--',color=cl[i])
        sdi = rd[par+'_pix'] >gdp[i]
        sdi &= rd[par+'_pix'] <gdp[i+1]
        facb = len(rd[sdi])#/len(obi_sel[sd])
        print(facb)
        #print(fac,facb,len(obi_sel[sd]))
        #plt.plot(be[:-1]+0.05,dndz_db/dndz_dt,':',color=cl[i])
        #plt.plot(be[:-1]+0.05,dndz_ob/dndz_ib*facb,label=str(i))
        plt.plot(be[:-1]+bs/2.,(dndz_ob*fac/facb-dndz_ot)/dndz_ot,label=str(round(gdp[i],roundfac))+'<'+par+'<'+str(round(gdp[i+1],roundfac)))
    plt.legend()
    ol = np.zeros(len(be[:-1]))
    plt.plot(be[:-1]+bs/2.,ol,':')
    plt.xlabel('redshift')
    plt.ylabel('relative change in n(z)')
    plt.title(survey+' '+sample+' '+reglab)
    plt.grid(True)
    plt.ylim(-.2,.2)
    plt.show()


def read_systematic_maps(data_ra, data_dec, rand_ra, rand_dec,sys_tab=None):
    
    #-- Dictionaries containing all different systematic values
    data_syst = {}
    rand_syst = {}

    pixmext = None
    if sys_tab is None:
        pixm = fitsio.read(pixfn)
        pixmext = fitsio.read(pixfn_ext)
        use_syst_names = syst_names
    else:
        pixm = sys_tab
        use_syst_names = list(pixm.dtype.names)
    data_pix = get_pix(nside, data_ra, data_dec,nest) 
    rand_pix = get_pix(nside, rand_ra, rand_dec,nest)

    for syst_name in use_syst_names:
        data_syst[syst_name] = pixm[syst_name][data_pix]
        rand_syst[syst_name] = pixm[syst_name][rand_pix]

    if pixmext is not None:
        sel = pixmext['EBVreconMEANF15'] < -1
        mno = np.mean(pixmext[~sel]['EBVreconMEANF15'])
        #pixmext[sel]['EBVreconMEANF15'] = mno #why did this not work !?
        for i in range(0,len(pixmext['EBVreconMEANF15'])):
            if pixmext['EBVreconMEANF15'][i] < -1:
                pixmext['EBVreconMEANF15'][i] = mno
        print(np.min(pixmext['EBVreconMEANF15']),np.min(pixmext['EBVreconMEANF15'][sel]),np.min(pixmext['EBVreconMEANF15'][~sel]))
        for syst_name in syst_names_ext:
            data_syst[syst_name] = pixmext[syst_name][data_pix]
            rand_syst[syst_name] = pixmext[syst_name][rand_pix]

        ebvd = pixm['EBV'] - pixmext['EBVreconMEANF15']
        data_syst['DELTA_EBV'] = ebvd[data_pix]
        rand_syst['DELTA_EBV'] = ebvd[rand_pix]

    return data_syst, rand_syst


def get_imweight(dd,rd,zmin,zmax,reg,fit_maps,use_maps,plotr=True,zcol='Z',sys_tab=None,wtmd='fracz',figname='temp.png',modoutname='temp.txt',logger=None):
    import LSS.common_tools as common
    sel = dd[zcol] > zmin
    sel &= dd[zcol] < zmax
    if reg == 'N' or reg == 'S':
        sel &= dd['PHOTSYS'] == reg
        selr = rd['PHOTSYS'] == reg
    elif 'DES' in reg:
        inDES = common.select_regressis_DES(dd)
        inDESr = common.select_regressis_DES(rd)
        if reg == 'DES':
            sel &= inDES
            selr = inDESr
        if reg == 'SnotDES':
            sel &= dd['PHOTSYS'] == 'S'
            sel &= ~inDES
            selr = rd['PHOTSYS'] == 'S'
            selr &= ~inDESr

    else:
        print('other regions not currently supported')
        return 'Exiting due to critical error with region'
 
    dds = dd[sel]

    if wtmd == 'clus':
        selr &= rd[zcol] > zmin
        selr &= rd[zcol] < zmax

    rd = rd[selr]

    #-- Dictionaries containing all different systematic values
    data_syst, rand_syst = read_systematic_maps(dds['RA'],dds['DEC'],rd['RA'],rd['DEC'],sys_tab=sys_tab)
    #print(data_syst.keys)
    cols = list(dd.dtype.names)
    weights_ran = np.ones(len(rd))
    if wtmd == 'fracz':
        common.printlog('using 1/FRACZ_TILELOCID based completeness weights',logger)
        wts = 1/dds['FRACZ_TILELOCID']
        if 'FRAC_TLOBS_TILES' in cols:
            common.printlog('using FRAC_TLOBS_TILES',logger)
            wts *= 1/dds['FRAC_TLOBS_TILES']
    if wtmd == 'wt':
        wts = dds['WEIGHT']
        weights_ran = rd['WEIGHT']
    if wtmd == 'wtfkp':
        wts = dds['WEIGHT']*dds['WEIGHT_FKP']
        weights_ran = rd['WEIGHT']*rd['WEIGHT_FKP']

    if wtmd == 'clus':
        wts = dds['WEIGHT']*dds['WEIGHT_FKP']/dds['WEIGHT_SYS']
        weights_ran = rd['WEIGHT']*rd['WEIGHT_FKP']/rd['WEIGHT_SYS']
    if wtmd == 'wt_comp':
        wts = dds['WEIGHT_COMP']

    if 'WEIGHT_ZFAIL' in cols and wtmd != 'clus':
        wts *= dds['WEIGHT_ZFAIL']

    data_we = wts
    rand_we = weights_ran
    #-- Create fitter object
    s = sf.Syst(data_we, rand_we)

    fit_maps = fit_maps
    use_maps = use_maps

    #-- Add the systematic maps we want 
    for syst_name in use_maps:
        s.add_syst(syst_name, data_syst[syst_name], rand_syst[syst_name])
    s.cut_outliers(p=0.5, verbose=1) 

    #-- Perform global fit
    nbins=10
    s.prepare(nbins=nbins)
    #for name in fit_maps:
    #    print(name,len(s.data_syst[name]),len(s.data_we))
    s.fit_minuit(fit_maps=fit_maps)
    common.printlog(str(s.best_pars),logger)
    common.printlog(str(list(s.best_pars)),logger)
    pars_dict = {}
    common.printlog('writing to '+modoutname,logger)
    fo = open(modoutname,'w')
    for par_name, p in zip(s.par_names, list(s.best_pars)):
        pars_dict[par_name] = p
        fo.write(str(par_name)+' '+str(p)+'\n')
    fo.close()
    if plotr:
        #s.plot_overdensity(pars=[None, s.best_pars], ylim=[0.7, 1.3])#, title=f'{sample_name}: global fit')
        common.printlog('saving figure to '+figname)
        s.plot_overdensity(pars=[None, pars_dict], ylim=[0.7, 1.3])
        plt.savefig(figname)
        plt.clf()
        #plt.show()
    #-- Get weights for global fit
    #data_weightsys_global = 1/s.get_model(s.best_pars, data_syst)
    data_weightsys_global = 1/s.get_model(pars_dict, data_syst)
    wsysl = np.ones(len(dd))
    wsysl[sel] = data_weightsys_global 
    del data_syst
    del rand_syst
    return wsysl


def radec2thphi(ra,dec):
    return (-dec+90.)*np.pi/180.,ra*np.pi/180.
    
def thphi2radec(theta,phi):
    return 180./np.pi*phi,-(180./np.pi*theta-90)

def obiELGvspar(reg,par,vmin=None,vmax=None,nbin=10,obidir='/global/cscratch1/sd/adematti/legacysim/dr9/ebv1000shaper/',elgandlrgbits = [1,5,6,7,8,9,11,12,13]):
    from desitarget import cuts
    NS = None
    if reg == 'N':
        NS = 'north'
        south = False
    if reg == 'DS' or reg == 'DN':
        NS = 'south' 
        south = True
    if NS == None:
        print('!!!NS not set, you must have chose an invalid reg; should be N, DS, or DN!!!')
        return(None)
    print(NS)
    obif = fitsio.read(obidir+NS+'/file0_rs0_skip0/merged/matched_input.fits')
    print(len(obif))
    obi_masked = masklc(obif,mb=elgandlrgbits)
    
    if reg != 'N':
        wr = sel_reg(obi_masked['ra'],obi_masked['dec'],reg)
        obi_masked = obi_masked[wr]    
    gflux = obi_masked['flux_g']/obi_masked['mw_transmission_g']
    rflux = obi_masked['flux_r']/obi_masked['mw_transmission_r']
    zflux = obi_masked['flux_z']/obi_masked['mw_transmission_z']
    ws = cuts.isELG_colors(gflux, rflux, zflux,south=south)
    print(len(obi_masked[ws])) 
    ws &= (obi_masked['ra']*0 == 0)
    print(len(obi_masked[ws])) 
    obi_elg = obi_masked[ws]             

    if vmin is None:
        vmin = np.min(rl[par])
    if vmax is None:
        vmax = np.max(rl[par])    
        
    rh,bn = np.histogram(obi_masked[par],bins=nbin,range=(vmin,vmax))
    dh,db = np.histogram(obi_elg[par],bins=bn)
    rf = len(obi_masked)/len(obi_elg)
    sv = dh/rh*rf
    ep = np.sqrt(dh)/rh*rf
    bc = []
    for i in range(0,len(bn)-1):
        bc.append((bn[i]+bn[i+1])/2.)
    plt.errorbar(bc,sv-1.,ep,fmt='ko')
    plt.hist(obi_masked[par],bins=nbin,range=(vmin,vmax),weights=0.2*np.ones(len(obi_masked))/np.max(rh))
    plt.ylim(-.3,.3)
    plt.xlabel(par)
    plt.ylabel('Ngal/<Ngal> - 1')
    plt.title('Obiwan ELGs in '+reg + ' footprint')
    plt.show()
    wv = (obi_masked[par]>vmin) & (obi_masked[par] < vmax)
    frac = len(obi_masked[~wv])/len(obi_masked)
    print('fraction of randoms not included in plot: '+str(frac))
    return bc,sv,ep

def obiLRGvspar(reg,par,vmin=None,vmax=None,syspix=False,md='sv3',nbin=10,obidir='/global/cfs/cdirs/desi/users/huikong/decals_ngc/subset/',elgandlrgbits = [1,5,6,7,8,9,11,12,13]):
    if md == 'sv3':
        from desitarget.sv3 import sv3_cuts as cuts
    else:
        from desitarget import cuts
    NS = None
    if reg == 'N':
        NS = 'north'
        south = False
    if reg == 'DS' or reg == 'DN' or reg =='S':
        NS = 'south' 
        south = True
    if NS == None:
        print('!!!NS not set, you must have chose an invalid reg; should be N, DS, or DN!!!')
        return(None)
    print(NS)
    obif = fitsio.read(obidir+'/subset_rs0.fits')
    sel = obif['matched'] == True
    sel &= obif['dec'] < 32.375
    obif = obif[sel]
    print(len(obif))
    obi_masked = masklc(obif,mb=elgandlrgbits)
    print(len(obi_masked))
    gflux = obi_masked['flux_g']/obi_masked['mw_transmission_g']
    rflux = obi_masked['flux_r']/obi_masked['mw_transmission_r']
    zflux = obi_masked['flux_z']/obi_masked['mw_transmission_z']
    w1flux = obi_masked['flux_w1']/obi_masked['mw_transmission_w1']
    zfiberflux = obi_masked['fiberflux_z']/obi_masked['mw_transmission_z']
    
    ws = cuts.isLRG_colors(gflux, rflux, zflux, w1flux,
                 zfiberflux,  south=south)
    
    if md == 'sv3':
        ws = ws[0]
    print(len(obi_masked[ws])) 
    ws &= (obi_masked['ra']*0 == 0)
    print(len(obi_masked[ws])) 
    obi_lrg = obi_masked[ws]             

        
    #if syspix:
    rth,rphi = radec2thphi(obi_masked['ra'],obi_masked['dec'])
    rpix = hp.ang2pix(nside,rth,rphi,nest=nest)
    pixlr = np.zeros(12*nside*nside)
    for pix in rpix:
        pixlr[pix] += 1.
    dth,dphi = radec2thphi(obi_lrg['ra'],obi_lrg['dec'])
    dpix = hp.ang2pix(nside,dth,dphi,nest=nest)
    pixld = np.zeros(12*nside*nside)
    for pix in dpix:
        pixld[pix] += 1.
    parv = fitsio.read(pixfn)[par.upper()]
    wp = pixlr > 0
    vmin = np.percentile(parv[wp],1)  
    vmax = np.percentile(parv[wp],99)  
    bcp,svp,epp = plot_pixdens1d(pixld[wp],pixlr[wp],parv[wp],vmin=vmin,vmax=vmax)   
    
    datf = fitsio.read(obidir+'subset_dr9_lrg_sv3.fits') 
    datf = masklc(datf)
    sel = datf['dec'] < 32.375
    datf = datf[sel]
    ranf = fitsio.read(obidir+'subset_random.fits')
    ranf = masklc(ranf)
    sel = ranf['dec'] < 32.375
    ranf = ranf[sel]

    rth,rphi = radec2thphi(ranf['ra'],ranf['dec'])
    rpix = hp.ang2pix(nside,rth,rphi,nest=nest)
    pixlr = np.zeros(12*nside*nside)
    for pix in rpix:
        pixlr[pix] += 1.
    dth,dphi = radec2thphi(datf['ra'],datf['dec'])
    dpix = hp.ang2pix(nside,dth,dphi,nest=nest)
    pixld = np.zeros(12*nside*nside)
    for pix in dpix:
        pixld[pix] += 1.
    parv = fitsio.read(pixfn)[par.upper()]
    wp = pixlr > 0
        
    bcpd,svpd,eppd = plot_pixdens1d(pixld[wp],pixlr[wp],parv[wp],vmin=vmin,vmax=vmax)   
    


    #else:
    if vmin is None:
        #vmin = np.min(obi_lrg[par])
        vmin = np.percentile(obi_lrg[par],1)
    if vmax is None:
        #vmax = np.max(obi_lrg[par]) 
        vmax = np.percentile(obi_lrg[par],99)   

    rh,bn = np.histogram(obi_masked[par],bins=nbin,range=(vmin,vmax))
    dh,db = np.histogram(obi_lrg[par],bins=bn)
    rf = len(obi_masked)/len(obi_lrg)
    sv = dh/rh*rf
    ep = np.sqrt(dh)/rh*rf
    rh,bn = np.histogram(ranf[par],bins=nbin,range=(vmin,vmax))
    dh,db = np.histogram(datf[par],bins=bn)
    rf = len(ranf)/len(datf)
    svd = dh/rh*rf
    epd = np.sqrt(dh)/rh*rf
    bc = []
    for i in range(0,len(bn)-1):
        bc.append((bn[i]+bn[i+1])/2.)
    plt.errorbar(bc,sv-1.,ep,fmt='r-',label='obiwan points')
    plt.plot(bc,svd-1.,'k-',label='data points')
    plt.plot(bcpd,svpd-1.,'k--',label='data pixelized')
    plt.errorbar(bcp,svp-1.,epp,fmt='r--',label='obiwan pixelized')
    plt.legend()
    plt.title('SV3 selection ')
    #plt.show()

    #plt.title('SV3 selection points')
    #plt.show()
   
    #plt.hist(obi_masked[par],bins=nbin,range=(vmin,vmax),weights=0.2*np.ones(len(obi_masked))/np.max(rh))
    #plt.ylim(-.3,.3)
    plt.xlabel(par)
    plt.ylabel('Ngal/<Ngal> - 1')
    #plt.title('Obiwan LRGs in '+reg + ' footprint')
    wv = (obi_masked[par]>vmin) & (obi_masked[par] < vmax)
    frac = len(obi_masked[~wv])/len(obi_masked)
    print('fraction of randoms not included in plot: '+str(frac))

    plt.show()
    return bc,sv,ep

def obiLRGvs_depthmag(reg,par,band,vmin=None,vmax=None,syspix=False,md='sv3',nbin=10,obidir='/global/cfs/cdirs/desi/users/huikong/decals_ngc/subset/',elgandlrgbits = [1,5,6,7,8,9,11,12,13]):
    if md == 'sv3':
        from desitarget.sv3 import sv3_cuts as cuts
    else:
        from desitarget import cuts
    NS = None
    if reg == 'N':
        NS = 'north'
        south = False
    if reg == 'DS' or reg == 'DN' or reg =='S':
        NS = 'south' 
        south = True
    if NS == None:
        print('!!!NS not set, you must have chose an invalid reg; should be N, DS, or DN!!!')
        return(None)
    print(NS)
    obif = fitsio.read(obidir+'/subset_rs0.fits')
    sel = obif['matched'] == True
    sel &= obif['dec'] < 32.375
    sel &= obif['dec'] > -10
    obif = obif[sel]
    print(len(obif))
    obi_masked = masklc(obif,mb=elgandlrgbits)
    print(len(obi_masked))
    gflux = obi_masked['flux_g']/obi_masked['mw_transmission_g']
    rflux = obi_masked['flux_r']/obi_masked['mw_transmission_r']
    zflux = obi_masked['flux_z']/obi_masked['mw_transmission_z']
    w1flux = obi_masked['flux_w1']/obi_masked['mw_transmission_w1']
    zfiberflux = obi_masked['fiberflux_z']/obi_masked['mw_transmission_z']
    
    ws = cuts.isLRG_colors(gflux, rflux, zflux, w1flux,
                 zfiberflux,  south=south)
    
    if md == 'sv3':
        ws = ws[0]
    print(len(obi_masked[ws])) 
    ws &= (obi_masked['ra']*0 == 0)
    print(len(obi_masked[ws])) 
    obi_lrg = obi_masked[ws]             

    par = par+'_'+band
    if band == 'g':
        ec = R_G   
    if band == 'r':
        ec = R_R
    if band ==  'z':
        ec = R_Z     
    #if syspix:
    rth,rphi = radec2thphi(obi_masked['ra'],obi_masked['dec'])
    rpix = hp.ang2pix(nside,rth,rphi,nest=nest)
    pixlr = np.zeros(12*nside*nside)
    for pix in rpix:
        pixlr[pix] += 1.
    dth,dphi = radec2thphi(obi_lrg['ra'],obi_lrg['dec'])
    dpix = hp.ang2pix(nside,dth,dphi,nest=nest)
    pixld = np.zeros(12*nside*nside)
    for pix in dpix:
        pixld[pix] += 1.
    parv = fitsio.read(pixfn)[par.upper()]
    parv = -2.5*(np.log10(5/np.sqrt(parv))-9)-ec*fitsio.read(pixfn)['EBV']
    wp = pixlr > 0
    wp &= parv > 0
    wp &= parv*0 == 0
    if vmin == None:
        vmin = np.percentile(parv[wp],2)  
    if vmax == None:
        vmax = np.percentile(parv[wp],99.5)  
    bcp,svp,epp = plot_pixdens1d(pixld[wp],pixlr[wp],parv[wp],vmin=vmin,vmax=vmax)   
    
    datf = fitsio.read(obidir+'subset_dr9_lrg_sv3.fits') 
    datf = masklc(datf,mb=elgandlrgbits)
    sel = datf['dec'] < 32.375
    sel &= datf['dec'] > -10
    datf = datf[sel]
    ranf = fitsio.read(obidir+'subset_random.fits')
    ranf = masklc(ranf,mb=elgandlrgbits)
    sel = ranf['dec'] < 32.375
    sel &= ranf['dec'] > -10
    ranf = ranf[sel]

    rth,rphi = radec2thphi(ranf['ra'],ranf['dec'])
    rpix = hp.ang2pix(nside,rth,rphi,nest=nest)
    pixlr = np.zeros(12*nside*nside)
    for pix in rpix:
        pixlr[pix] += 1.
    dth,dphi = radec2thphi(datf['ra'],datf['dec'])
    dpix = hp.ang2pix(nside,dth,dphi,nest=nest)
    pixld = np.zeros(12*nside*nside)
    for pix in dpix:
        pixld[pix] += 1.
    
        
    bcpd,svpd,eppd = plot_pixdens1d(pixld[wp],pixlr[wp],parv[wp],vmin=vmin,vmax=vmax)   
    


    #else:

    omp = -2.5*(np.log10(5/np.sqrt(obi_masked[par]))-9)-ec*obi_masked['ebv']
    rh,bn = np.histogram(omp,bins=nbin,range=(vmin,vmax))
    olp = -2.5*(np.log10(5/np.sqrt(obi_lrg[par]))-9)-ec*obi_lrg['ebv']
    dh,db = np.histogram(olp,bins=bn)
    rf = len(obi_masked)/len(obi_lrg)
    sv = dh/rh*rf
    ep = np.sqrt(dh)/rh*rf
    rp = -2.5*(np.log10(5/np.sqrt(ranf[par]))-9)-ec*ranf['ebv']
    rh,bn = np.histogram(rp,bins=nbin,range=(vmin,vmax))
    dp = -2.5*(np.log10(5/np.sqrt(datf[par]))-9)-ec*datf['ebv']
    dh,db = np.histogram(dp,bins=bn)
    rf = len(ranf)/len(datf)
    svd = dh/rh*rf
    epd = np.sqrt(dh)/rh*rf
    bc = []
    for i in range(0,len(bn)-1):
        bc.append((bn[i]+bn[i+1])/2.)
    plt.errorbar(bc,sv-1.,ep,fmt='r-',label='obiwan points')
    plt.plot(bc,svd-1.,'k-',label='data points')
    plt.plot(bcpd,svpd-1.,'k--',label='data pixelized')
    plt.errorbar(bcp,svp-1.,epp,fmt='r--',label='obiwan pixelized')
    plt.legend()
    plt.title('SV3 selection ')
    #plt.show()

    #plt.title('SV3 selection points')
    #plt.show()
   
    #plt.hist(obi_masked[par],bins=nbin,range=(vmin,vmax),weights=0.2*np.ones(len(obi_masked))/np.max(rh))
    #plt.ylim(-.3,.3)
    plt.xlabel(par+ ' magnitudes (ext corrected)')
    plt.ylabel('Ngal/<Ngal> - 1')
    #plt.title('Obiwan LRGs in '+reg + ' footprint')
    wv = (omp>vmin) & (omp < vmax)
    frac = len(obi_masked[~wv])/len(obi_masked)
    print('fraction of randoms not included in plot: '+str(frac))

    plt.show()
    plt.hist(parv[wp],histtype='step',density=True,bins=30)
    plt.hist(rp,histtype='step',density=True,bins=30)
    plt.show()
    return bc,sv,ep



def gethpmap(dl,reg=False,weights=None):
    if reg:
        if reg == 'S' or reg == 'N':
            wr = dl['PHOTSYS'] == reg
        else:
            wr = sel_reg(dl['RA'],dl['DEC'],reg)
        dl = dl[wr]
    rth,rphi = radec2thphi(dl['RA'],dl['DEC'])
    rpix = hp.ang2pix(nside,rth,rphi,nest=nest)
    wts = np.ones(len(rth))
    if weights is not None:
        wts = dl[weights]
    pixlr = np.zeros(12*nside*nside)
    for pix,wt in zip(rpix,wts):
        
        pixlr[pix] += wt
    return pixlr

def gethpmap_var(dl,reg=False):
    if reg:
        if reg == 'S' or reg == 'N':
            wr = dl['PHOTSYS'] == reg
        else:
            wr = sel_reg(dl['RA'],dl['DEC'],reg)
        dl = dl[wr]
    rth,rphi = radec2thphi(dl['RA'],dl['DEC'])
    rpix = hp.ang2pix(nside,rth,rphi,nest=nest)
    pixlp = np.zeros(12*nside*nside)
    pixlv = np.zeros(12*nside*nside)
    for i in range(0,len(rpix)): 
        pix = rpix[i]
        pixlp[pix] += dl[i][par.split('-')[1]]
        pixlv[pix] += dl[i][par.split('-')[1]]**2.
    return pixlp,pixlv


def plot_hpmap(wp,od,reg=False,sz=.1,vx=1.5,vm=.5,titl=''):
    pixls = np.arange(12*nside*nside,dtype=int)
    th,phi = hp.pix2ang(nside,pixls[wp],nest=nest)
    ra,dec = thphi2radec(th,phi)
    if reg == 'DS':
        wr = ra > 250
        ra[wr] -=360
    if vx == None:
        vx = np.max(od)
    if vm == None:
        vm = np.min(od)    

    plt.scatter(ra,np.sin(dec*np.pi/180),c=od,s=sz,edgecolor='none',vmax=vx,vmin=vm)#,vmin=1.,vmax=2)
    plt.xlabel('RA')
    plt.ylabel('sin(DEC)')
    plt.colorbar()
    plt.title(titl)
    plt.show()
   

def plot_hpdens(rl,ft,reg=False,fnc=None,sz=.2,vx=1.5,vm=.5,datweights=None,weights=None,wsel=None,titl=''):
    pixlr = gethpmap(rl,reg)
    print('randoms done')
    pixlg = gethpmap(ft,reg,datweights)
    print('data done')
    
    if weights is None:
        weights = np.ones(len(pixlr))
    if wsel is not None:
        wp = wsel
        wp &= (pixlr > 0)
    else:
        wp = (pixlr > 0) 
    wp &= (weights*0 == 0)
    wp &= pixlg*0 == 0
    od = pixlg[wp]/pixlr[wp]*weights[wp]
    od = od/np.mean(od)
    plot_hpmap(wp,od,reg,sz,vx,vm,titl)

def get_hpdens(rl,ft,reg=False,fnc=None,sz=.2,vx=1.5,vm=.5,datweights=None,weights=None,wsel=None,titl=''):
    pixlr = gethpmap(rl,reg)
    print('randoms done')
    pixlg = gethpmap(ft,reg,datweights)
    print('data done')
    
    if weights is None:
        weights = np.ones(len(pixlr))
    if wsel is not None:
        wp = wsel
        wp &= (pixlr > 0)
    else:
        wp = (pixlr > 0) 
    wp &= (weights*0 == 0)
    wp &= pixlg*0 == 0
    od = pixlg[wp]/pixlr[wp]*weights[wp]
    od = od/np.mean(od)
    return wp,od


def plot_hpprop(rl,par,reg=False,fnc=None,sz=.2,vx=None,vm=None,weights=None):
    pixlr = gethpmap(rl,reg)
    print('randoms done')

    if weights is None:
        weights = np.ones(len(pixlr))
    wp = (pixlr > 0) 
    wp &= (weights*0 == 0)
    parv = fitsio.read(pixfn)
    if par == 'PSFTOT':
        parv = (parv[wp]['PSFSIZE_G'])*(parv[wp]['PSFSIZE_R'])*(parv[wp]['PSFSIZE_Z'])
    elif par == 'SN2TOT_FLAT':
        ebv = parv[wp]['EBV']
        parv = 10.**(-0.4*R_G*ebv*2.)*parv[wp]['PSFDEPTH_G'] + 10.**(-0.4*R_R*ebv*2.)*parv[wp]['PSFDEPTH_R'] + 10.**(-0.4*R_Z*ebv*2.)*parv[wp]['PSFDEPTH_Z']

    elif par == 'fracPSF':
        wpsf = ft['MORPHTYPE'] == 'PSF'
        pixlgp = np.zeros(12*nside*nside)
        dpixp = dpix[wpsf]
        for i in range(0,len(dpixp)): 
            pix = dpixp[i]
            pixlgp[pix] += 1.
        parv = pixlgp[wp]/pixlg[wp]

    else:    
        parv = parv[wp][par]
    od = parv
    plot_hpmap(wp,od,reg,sz,vx,vm,titl=par)

def densvsinput_pix(rl,ft,parl,xlab='',wsel=None,reg=None,fnc=None,vmin=None,vmax=None,ebvcut=None,edscut=None,sn2cut=None,fpsfcut=None,gfluxcut=None,rfluxcut=None,gbcut=None,nbin=10,weights=None,titl=''):        
    pixlr = gethpmap(rl,reg)
    print('randoms done')
    pixlg = gethpmap(ft,reg)
    print('data done')

    if weights is None:
        weights = np.ones(len(pixlr))

    if wsel is not None:
        wp = wsel
        wp &= (pixlr > 0) 
    else:
        wp = (pixlr > 0) 
    wp &= (weights*0 == 0)

    parv = fitsio.read(pixfn)
    ebv = parv['EBV']
    sn2tf = 10.**(-0.4*R_G*ebv*2.)*parv['PSFDEPTH_G'] + 10.**(-0.4*R_R*ebv*2.)*parv['PSFDEPTH_R'] + 10.**(-0.4*R_Z*ebv*2.)*parv['PSFDEPTH_Z']
    print(len(parv[wp]))
    if sn2cut:
        wp &= (sn2tf > sn2cut)
    
    if fpsfcut:
        wpsf = ft['MORPHTYPE'] == 'PSF'
        pixlgp = np.zeros(12*nside*nside)
        dpixp = dpix[wpsf]
        for i in range(0,len(dpixp)): 
            pix = dpixp[i]
            pixlgp[pix] += 1.
        fpsf = pixlgp/pixlg
        wp &= (fpsf < fpsfcut)
    if ebvcut:
        wp &= (parv['EBV'] < ebvcut)

    if edscut:
        eds = parv['EBV']/parv['STARDENS']
        wp &= (eds < edscut)

    parv = parl 

    wp &= parv !=0
    wp &= parv*0 == 0
    print(len(parv[wp]))
    bc,sv,ep = plot_pixdens1d(pixlg[wp],pixlr[wp],parv[wp],weights[wp],vmin,vmax,titl=titl,xlab=xlab)
    return bc,sv,ep


def plot_pixdens1d(pixlg,pixlr,parv,weights=None,vmin=None,vmax=None,smean=True,addhist=True,rng=0.3,titl='',nbin=10,xlab=''):
    if vmin is None:
        vmin = np.min(parv)
    if vmax is None:
        vmax = np.max(parv)
    if weights is None:
        weights = np.ones(len(pixlg))
    rh,bn = np.histogram(parv,bins=nbin,range=(vmin,vmax),weights=pixlr)
    dh,db = np.histogram(parv,bins=bn,weights=pixlg*weights)
    norm = sum(rh)/sum(dh)
    sv = dh/rh*norm
    ep = np.sqrt(dh)/rh*norm
    bc = []
    for i in range(0,len(bn)-1):
        bc.append((bn[i]+bn[i+1])/2.)
    sb = 0
    if smean:
        sb = 1
        plt.ylabel('Ngal/<Ngal> - 1')
    else:
        plt.ylabel('Ngal/<Ngal> ')    
    plt.errorbar(bc,sv-sb,ep,fmt='ko')
    if addhist:
        plt.hist(parv,bins=nbin,range=(vmin,vmax),weights=pixlr*0.66*rng*np.ones(len(pixlr))/np.max(rh))
    plt.ylim(1-rng-sb,1+rng-sb)
    plt.xlabel(xlab)
    
    plt.title(titl)
    coeff = np.polyfit(bc,sv,1,w=1/ep)
    plt.plot(bc,coeff[1]+np.array(bc)*coeff[0]-sb,'r--')
    print(coeff)
    plt.show()
    wv = (parv>=vmin) & (parv <=vmax)
    frac = sum(pixlr[~wv])/sum(pixlr)
    print('fraction of randoms not included in plot: '+str(frac))
    return bc,sv,ep 

def densvsimpar_pix(rl,ft,par,reg=None,wsel=None,xlab='',datweights=None,bl=None,fnc=None,vmin=None,vmax=None,ebvcut=None,edscut=None,sn2cut=None,fpsfcut=None,gfluxcut=None,rfluxcut=None,gbcut=None,nbin=10,weights=None,titl='',rng=0.3):        
    if bl is not None:
        wr = np.isin(rl['BRICKID'],bl)
        rl = rl[wr]
        wd = np.isin(ft['BRICKID'],bl)
        ft = ft[wd]
            #reg += ' Obiwan bricks'

    pixlr = gethpmap(rl,reg)
    print('randoms done')
    pixlg = gethpmap(ft,reg,datweights)
    print('data done')

    if weights is None:
        weights = np.ones(len(pixlr))

    if wsel is not None:
        wp = wsel
        wp &= (pixlr > 0) 
    else:
        wp = (pixlr > 0) 
    wp &= (weights*0 == 0)



    parv = fitsio.read(pixfn)
    ebv = parv['EBV']
    sn2tf = 10.**(-0.4*R_G*ebv*2.)*parv['PSFDEPTH_G'] + 10.**(-0.4*R_R*ebv*2.)*parv['PSFDEPTH_R'] + 10.**(-0.4*R_Z*ebv*2.)*parv['PSFDEPTH_Z']
    print(len(parv[wp]))
    if sn2cut:
        wp &= (sn2tf > sn2cut)
    
    if fpsfcut:
        wpsf = ft['MORPHTYPE'] == 'PSF'
        pixlgp = gethpmap(ft[wpsf],reg)
        fpsf = pixlgp/pixlg
        wp &= (fpsf < fpsfcut)
    if ebvcut:
        wp &= (parv['EBV'] < ebvcut)

    if edscut:
        eds = parv['EBV']/parv['STARDENS']
        wp &= (eds < edscut)

    if gbcut is not None:
    

        print('applying background cut of '+str(gbcut))
        rf = fitsio.read('/global/u2/r/rongpu/share/desi/sky_residual_dr9_partial/sky_residual_dr9_north_256.fits')
        gb = np.zeros(12*nside*nside)
        for i in range(0,len(rf)):
            px = rf['hp_idx'][i]
            gb[px] = rf['g_blobsky'][i]  
        gb = hp.reorder(gb,r2n=True)    
        wp &= (gb != 0)  
        wp &= (gb < gbcut)    

    
    print(len(parv[wp]))
    if type(par) == str:
        if par.split('-')[0] == 'VAR' or par.split('-')[0] == 'STDPER':
            pixlp,pixlv = gethpmap_var(ft,reg)
            if par.split('-')[0] == 'VAR':
                parv = pixlv[wp]/pixlg[wp]-(pixlp[wp]/pixlg[wp])**2.  
            elif par.split('-')[0] == 'STDPER':
                var = pixlv[wp]/pixlg[wp]-(pixlp[wp]/pixlg[wp])**2. 
                parv = var**.5/(pixlp[wp]/pixlg[wp])
        elif par == 'fracPSF':
            wpsf = ft['MORPHTYPE'] == 'PSF'
            pixlgp = np.zeros(12*nside*nside)
            dpixp = dpix[wpsf]
            for i in range(0,len(dpixp)): 
                pix = dpixp[i]
                pixlgp[pix] += 1.
            parv = pixlgp[wp]/pixlg[wp]
        else:
            parv = get_prop_map(par)[wp]
#         parsp = par.split('-')
#         if len(par.split('-')) > 1: 
#             if parsp[1] == 'EBV':
#                 ebv = parv[wp]['EBV']
#                 if '_R' in par:
#                     R_v = R_R
#                 if '_G' in par:
#                     R_v = R_g
#                 if '_Z' in par:
#                     R_v = R_Z
#                 parv = 10.**(-0.4*R_v*ebv*2.)*parv[wp][parsp[0]]
# 
#         
#             #print(par.split('-'))   
#             elif par.split('-')[1] == 'X':
#                 parv = parv[wp][par.split('-')[0]]*parv[wp][par.split('-')[2]]
#             elif par.split('-')[1] == 'DIV':
#                 parv = parv[wp][par.split('-')[0]]/parv[wp][par.split('-')[2]]
#         elif par == 'PSFTOT':
#             parv = (parv[wp]['PSFSIZE_G'])*(parv[wp]['PSFSIZE_R'])*(parv[wp]['PSFSIZE_Z'])
#         elif par == 'SN2TOT_FLAT':
#             ebv = parv[wp]['EBV']
#             parv = 10.**(-0.4*R_G*ebv*2.)*parv[wp]['PSFDEPTH_G'] + 10.**(-0.4*R_R*ebv*2.)*parv[wp]['PSFDEPTH_R'] + 10.**(-0.4*R_Z*ebv*2.)*parv[wp]['PSFDEPTH_Z']
#         elif par == 'SNTOT_FLAT':
#             ebv = parv[wp]['EBV']
#             parv = 10.**(-0.4*R_G*ebv*2.)*parv[wp]['PSFDEPTH_G'] + 10.**(-0.4*R_R*ebv*2.)*parv[wp]['PSFDEPTH_R'] + 10.**(-0.4*R_Z*ebv*2.)*parv[wp]['PSFDEPTH_Z']
#             parv = np.sqrt(parv)
#         elif par == 'SN2TOT_G':
#             ebv = parv[wp]['EBV']
#             parv = 10.**(-0.4*R_G*ebv*2.)*parv[wp]['PSFDEPTH_G']
# 
#         else:
#             parv = parv[wp][par]
    else:
        parv  = par[wp]
    
    bc,sv,ep = plot_pixdens1d(pixlg[wp],pixlr[wp],parv,weights[wp],vmin,vmax,titl=titl,xlab=xlab,rng=rng)
    return bc,sv,ep



class densvar:
    def __init__(self,type,sdir='',tv='0.49.0',rel='DR9',ti=None,elgandlrgbits = [1,5,6,7,8,9,11,12,13],columns=['RA','DEC','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z','MASKBITS','EBV','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','BRICKID']):
        df = sdir+type+'targets'+rel+'v'+tv+'.fits'
        columnsg = columns.copy()
        columnsr = columns.copy()
        
        if ti is not None:
            columnsg.append('MORPHTYPE')
        print(columnsr)
        ft = fitsio.read(df,columns=columnsg)
        print(len(ft))
        self.ft = mask(ft,mb=elgandlrgbits)
        if ti is not None:
            wt = self.ft['MORPHTYPE'] != ti[0]
            for i in range(1,len(ti)):
                wt &= self.ft['MORPHTYPE'] != ti[i]
            self.ft = self.ft[wt]    

        print(len(self.ft))
        del ft
        rl = fitsio.read(ranf,columns=columnsr)
        print(len(rl))
        self.rl = mask(rl,mb=elgandlrgbits)
        print(len(self.rl))
        del rl
        self.type = type




    def plot_brickdens(self,reg=False,sz=.2,vx=2):
        brickf = fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/survey-bricks.fits.gz') 
        brickdictrd = {}
        for i in range(0,len(brickf)):
            brickdictrd[brickf[i]['BRICKID']] = (brickf[i]['RA'],brickf[i]['DEC'])
        if reg:
            if reg == 'S' or reg == 'N':
                wr = self.rl['PHOTSYS'] == reg
                wd = self.ft['PHOTSYS'] == reg
            else:
                wr = sel_reg(self.rl['RA'],self.rl['DEC'],reg)
                wd = sel_reg(self.ft['RA'],self.ft['DEC'],reg)
            
            rl = self.rl[wr]        
            ft = self.ft[wd]
        else:
            rl = self.rl       
            ft = self.ft
        nbr = np.max(rl['BRICKID'])
        nbd = np.max(ft['BRICKID'])
        nbx = np.max([nbr,nbd])+1
        print('maximum brickid is '+str(nbx))
        pixlr = np.zeros(nbx)
        pixlg = np.zeros(nbx)
        for i in range(0,len(rl)):
            id = rl[i]['BRICKID']
            pixlr[id] += 1.
        print('randoms done')
        for i in range(0,len(ft)):
            id = ft[i]['BRICKID']
            pixlg[id] += 1.
        wp = pixlr > 0
        pixls = []
        for i in range(0,len(pixlr)):
            if pixlr[i] > 0:
                pixls.append(i)
        pixls = np.array(pixls).astype(int) 
        rap = []
        decp = []
        for id in pixls:
            rai,deci = brickdictrd[id]
            rap.append(rai)
            decp.append(deci)   
        od = pixlg[wp]/pixlr[wp]
        od = od/np.mean(od)
        decp = np.array(decp)
        if reg == 'DS':
            wr = rap > 250
            rap[wr] -=360
        if vx == None:
            vx = np.max(od)
        if vm == None:
            vm = np.min(od)    

        plt.scatter(rap,np.sin(decp*np.pi/180),c=od,s=sz,vmax=vx,vmin=vm)#,vmin=1.,vmax=2)
        plt.xlabel('RA')
        plt.ylabel('sin(DEC)')

        plt.show()

    def plot_brickprop(self,prop,reg=False,fnc=None,sz=.2,vx=None,vm=None,decmin=-90,decmax=91,decsp=30):
        brickf = fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/survey-bricks.fits.gz') 
        brickdictrd = {}
        for i in range(0,len(brickf)):
            brickdictrd[brickf[i]['BRICKID']] = (brickf[i]['RA'],brickf[i]['DEC'])
        if reg:
            if reg == 'S' or reg == 'N':
                wr = self.rl['PHOTSYS'] == reg
                wd = self.ft['PHOTSYS'] == reg
            else:
                wr = sel_reg(self.rl['RA'],self.rl['DEC'],reg)
                wd = sel_reg(self.ft['RA'],self.ft['DEC'],reg)
            
            rl = self.rl[wr]        
            ft = self.ft[wd]
        else:
            rl = self.rl       
            ft = self.ft
        nbd = np.max(ft['BRICKID'])
        nbx = nbd+1
        print('maximum brickid is '+str(nbx))
        pixlr = np.zeros(nbx)
        pixlg = np.zeros(nbx)

        for i in range(0,len(ft)):
            id = ft[i]['BRICKID']
            pixlr[id] += 1.
            pixlg[id] += ft[i][prop]
        wp = pixlr > 0
        pixls = []
        for i in range(0,len(pixlr)):
            if pixlr[i] > 0:
                pixls.append(i)
        pixls = np.array(pixls).astype(int) 
        rap = []
        decp = []
        for id in pixls:
            rai,deci = brickdictrd[id]
            rap.append(rai)
            decp.append(deci)   
        od = pixlg[wp]/pixlr[wp]
        if vx == None:
            vx = np.max(od)
        if vm == None:
            vm = np.min(od)    
        #od = od/np.mean(od)
        print(vm,vx)
        decp = np.array(decp)

        if reg == 'DS':
            wr = rap > 250
            rap[wr] -=360
    
    
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, .8, .8])
        ax.scatter(rap,np.sin(decp*np.pi/180),c=od,s=sz,vmax=vx,vmin=vm)
        ax.set_title(prop +' averaged in bricks')
        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')
        yt = np.arange(decmin,decmax,decsp)
        yl = []
        for i in range(0,len(yt)):
            yl.append(str(yt[i]))
        ax.set_yticks(np.sin(yt*np.pi/180.))
        ax.set_yticklabels(yl)
        rarn = np.max(rap)/90.-np.min(rap)/90.
        ar = (np.sin(np.pi/180.*decmax)-np.sin(np.pi/180.*decmin))/rarn
        print(ar,rarn,np.sin(np.pi/180.*decmax)-np.sin(np.pi/180.*decmin))
        ax.set_aspect(ar*180)
        fig.colorbar()

        plt.show()

    def plot_brickpropvar(self,prop,reg=False,fnc=None,sz=.2,vx=None,vm=None):
        brickf = fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/survey-bricks.fits.gz') 
        brickdictrd = {}
        for i in range(0,len(brickf)):
            brickdictrd[brickf[i]['BRICKID']] = (brickf[i]['RA'],brickf[i]['DEC'])
        if reg:
            if reg == 'S' or reg == 'N':
                wr = self.rl['PHOTSYS'] == reg
                wd = self.ft['PHOTSYS'] == reg
            else:
                wr = sel_reg(self.rl['RA'],self.rl['DEC'],reg)
                wd = sel_reg(self.ft['RA'],self.ft['DEC'],reg)
            
            rl = self.rl[wr]        
            ft = self.ft[wd]
        else:
            rl = self.rl       
            ft = self.ft
        nbd = np.max(ft['BRICKID'])
        nbx = nbd+1
        print('maximum brickid is '+str(nbx))
        pixlr = np.zeros(nbx)
        pixlg = np.zeros(nbx)
        pixlv = np.zeros(nbx)

        for i in range(0,len(ft)):
            id = ft[i]['BRICKID']
            pixlr[id] += 1.
            pixlg[id] += ft[i][prop]
            pixlv[id] += ft[i][prop]**2.
        wp = pixlr > 0
        pixls = []
        for i in range(0,len(pixlr)):
            if pixlr[i] > 0:
                pixls.append(i)
        pixls = np.array(pixls).astype(int) 
        rap = []
        decp = []
        for id in pixls:
            rai,deci = brickdictrd[id]
            rap.append(rai)
            decp.append(deci)   
        od = pixlv[wp]/pixlr[wp]-(pixlg[wp]/pixlr[wp])**2.
        if vx == None:
            vx = np.max(od)
        if vm == None:
            vm = np.min(od)    
        #od = od/np.mean(od)
        print(vm,vx)
        decp = np.array(decp)
        plt.scatter(rap,np.sin(decp*np.pi/180),c=od,s=sz,vmax=vx,vmin=vm)
        plt.title('variance in ' +prop +' per brick')
        plt.xlabel('RA')
        plt.ylabel('sin(DEC)')
        plt.colorbar()

        plt.show()

    def plot_brickprop_stdper(self,prop,reg=False,fnc=None,sz=.2,vx=None,vm=None,minn = 10):
        brickf = fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/survey-bricks.fits.gz') 
        brickdictrd = {}
        for i in range(0,len(brickf)):
            brickdictrd[brickf[i]['BRICKID']] = (brickf[i]['RA'],brickf[i]['DEC'])
        if reg:
            if reg == 'S' or reg == 'N':
                wr = self.rl['PHOTSYS'] == reg
                wd = self.ft['PHOTSYS'] == reg
            else:
                wr = sel_reg(self.rl['RA'],self.rl['DEC'],reg)
                wd = sel_reg(self.ft['RA'],self.ft['DEC'],reg)
            
            rl = self.rl[wr]        
            ft = self.ft[wd]
        else:
            rl = self.rl       
            ft = self.ft
        nbd = np.max(ft['BRICKID'])
        nbx = nbd+1
        print('maximum brickid is '+str(nbx))
        pixlr = np.zeros(nbx)
        pixlg = np.zeros(nbx)
        pixlv = np.zeros(nbx)

        for i in range(0,len(ft)):
            id = ft[i]['BRICKID']
            pixlr[id] += 1.
            pixlg[id] += ft[i][prop]
            pixlv[id] += ft[i][prop]**2.
    
        wp = pixlr > minn
        pixls = []
        for i in range(0,len(pixlr)):
            if pixlr[i] > minn:
                pixls.append(i)
        pixls = np.array(pixls).astype(int) 
        rap = []
        decp = []
        for id in pixls:
            rai,deci = brickdictrd[id]
            rap.append(rai)
            decp.append(deci)   
        od = (pixlv[wp]/pixlr[wp]-(pixlg[wp]/pixlr[wp])**2.)**.5/(pixlg[wp]/pixlr[wp])
        wo = od*0 == 0
        od = od[wo]
        if vx == None:
            vx = np.max(od)
        if vm == None:
            vm = np.min(od)    
        #od = od/np.mean(od)
        print(vm,vx)
        decp = np.array(decp)
        rap = np.array(rap)
        plt.scatter(rap[wo],np.sin(decp[wo]*np.pi/180),c=od,s=sz,vmax=vx,vmin=vm)
        plt.title('variance/mean in ' +prop +' per brick')
        plt.xlabel('RA')
        plt.ylabel('sin(DEC)')
        plt.colorbar()

        plt.show()



    def densvsimpar_ran(self,par,reg=None,fnc=None,vmin=None,vmax=None,nbin=10,bl=None):
        if reg:
            if reg == 'S' or reg == 'N':
                wr = self.rl['PHOTSYS'] == reg
                wd = self.ft['PHOTSYS'] == reg
            else:
                wr = sel_reg(self.rl['RA'],self.rl['DEC'],reg)
                wd = sel_reg(self.ft['RA'],self.ft['DEC'],reg)
            
            rl = self.rl[wr]        
            ft = self.ft[wd]
        else:
            rl = self.rl       
            ft = self.ft
        if vmin is None:
            vmin = np.min(rl[par])
        if vmax is None:
            vmax = np.max(rl[par])    
        
        if bl is not None:
            wr = np.isin(rl['BRICKID'],bl)
            rl = rl[wr]
            wd = np.isin(ft['BRICKID'],bl)
            ft = ft[wd]
            reg += ' Obiwan bricks'
        rh,bn = np.histogram(rl[par],bins=nbin,range=(vmin,vmax))
        dh,db = np.histogram(ft[par],bins=bn)
        rf = len(rl)/len(ft)
        sv = dh/rh*rf
        ep = np.sqrt(dh)/rh*rf
        bc = []
        for i in range(0,len(bn)-1):
            bc.append((bn[i]+bn[i+1])/2.)
        plt.errorbar(bc,sv-1.,ep,fmt='ko')
        plt.hist(rl[par],bins=nbin,range=(vmin,vmax),weights=0.2*np.ones(len(rl))/np.max(rh))
        plt.ylim(-.3,.3)
        plt.xlabel(par)
        plt.ylabel('Ngal/<Ngal> - 1')
        plt.title(self.type+' in '+reg + ' footprint')
        plt.show()
        wv = (rl[par]>vmin) & (rl[par] < vmax)
        frac = len(rl[~wv])/len(rl)
        print('fraction of randoms not included in plot: '+str(frac))
        return bc,sv,ep



    def densvsskyres_pix(self,par,reg=None,fnc=None,vmin=None,vmax=None,ebvcut=None,edscut=None,sn2cut=None,fpsfcut=None,gfluxcut=None,rfluxcut=None,gbcut=None,nbin=10,weights=None,titl=''):        
        #test against Rongpu's residuals
        print('SOMETHING SEEMS WRONG, DID YOU FIX IT??')
        if reg:
            if reg == 'S' or reg == 'N':
                wr = self.rl['PHOTSYS'] == reg
                wd = self.ft['PHOTSYS'] == reg
            else:
                wr = sel_reg(self.rl['RA'],self.rl['DEC'],reg)
                wd = sel_reg(self.ft['RA'],self.ft['DEC'],reg)
            
            rl = self.rl[wr]        
            ft = self.ft[wd]
        else:
            rl = self.rl       
            ft = self.ft
        rth,rphi = radec2thphi(rl['RA'],rl['DEC'])
        rpix = hp.ang2pix(nside,rth,rphi,nest=nest)
        dth,dphi = radec2thphi(ft['RA'],ft['DEC'])
        dpix = hp.ang2pix(nside,dth,dphi,nest=nest)
        pixlr = np.zeros(12*nside*nside)
        pixlg = np.zeros(12*nside*nside)

        if weights is None:
            weights = np.ones(len(pixlr))
        for pix in rpix:
            pixlr[pix] += 1.
        print('randoms done')
        for i in range(0,len(dpix)): 
            pix = dpix[i]
            pixlg[pix] += 1.
    
        wp = (pixlr > 0) & (weights*0 == 0)

        parv = fitsio.read(pixfn)
        ebv = parv['EBV']
        sn2tf = 10.**(-0.4*R_G*ebv*2.)*parv['PSFDEPTH_G'] + 10.**(-0.4*R_R*ebv*2.)*parv['PSFDEPTH_R'] + 10.**(-0.4*R_Z*ebv*2.)*parv['PSFDEPTH_Z']
        print(len(parv[wp]))
        if sn2cut:
            wp &= (sn2tf > sn2cut)
        
        if fpsfcut:
            wpsf = ft['MORPHTYPE'] == 'PSF'
            pixlgp = np.zeros(12*nside*nside)
            dpixp = dpix[wpsf]
            for i in range(0,len(dpixp)): 
                pix = dpixp[i]
                pixlgp[pix] += 1.
            fpsf = pixlgp/pixlg
            wp &= (fpsf < fpsfcut)
        if ebvcut:
            wp &= (parv['EBV'] < ebvcut)

        if edscut:
            eds = parv['EBV']/parv['STARDENS']
            wp &= (eds < edscut)
    
        rf = fitsio.read('/global/u2/r/rongpu/share/desi/sky_residual_dr9_partial/sky_residual_dr9_north_256.fits')
        parv = np.zeros(12*nside*nside)
        for i in range(0,len(rf)):
            px = rf['hp_idx'][i]
            parv[px] = rf[par][i]  
        #parv = hp.reorder(parv,r2n=True)      

        wp &= parv !=0
        wp &= parv*0 == 0
        print(len(parv[wp]))
    
        if vmin is None:
            vmin = np.min(parv[wp])
        if vmax is None:
            vmax = np.max(parv[wp])
        parv = parv[wp]
        rh,bn = np.histogram(parv,bins=nbin,range=(vmin,vmax),weights=pixlr[wp])
        dh,db = np.histogram(parv,bins=bn,weights=pixlg[wp]*weights[wp])
        norm = sum(rh)/sum(dh)
        sv = dh/rh*norm
        ep = np.sqrt(dh)/rh*norm
        bc = []
        for i in range(0,len(bn)-1):
            bc.append((bn[i]+bn[i+1])/2.)

        plt.errorbar(bc,sv-1.,ep,fmt='ko')
        plt.hist(parv,bins=nbin,range=(vmin,vmax),weights=pixlr[wp]*0.2*np.ones(len(pixlr[wp]))/np.max(rh))
        plt.ylim(-.3,.3)
        plt.xlabel(par)
        plt.ylabel('Ngal/<Ngal> - 1')
        plt.title(self.type+' in '+reg + ' footprint, using pixelized map'+titl)
        plt.show()
        wv = (parv>=vmin) & (parv <=vmax)
        frac = sum(pixlr[wp][~wv])/sum(pixlr[wp])
        print('fraction of randoms not included in plot: '+str(frac))
        return bc,sv,ep

   
    #plot density vs depth with healpix values
    def plotvshp_compmc(self,sys,rng,mcl=None,ws=None,reg=None,fnc=None,gdzm=0,ebvm=100,title='',effac=1.,mingd=0,maxgd=1.e6,minpsfg=0,maxpsfg=100,south=True):
        if reg:
            if reg == 'S' or reg == 'N':
                wr = self.rl['PHOTSYS'] == reg
                wd = self.ft['PHOTSYS'] == reg
            else:
                wr = sel_reg(self.rl['RA'],self.rl['DEC'],reg)
                wd = sel_reg(self.ft['RA'],self.ft['DEC'],reg)
            
            rl = self.rl[wr]        
            ft = self.ft[wd]
        else:
            rl = self.rl       
            ft = self.ft
        rth,rphi = radec2thphi(rl['RA'],rl['DEC'])
        rpix = hp.ang2pix(nside,rth,rphi,nest=nest)
        dth,dphi = radec2thphi(ft['RA'],ft['DEC'])
        dpix = hp.ang2pix(nside,dth,dphi,nest=nest)
        r1 = np.zeros(12*nside*nside)
        d1= np.zeros(12*nside*nside)
        for pix in rpix:
            r1[pix] += 1.
        print('randoms done')
        for pix in dpix:
            d1[pix] += 1.

        hpq = fitsio.read(pixfn)
        #hpq = parv[par]

    
        w = r1 > 0
        print(len(hpq[w]))
        w &= hpq['GALDEPTH_Z'] > gdzm
        w &= hpq['GALDEPTH_G'] > mingd
        w &= hpq['GALDEPTH_G'] < maxgd
        w &= hpq['EBV'] < ebvm
        w &= hpq['PSFSIZE_G'] > minpsfg
        w &= hpq['PSFSIZE_G'] < maxpsfg
        if ws is not None:
            w &= ws*0 == 0
        if mcl is not None:
            w &= mcl*0 == 0
            w &= mcl > 0
        print(len(hpq[w]))
        #w 
        if sys != 'gdc' and sys != 'rdc' and sys != 'zdc' and sys != 'dg' and sys != 'dr' and sys != 'dz' and sys != 'dgr' and sys != 'drz' and sys != 'dgz':
            sm = hpq[w][sys]
            xlab = sys
        else:
            if sys == 'gdc':
                print('g band depth, extinction corrected')
                sm = hpq[w]['GALDEPTH_G']*10.**(-0.4*R_G*hpq[w]['EBV'])
                xlab = 'g band depth, extinction corrected'
            if sys == 'rdc':
                sm = hpq[w]['GALDEPTH_R']*10.**(-0.4*R_R*hpq[w]['EBV'])
                xlab = 'r band depth, extinction corrected'
            if sys == 'zdc':
                sm = hpq[w]['GALDEPTH_Z']*10.**(-0.4*R_Z*hpq[w]['EBV'])
                xlab = 'z band depth, extinction corrected'
            if sys == 'dg':
                sm = dg[w]
                xlab = 'g band PS1 residual'
            if sys == 'dr':
                sm = dr[w]
                xlab = 'r band PS1 residual'
            if sys == 'dz':
                sm = dz[w]
                xlab = 'z band PS1 residual'
            if sys == 'dgr':
                sm = dg[w]-dr[w]
                xlab = 'g-r band PS1 residual'
            if sys == 'drz':
                sm = dr[w]-dz[w]
                xlab = 'r-z band PS1 residual'
            if sys == 'dgz':
                sm = dg[w]-dz[w]
                xlab = 'g-z band PS1 residual'

        ds = np.ones(len(d1))
        print(len(ds),len(d1),len(w),len(sm))
        hdnoc = np.histogram(sm,weights=d1[w],range=rng)
        #print(hd1)
        hr1 = np.histogram(sm,weights=r1[w],bins=hdnoc[1],range=rng)



        xl = []
        for i in range(0,len(hr1[0])):
            xl.append((hr1[1][i]+hr1[1][i+1])/2.)

        plt.errorbar(xl,hdnoc[0]/hr1[0]/(sum(d1[w])/sum(r1[w])),np.sqrt(hdnoc[0])/hr1[0]/(sum(d1[w])/sum(r1[w])),fmt='ko',label='raw')
        if ws is not None:
            ds = ws
            hd1 = np.histogram(sm,weights=d1[w]*ds[w],bins=hdnoc[1],range=rng)
            plt.plot(xl,hd1[0]/hr1[0]/(sum(d1[w]*ds[w])/sum(r1[w])),'b-',label='+ EBV weights')

        #hd1 = np.histogram(sm,weights=d1[w]*ds[w],bins=hdnoc[1],range=rng)
        #plt.plot(xl,hd1[0]/hr1[0]/(sum(d1[w]*ds[w])/sum(r1[w])),'k--',label='with stellar density weights')
        if mcl is not None:
            dmcse = mcl**effac
            hd1 = np.histogram(sm,weights=d1[w]/dmcse[w],bins=hdnoc[1],range=rng)
            plt.plot(xl,hd1[0]/hr1[0]/(sum(d1[w]/dmcse[w])/sum(r1[w])),'r-',label='+MC weights')
        if ws is not None and mcl is not None:
            hd1 = np.histogram(sm,weights=d1[w]*ds[w]/dmcse[w],bins=hdnoc[1],range=rng)
            plt.plot(xl,hd1[0]/hr1[0]/(sum(d1[w]*ds[w]/dmcse[w])/sum(r1[w])),'-',color='purple',label='+MC weights + EBV weights')
        #dmcs = mcls**effac
        #hd1 = np.histogram(sm,weights=d1[w]*ds[w]/dmcs[w],bins=hdnoc[1],range=rng)
        #plt.plot(xl,hd1[0]/hr1[0]/(sum(d1[w]*ds[w]/dmcs[w])/sum(r1[w])),'b-',label='+MC; sed w ext sigma')
        #dmco = mclo**effac
        #hd1 = np.histogram(sm,weights=d1[w]*ds[w]/dmco[w],bins=hdnoc[1],range=rng)
        #plt.plot(xl,hd1[0]/hr1[0]/(sum(d1[w]*ds[w]/dmco[w])/sum(r1[w])),'-',color='purple',label='old MC')
    
        #plt.title(str(mp)+reg)
        plt.plot(xl,np.ones(len(xl)),'k:',label='null')
        plt.legend()#(['raw','with stellar density weights','+sed ext MC','just sed MC','old MC','null']))
        plt.ylabel('relative density')
        plt.xlabel(xlab)
        plt.ylim(0.7,1.3)
        plt.title(title)
        plt.show()    

class cell:
    def __init__(self,dat,ran,randir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/',ranallf='randoms-allsky-1-0.fits'):
        self.dat = dat
        self.ran = ran
        ranall = fitsio.read(randir+ranallf,columns=['RA','DEC'])
        th,phi = radec2thphi(ranall['RA'],ranall['DEC'])
        ranpix = hp.ang2pix(256,th,phi)
        ranpall = np.zeros(12*256*256)
        for pix in ranpix:
            ranpall[pix] += 1.
        self.ranpall = ranpall
 
    def get_delta(self,reg,racol='RA',decol='DEC',wts=None,wtspix=None,thresh=0,nest=False,appfrac=True):#,ranpall=None
        dat = sel_reg(dat[racol],dat[decol],reg)
        ran = sel_reg(ran[racol],ran[decol],reg)
        th,phi = radec2thphi(dat[racol],dat[decol])
        datpix = hp.ang2pix(256,th,phi,nest=nest)
        datp = np.zeros(12*256*256)
        for i in range(0,len(datpix)):
            pix = datpix[i]
            if wts is not None:
                datp[pix] += wts[i]
            else:
                datp[pix] += 1.
        if wtspix is not None:
            datp *= wtspix
        th,phi = radec2thphi(ran[racol],ran[decol])
        ranpix = hp.ang2pix(256,th,phi,nest=nest)
        ranp = np.zeros(12*256*256)
        for pix in ranpix:
            ranp[pix] += 1.
        #ranp /= rannorm
    
        sel = ranp > thresh
        mnr = np.mean(datp[sel]/ranp[sel])
        print(mnr)
        delta = (datp/ranp/mnr -1)
        #if ranpall is not None:
        if appfrac:
            if nest:
                frac = ranp/self.ranpall
            #else:
            #    frac = ranp/ranpall_nest
            delta *= frac
        delta[~sel] = hp.UNSEEN
        fsky = np.sum(ranp[sel])/np.sum(ranpall)
        return delta,fsky 

            
if __name__ == "__main__":
    obiLRGvs_depthmag('S','psfdepth','g',md='sv3',syspix=True)#,vmin=24.3,vmax=25)
    obiLRGvs_depthmag('S','psfdepth','r',md='sv3',syspix=True)#,vmin=23.7,vmax=24.5)
    obiLRGvs_depthmag('S','psfdepth','z',md='sv3',syspix=True)#,vmin=22.7,vmax=23.7)
    
    