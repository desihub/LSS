import fitsio
import astropy.io.fits as fits
import healpy as hp
import numpy as np
from matplotlib import pyplot as plt

pixfn      = '/global/cfs/cdirs/desi/target/catalogs/dr9m/0.44.0/pixweight/main/resolve/dark/pixweight-dark.fits'
hdr        = fits.getheader(pixfn,1)
nside,nest = hdr['HPXNSIDE'],hdr['HPXNEST']
print(nside,nest)

R_G=3.214 # http://legacysurvey.org/dr8/catalogs/#galactic-extinction-coefficients
R_R=2.165
R_Z=1.211

dr = '9'

fidf = 'targetDR9m44.fits'
ranf = '/global/cfs/cdirs/desi/target/catalogs/dr9m/0.44.0/randoms/resolve/randoms-1-0.fits'

if dr == '9':
    #this will be needed no matter the sample, might want more
    #rall = fitsio.read('/global/cfs/cdirs/desi/target/catalogs/dr9m/0.44.0/randoms/resolve/randoms-1-0.fits')
    #print(len(rall))
    sdir = '/project/projectdirs/desi/users/ajross/dr9/'

def mask(dd,mb=[1,5,6,7,11,12,13]):
    keep = (dd['NOBS_G']>0) & (dd['NOBS_R']>0) & (dd['NOBS_Z']>0)
    print(len(dd[keep]))
    
    keepelg = keep
    for bit in mb:
        keepelg &= ((dd['MASKBITS'] & 2**bit)==0)
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

#def split_reg():

#rall = mask(rall)    

def radec2thphi(ra,dec):
    return (-dec+90.)*np.pi/180.,ra*np.pi/180.
    
def thphi2radec(theta,phi):
    return 180./np.pi*phi,-(180./np.pi*theta-90)


def plot_hpdens(type,reg=False,fnc=None,sz=.2,vx=1.5,vm=.5,weights=None):
    if fnc is None:
        ff = fidf
    else:
        ff = fnc    
    ft = fitsio.read(sdir+type+ff,columns=['RA','DEC','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z','MASKBITS'])
    
    print(len(ft))
    ft = mask(ft)
    print(len(ft))
    rl = fitsio.read(ranf,columns=['RA','DEC','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z','MASKBITS'])
    print(len(rl))
    rl = mask(rl)
    print(len(rl))
    if reg:
        if reg == 'S' or reg == 'N':
            wr = rl['PHOTSYS'] == reg
            wd = ft['PHOTSYS'] == reg
        else:
            wr = sel_reg(rl['RA'],rl['DEC'],reg)
            wd = sel_reg(ft['RA'],ft['DEC'],reg)
            
        rl = rl[wr]        
        ft = ft[wd]
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
    for pix in dpix:
        pixlg[pix] += 1.
    wp = (pixlr > 0) & (weights*0 == 0)
    pixls = []
    for i in range(0,len(pixlr)):
        if pixlr[i] > 0 and weights[i]*0 == 0:
            pixls.append(i)
    pixls = np.array(pixls).astype(int)        
    th,phi = hp.pix2ang(nside,pixls,nest=nest)
    od = pixlg[wp]/pixlr[wp]*weights[wp]
    od = od/np.mean(od)
    ra,dec = thphi2radec(th,phi)
    if reg == 'DS':
        wr = ra > 250
        ra[wr] -=360
    if vx == None:
        vx = np.max(od)
    if vm == None:
        vm = np.min(od)    

    plt.scatter(ra,np.sin(dec*np.pi/180),c=od,s=sz,vmax=vx,vmin=vm)#,vmin=1.,vmax=2)
    plt.xlabel('RA')
    plt.ylabel('sin(DEC)')
    plt.colorbar()
    plt.title('relative '+type+' density')
    plt.show()

def plot_hpprop(par,type='ELG',reg=False,fnc=None,sz=.2,vx=None,vm=None,weights=None):
    if fnc is None:
        ff = fidf
    else:
        ff = fnc    
    ft = fitsio.read(sdir+type+ff,columns=['RA','DEC','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z','MASKBITS'])
    
    print(len(ft))
    ft = mask(ft)
    print(len(ft))
    rl = fitsio.read(ranf,columns=['RA','DEC','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z','MASKBITS'])
    print(len(rl))
    rl = mask(rl)
    print(len(rl))
    if reg:
        if reg == 'S' or reg == 'N':
            wr = rl['PHOTSYS'] == reg
            wd = ft['PHOTSYS'] == reg
        else:
            wr = sel_reg(rl['RA'],rl['DEC'],reg)
            wd = sel_reg(ft['RA'],ft['DEC'],reg)
            
        rl = rl[wr]        
        ft = ft[wd]
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
    for pix in dpix:
        pixlg[pix] += 1.
    wp = (pixlr > 0) & (weights*0 == 0)
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
    pixls = []
    for i in range(0,len(pixlr)):
        if pixlr[i] > 0 and weights[i]*0 == 0:
            pixls.append(i)
    pixls = np.array(pixls).astype(int)        
    th,phi = hp.pix2ang(nside,pixls,nest=nest)
    od = parv
    if vx == None:
        vx = np.max(od)
    if vm == None:
        vm = np.min(od)    
    
    ra,dec = thphi2radec(th,phi)
    if reg == 'DS':
        wr = ra > 250
        ra[wr] -= 360
    plt.scatter(ra,np.sin(dec*np.pi/180),c=od,s=sz,vmax=vx,vmin=vm)#,vmin=1.,vmax=2)
    plt.xlabel('RA')
    plt.ylabel('sin(DEC)')
    plt.colorbar()
    plt.title(par)

    plt.show()


def plot_brickdens(type,reg=False,sz=.2,vx=2):
    brickf = fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/survey-bricks.fits.gz') 
    brickdictrd = {}
    for i in range(0,len(brickf)):
        brickdictrd[brickf[i]['BRICKID']] = (brickf[i]['RA'],brickf[i]['DEC'])
    if fnc is None:
        ff = fidf
    else:
        ff = fnc    
    ft = fitsio.read(sdir+type+ff)
    print(len(ft))
    rl = fitsio.read(ranf,columns=['RA','DEC','PHOTSYS'])
    if reg:
        wr = rl['PHOTSYS'] == reg
        rl = rl[wr]
        wd = ft['PHOTSYS'] == reg
        ft = ft[wd]
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
    plt.scatter(rap,np.sin(decp*np.pi/180),c=od,s=sz,vmax=vx)#,vmin=1.,vmax=2)
    plt.xlabel('RA')
    plt.ylabel('sin(DEC)')

    plt.show()

def plot_brickprop(type,prop,reg=False,fnc=None,sz=.2,vx=None,vm=None,decmin=-90,decmax=91,decsp=30):
    brickf = fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/survey-bricks.fits.gz') 
    brickdictrd = {}
    for i in range(0,len(brickf)):
        brickdictrd[brickf[i]['BRICKID']] = (brickf[i]['RA'],brickf[i]['DEC'])
    if fnc is None:
        ff = fidf
    else:
        ff = fnc    
    ft = fitsio.read(sdir+type+ff,columns=['RA','DEC','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z','MASKBITS','BRICKID',prop])
    print(len(ft))
    if reg:
        wd = ft['PHOTSYS'] == reg
        ft = ft[wd]
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

def plot_brickpropvar(type,prop,reg=False,fnc=None,sz=.2,vx=None,vm=None):
    brickf = fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/survey-bricks.fits.gz') 
    brickdictrd = {}
    for i in range(0,len(brickf)):
        brickdictrd[brickf[i]['BRICKID']] = (brickf[i]['RA'],brickf[i]['DEC'])
    if fnc is None:
        ff = fidf
    else:
        ff = fnc    
    ft = fitsio.read(sdir+type+ff,columns=['RA','DEC','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z','MASKBITS','BRICKID',prop])
    print(len(ft))
    if reg:
        wd = ft['PHOTSYS'] == reg
        ft = ft[wd]
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

def plot_brickprop_stdper(type,prop,reg=False,fnc=None,sz=.2,vx=None,vm=None,minn = 10):
    brickf = fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/survey-bricks.fits.gz') 
    brickdictrd = {}
    for i in range(0,len(brickf)):
        brickdictrd[brickf[i]['BRICKID']] = (brickf[i]['RA'],brickf[i]['DEC'])
    if fnc is None:
        ff = fidf
    else:
        ff = fnc    
    ft = fitsio.read(sdir+type+ff,columns=['RA','DEC','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z','MASKBITS','BRICKID',prop])
    print(len(ft))
    if reg:
        wd = ft['PHOTSYS'] == reg
        ft = ft[wd]
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



def densvsimpar_ran(type,par,reg=None,fnc=None,vmin=None,vmax=None,nbin=10):
    if fnc is None:
        ff = fidf
    else:
        ff = fnc    
    print('using '+sdir+type+ff)    

    ft = fitsio.read(sdir+type+ff)
    ft = mask(ft)
    print(len(ft))
    rl = fitsio.read(ranf,columns=['RA','DEC','PHOTSYS',par,'NOBS_G','NOBS_R','NOBS_Z','MASKBITS'])
    rl = mask(rl)
    print(len(rl))
    if reg == None:
        reg = 'All'
    else:    
        wr = rl['PHOTSYS'] == reg
        rl = rl[wr]
        wd = ft['PHOTSYS'] == reg
        ft = ft[wd]
    if vmin is None:
        vmin = np.min(rl[par])
    if vmax is None:
        vmax = np.max(rl[par])    
        
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
    plt.title(type+' in '+reg + ' footprint')
    plt.show()
    wv = (rl[par]>vmin) & (rl[par] < vmax)
    frac = len(rl[~wv])/len(rl)
    print('fraction of randoms not included in plot: '+str(frac))

def densvsinput_pix(type,parl,wsel,reg=None,fnc=None,xlab='',vmin=None,vmax=None,ebvcut=None,edscut=None,sn2cut=None,fpsfcut=None,gfluxcut=None,rfluxcut=None,gbcut=None,nbin=10,weights=None,titl=''):        
    #input custom map/mask
    if fnc is None:
        ff = fidf
    else:
        ff = fnc    

    ft = fitsio.read(sdir+type+ff)
    print(len(ft))
    rl = fitsio.read(ranf,columns=['RA','DEC','PHOTSYS'])
    if reg:
        wr = rl['PHOTSYS'] == reg
        rl = rl[wr]
        wd = ft['PHOTSYS'] == reg
        ft = ft[wd]
    if gfluxcut:
        wg = ft['FLUX_G']/ft['MW_TRANSMISSION_G'] > gfluxcut
        print(len(ft))      
        ft = ft[wg]
        print(len(ft))
    if rfluxcut:
        wg = ft['FLUX_R']/ft['MW_TRANSMISSION_R'] > rfluxcut
        ft = ft[wg]
        

    
            
        
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
    
    wp = wsel
    wp &= (pixlr > 0) & (weights*0 == 0)

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
    plt.xlabel(xlab)
    plt.ylabel('Ngal/<Ngal> - 1')
    plt.title(type+' in '+reg + ' footprint, using pixelized map'+titl)
    plt.show()
    wv = (parv>=vmin) & (parv <=vmax)
    frac = sum(pixlr[wp][~wv])/sum(pixlr[wp])
    print('fraction of randoms not included in plot: '+str(frac))
    return bc,sv,ep


def densvsskyres_pix(type,par,reg=None,fnc=None,vmin=None,vmax=None,ebvcut=None,edscut=None,sn2cut=None,fpsfcut=None,gfluxcut=None,rfluxcut=None,gbcut=None,nbin=10,weights=None,titl=''):        
    #test against Rongpu's residuals
    if fnc is None:
        ff = fidf
    else:
        ff = fnc    

    ft = fitsio.read(sdir+type+ff)
    print(len(ft))
    rl = rall
    if reg:
        wr = rall['PHOTSYS'] == reg
        rl = rl[wr]
        wd = ft['PHOTSYS'] == reg
        ft = ft[wd]
    if gfluxcut:
        wg = ft['FLUX_G']/ft['MW_TRANSMISSION_G'] > gfluxcut
        print(len(ft))      
        ft = ft[wg]
        print(len(ft))
    if rfluxcut:
        wg = ft['FLUX_R']/ft['MW_TRANSMISSION_R'] > rfluxcut
        ft = ft[wg]
        

    
            
        
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
    parv = hp.reorder(parv,r2n=True)      

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
    plt.title(type+' in '+reg + ' footprint, using pixelized map'+titl)
    plt.show()
    wv = (parv>=vmin) & (parv <=vmax)
    frac = sum(pixlr[wp][~wv])/sum(pixlr[wp])
    print('fraction of randoms not included in plot: '+str(frac))
    return bc,sv,ep

def densvsimpar_pix(type,par,reg=None,fnc=None,vmin=None,vmax=None,ebvcut=None,edscut=None,sn2cut=None,fpsfcut=None,gfluxcut=None,rfluxcut=None,gbcut=None,nbin=10,weights=None,titl=''):        
    if fnc is None:
        ff = fidf
    else:
        ff = fnc    
    ft = fitsio.read(sdir+type+ff)
    print(len(ft))
    ft = mask(ft)
    print(len(ft))
    rl = fitsio.read(ranf,columns=['RA','DEC','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z','MASKBITS'])
    print(len(rl))
    rl = mask(rl)
    print(len(rl))
    if reg:
        if reg == 'S' or reg == 'N':
            wr = rl['PHOTSYS'] == reg
            wd = ft['PHOTSYS'] == reg
        else:
            wr = sel_reg(rl['RA'],rl['DEC'],reg)
            wd = sel_reg(ft['RA'],ft['DEC'],reg)
            
        rl = rl[wr]        
        ft = ft[wd]
    if gfluxcut:
        wg = ft['FLUX_G']/ft['MW_TRANSMISSION_G'] > gfluxcut
        print(len(ft))      
        ft = ft[wg]
        print(len(ft))
    if rfluxcut:
        wg = ft['FLUX_R']/ft['MW_TRANSMISSION_R'] > rfluxcut
        ft = ft[wg]

    
            
        
    rth,rphi = radec2thphi(rl['RA'],rl['DEC'])
    rpix = hp.ang2pix(nside,rth,rphi,nest=nest)
    dth,dphi = radec2thphi(ft['RA'],ft['DEC'])
    dpix = hp.ang2pix(nside,dth,dphi,nest=nest)
    pixlr = np.zeros(12*nside*nside)
    pixlg = np.zeros(12*nside*nside)
    if par.split('-')[0] == 'VAR' or par.split('-')[0] == 'STDPER':
        pixlp = np.zeros(12*nside*nside)
        pixlv = np.zeros(12*nside*nside)
    if weights is None:
        weights = np.ones(len(pixlr))
    for pix in rpix:
        pixlr[pix] += 1.
    print('randoms done')
    for i in range(0,len(dpix)): 
        pix = dpix[i]
        pixlg[pix] += 1.
        if par.split('-')[0] == 'VAR' or par.split('-')[0] == 'STDPER':
            pixlp[pix] += ft[i][par.split('-')[1]]
            pixlv[pix] += ft[i][par.split('-')[1]]**2.
    
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
    if len(par.split('-')) > 1: 
        
        if par.split('-')[0] == 'VAR':
            parv = pixlv[wp]/pixlg[wp]-(pixlp[wp]/pixlg[wp])**2.  
        elif par.split('-')[0] == 'STDPER':
            var = pixlv[wp]/pixlg[wp]-(pixlp[wp]/pixlg[wp])**2. 
            parv = var**.5/(pixlp[wp]/pixlg[wp])
        elif par.split('-')[1] == 'X':
            parv = parv[wp][par.split('-')[0]]*parv[wp][par.split('-')[2]]
        elif par.split('-')[1] == 'DIV':
            parv = parv[wp][par.split('-')[0]]/parv[wp][par.split('-')[2]]
    elif par == 'PSFTOT':
        parv = (parv[wp]['PSFSIZE_G'])*(parv[wp]['PSFSIZE_R'])*(parv[wp]['PSFSIZE_Z'])
    elif par == 'SN2TOT_FLAT':
        ebv = parv[wp]['EBV']
        parv = 10.**(-0.4*R_G*ebv*2.)*parv[wp]['PSFDEPTH_G'] + 10.**(-0.4*R_R*ebv*2.)*parv[wp]['PSFDEPTH_R'] + 10.**(-0.4*R_Z*ebv*2.)*parv[wp]['PSFDEPTH_Z']

    elif par == 'SN2TOT_G':
        ebv = parv[wp]['EBV']
        parv = 10.**(-0.4*R_G*ebv*2.)*parv[wp]['PSFDEPTH_G']

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

    wo = parv*0 == 0
    if vmin is None:
        vmin = np.min(parv[wo])
    if vmax is None:
        vmax = np.max(parv[wo])
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
    plt.title(type+' in '+reg + ' footprint, using pixelized map'+titl)
    plt.show()
    wv = (parv>=vmin) & (parv <=vmax)
    frac = sum(pixlr[wp][~wv])/sum(pixlr[wp])
    print('fraction of randoms not included in plot: '+str(frac))
    return bc,sv,ep
   
#plot density vs depth with healpix values
def plotvshp_compmc(type,sys,rng,mcl=None,ws=None,reg=None,fnc=None,gdzm=0,ebvm=100,title='',effac=1.,mingd=0,maxgd=1.e6,minpsfg=0,maxpsfg=100,south=True):
    if fnc is None:
        ff = fidf
    else:
        ff = fnc    
    ft = fitsio.read(sdir+type+ff)
    print(len(ft))
    rl = rall
    if reg:
        wr = rall['PHOTSYS'] == reg
        rl = rl[wr]
        wd = ft['PHOTSYS'] == reg
        ft = ft[wd]
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
            

    
    