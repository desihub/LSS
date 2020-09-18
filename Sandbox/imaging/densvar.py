import fitsio
import astropy.io.fits as fits
import healpy as hp
import numpy as np
from matplotlib import pyplot as plt

pixfn      = '/global/cfs/cdirs/desi/target/catalogs/dr9m/0.42.0/pixweight/main/resolve/dark/pixweight-dark.fits'
hdr        = fits.getheader(pixfn,1)
nside,nest = hdr['HPXNSIDE'],hdr['HPXNEST']
print(nside,nest)


dr = '9'

if dr == '9':
    #this will be needed no matter the sample, might want more
    rall = fitsio.read('/global/cfs/cdirs/desi/target/catalogs/dr9m/0.42.0/randoms/resolve/randoms-randomized-1.fits')
    print(len(rall))
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

rall = mask(rall)    

def radec2thphi(ra,dec):
    return (-dec+90.)*np.pi/180.,ra*np.pi/180.
    
def thphi2radec(theta,phi):
    return 180./np.pi*phi,-(180./np.pi*theta-90)


def plot_hpdens(type,reg=False,ff='targetDR9m42.fits',sz=.2,vx=2,weights=None):
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
        if pixlr[i] > 0:
            pixls.append(i)
    pixls = np.array(pixls).astype(int)        
    th,phi = hp.pix2ang(nside,pixls,nest=nest)
    od = pixlg[wp]/pixlr[wp]*weights[wp]
    od = od/np.mean(od)
    ra,dec = thphi2radec(th,phi)
    plt.scatter(ra,np.sin(dec*np.pi/180),c=od,s=sz,vmax=vx)#,vmin=1.,vmax=2)
    plt.show()

def plot_brickdens(type,reg=False,ff='targetDR9m42.fits',sz=.2,vx=2):
    brickf = fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/survey-bricks.fits.gz') 
    brickdictrd = {}
    for i in range(0,len(brickf)):
        brickdictrd[brickf[i]['BRICKID']] = (brickf[i]['RA'],brickf[i]['DEC'])
    ft = fitsio.read(sdir+type+ff)
    print(len(ft))
    rl = rall
    if reg:
        wr = rall['PHOTSYS'] == reg
        rl = rl[wr]
        wd = ft['PHOTSYS'] == reg
        ft = ft[wd]
    nbr = np.max(rl['BRICKID'])
    nbd = np.max(ft['BRICKID'])
    nbx = np.max([nbr,nbd])
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
    plt.show()


def densvsimpar_ran(type,par,reg=None,ff='targetDR9m42.fits',vmin=None,vmax=None,nbin=10):
    ft = fitsio.read(sdir+type+ff)
    print(len(ft))
    rl = rall
    if reg == None:
        reg = 'All'
    else:    
        wr = rall['PHOTSYS'] == reg
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

def densvsimpar_pix(type,par,reg=None,ff='targetDR9m42.fits',vmin=None,vmax=None,nbin=10):        
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
    pixlr = np.zeros(12*nside*nside)
    pixlg = np.zeros(12*nside*nside)
    for pix in rpix:
        pixlr[pix] += 1.
    print('randoms done')
    for pix in dpix:
        pixlg[pix] += 1.
    parv = fitsio.read(pixfn)
    parv = parv[par]
    wp = pixlr > 0
    if vmin is None:
    	vmin = np.min(parv[wp])
    if vmax is None:
        vmax = np.max(parv[wp])
    rh,bn = np.histogram(parv[wp],bins=nbin,range=(vmin,vmax),weights=pixlr[wp])
    dh,db = np.histogram(parv[wp],bins=bn,weights=pixlg[wp])
    norm = sum(rh)/sum(dh)
    sv = dh/rh*norm
    ep = np.sqrt(dh)/rh*norm
    bc = []
    for i in range(0,len(bn)-1):
        bc.append((bn[i]+bn[i+1])/2.)

    plt.errorbar(bc,sv-1.,ep,fmt='ko')
    plt.hist(parv[wp],bins=nbin,range=(vmin,vmax),weights=pixlr[wp]*0.2*np.ones(len(pixlr[wp]))/np.max(rh))
    plt.ylim(-.3,.3)
    plt.xlabel(par)
    plt.ylabel('Ngal/<Ngal> - 1')
    plt.title(type+' in '+reg + ' footprint, using pixelized map')
    plt.show()
    wv = (parv[wp]>=vmin) & (parv[wp] <=vmax)
    frac = sum(pixlr[wp][~wv])/sum(pixlr[wp])
    print('fraction of randoms not included in plot: '+str(frac))
   
            

    
    