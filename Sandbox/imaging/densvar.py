import fitsio
import astropy.io.fits as fits
import healpy as hp
import numpy as np
from matplotlib import pyplot as plt

pixfn      = '/project/projectdirs/desi/target/catalogs/dr8/0.31.1/pixweight/pixweight-dr8-0.31.1.fits'
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


def plot_hpdens(type,reg=False,ff='targetDR9m42.fits',sz=.2,vx=2):
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
    wp = pixlr > 0
    pixls = []
    for i in range(0,len(pixlr)):
        if pixlr[i] > 0:
            pixls.append(i)
    pixls = np.array(pixls).astype(int)        
    th,phi = hp.pix2ang(nside,pixls,nest=nest)
    od = pixlg[wp]/pixlr[wp]
    od = od/np.mean(od)
    ra,dec = thphi2radec(th,phi)
    plt.scatter(ra,np.sin(dec*np.pi/180),c=od,s=sz,vmax=vx)#,vmin=1.,vmax=2)
    plt.show()


    
    