#commented packages are not used, yet
import sys,os
import numpy as np
import matplotlib.pyplot as plt
#from pkg_resources import resource_filename
#from desiutil.log import get_logger
from desitarget import cuts
#import astropy.io.fits as pyfits
import fitsio
#import healpy as hp

colorcuts_function = cuts.isELG_colors

#deep DECaLS imaging, with photozs from HSC
truthf = '/global/cscratch1/sd/raichoor/desi_mcsyst/desi_mcsyst_truth.dr7.34ra38.-7dec-3.fits' 
truth  = fitsio.read(truthf,1)
gmag  = truth["g"]
w = gmag < 24.5
truth = truth[w]
gmag  = truth["g"]
rmag  = truth["r"]
zmag  = truth["z"]
photz = truth['hsc_mizuki_photoz_best']

def mag2flux(mag) :
    return 10**(-0.4*(mag-22.5))

def flux2mag(flux) :
    mag = -2.5*np.log10(flux*(flux>0)+0.001*(flux<=0)) + 22.5
    mag[(flux<=0)] = 0.
    return mag

gflux = mag2flux(truth["g"])
rflux = mag2flux(truth["r"])
zflux = mag2flux(truth["z"])
w1flux = np.zeros(gflux.shape)#WISE not used in ELG selection, but still needed for code
w2flux = np.zeros(gflux.shape)
grand = np.random.normal(size=gflux.shape)
rrand = np.random.normal(size=rflux.shape)
zrand = np.random.normal(size=zflux.shape)


def ELGsel(gsig,rsig,zsig,south=True,snrc=True):
    '''
    calculate the ELG selection for given g,r,z flux uncertainties and a given region's selection
    snrc can be toggled for the total signal to noise cut; tests suggest it is indeed necessary
    '''
    mgflux = gflux + grand*gsig
    mrflux = rflux + rrand*rsig
    mzflux = zflux + zrand*zsig
    
    combined_snr = np.sqrt(mgflux**2/gsig**2+mrflux**2/rsig**2+mzflux**2/zsig**2) #combined signal to noise; tractor ignores anything < 6
    
    selection = colorcuts_function(gflux=mgflux, rflux=mrflux, zflux=mzflux, w1flux=w1flux, w2flux=w2flux, south=south) 
    if snrc:
        selection *= ( combined_snr > 6 ) * ( mgflux > 4*gsig ) * ( mrflux > 3.5*rsig ) * ( mzflux > 2.*zsig )
    return selection

selmed = ELGsel(0.023,0.041,.06) #slightly worse than median, one could improve this

def getrelnz(sigg,sigr,sigz,bs=0.05,zmin=0.6,zmax=1.4,south=True):
	nb = int((zmax*1.001-zmin)/bs)
	print(nb)
	zbe = np.linspace(zmin,zmax,nb+1) #bin edges for histogram, so nb +1
	zl = np.arange(zmin+bs/2.,zmax,bs) #bin centers
	nzmed = np.histogram(photz[selmed],range=(zmin,zmax),bins=zbe)
	sel = ELGsel(sigg,sigr,sigz,south=south)
	nztest = np.histogram(photz[sel],range=(zmin,zmax),bins=zbe)
	nzrel = nztest[0]/nzmed[0]
	
	return zl,nzrel

if __name__ == '__main__':
	'''
	For now, just a test to show it works
	'''
	
	sigg,sigr,sigz = float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3])
	zl,nzrel = getrelnz(sigg,sigr,sigz)
	plt.plot(zl,nzrel,'k-')
	plt.show()
	'''
	AJR suggests taking this curve, fitting a smooth function to it, and then using that to modulate the mocks
	'''