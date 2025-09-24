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
import astropy.io.fits as fits
import healpy as hp

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

true_selection = colorcuts_function(gflux=gflux, rflux=rflux, zflux=zflux, w1flux=w1flux, w2flux=w2flux,south=True)
true_mean=np.mean(true_selection.astype(float))
print(true_mean)



grand = np.random.normal(size=gflux.shape)
rrand = np.random.normal(size=rflux.shape)
zrand = np.random.normal(size=zflux.shape)
R_G=3.214 # http://legacysurvey.org/dr8/catalogs/#galactic-extinction-coefficients
R_R=2.165
R_Z=1.211



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

def ELGeffcalcExt(gsig,rsig,zsig,ebv,south=True,snrc=True,zmin=-1,zmax=20):
    '''
    calculate the ELG efficiency for given g,r,z flux uncertainties and a given region's selection
    '''
    wz = (photz > zmin) & (photz <= zmax)
    wtg = 10.**(-0.4*R_G*ebv)
    wtr = 10.**(-0.4*R_R*ebv)
    wtz = 10.**(-0.4*R_Z*ebv)
    mgflux = gflux[wz]*wtg + grand[wz]*gsig
    mrflux = rflux[wz]*wtr + rrand[wz]*rsig
    mzflux = zflux[wz]*wtz + zrand[wz]*zsig
    
    combined_snr = np.sqrt(mgflux**2/gsig**2+mrflux**2/rsig**2+mzflux**2/zsig**2) #combined signal to noise; tractor ignores anything < 6
    snrcut = ( combined_snr > 6 ) * ( mgflux > 4*gsig ) * ( mrflux > 3.5*rsig ) * ( mzflux > 2.*zsig )
    selection = colorcuts_function(gflux=mgflux/wtg, rflux=mrflux/wtr, zflux=mzflux/wtz, w1flux=w1flux, w2flux=w2flux, south=south) 
    if snrc:
        selection *= snrcut#( combined_snr > 6 ) * ( mgflux > 4*gsig ) * ( mrflux > 3.5*rsig ) * ( mzflux > 2.*zsig )
        
    efficiency=np.mean(selection.astype(float))/true_mean
    #efficiency_of_true_elgs=np.mean((true_selection*selection).astype(float))/true_mean
    #efficiency_of_stars=np.mean((star_selection*selection).astype(float))/true_mean  
    
    return efficiency#,efficiency_of_true_elgs,efficiency_of_stars

def alt():
    selection_snr   = np.zeros_like(mgflux, dtype=bool)
    snrg = mgflux/gsig
    selection_snr = selection_snr | (snrg > 6.)
    snrr = mrflux/rsig
    selection_snr = selection_snr | (snrr > 6.)
    snrz = mzflux/zsig
    selection_snr = selection_snr | (snrz > 6.)
    ivars = [1./gsig**2,1./rsig**2,1./zsig**2]
    combined_snr = np.sqrt(mgflux**2/gsig**2+mrflux**2/rsig**2+mzflux**2/zsig**2) #combined signal to noise; tractor ignores anything < 6
    selection_snr selection_snr | (combined_snr > 6)
    redmap = mgflux**2/(gsig)**2/2.5+mrflux**2/rsig**2+mzflux**2/(zsig)**2/0.4
    sediv = 1./(gsig*2.5)**2+1./rsig**2+1./(zsig*0.4)**2
    sedmap   /= np.maximum(1.e-16, sediv)
    combined_snrred = sedmap * np.sqrt(sediv) #combined signal to noise; red sed
    selection_snr = selection_snr | (combined_snrred>6.)
    selection = colorcuts_function(gflux=mgflux, rflux=mrflux, zflux=mzflux, w1flux=w1flux, w2flux=w2flux, south=south) 
    if snrc:
        selection *= selection_snr#( combined_snr > 6 ) * ( mgflux > 4*gsig ) * ( mrflux > 3.5*rsig ) * ( mzflux > 2.*zsig )


def thphi2radec(theta,phi):
    return 180./np.pi*phi,-(180./np.pi*theta-90)


def mkeffmap(south=True):
	#read in DR8 properties map
	pixfn      = '/project/projectdirs/desi/target/catalogs/dr8/0.31.1/pixweight/pixweight-dr8-0.31.1.fits'
	pix,header=fitsio.read(pixfn,header=True)
	HPXNSIDE=header["HPXNSIDE"]
	print(pix.dtype.names)
	print("number of pixels=",pix.size)
	ii=np.where((pix["GALDEPTH_G"]>0)&(pix["GALDEPTH_R"]>0)&(pix["GALDEPTH_Z"]>0)&(pix["FRACAREA"]>0.01))[0]
	npix=ii.size
	print(npix)
	pix = pix[ii]
	print(len(pix))
	th,phi = hp.pix2ang(HPXNSIDE,pix['HPXPIXEL'])
	ra,dec = thphi2radec(th,phi)
	depth_keyword="PSFDEPTH"
	#gdepth=-2.5*np.log10(5/np.sqrt(pix[depth_keyword+"_G"]))+22.5
	#rdepth=-2.5*np.log10(5/np.sqrt(pix[depth_keyword+"_R"]))+22.5
	#zdepth=-2.5*np.log10(5/np.sqrt(pix[depth_keyword+"_Z"]))+22.5
	gsigma=1./np.sqrt(pix[depth_keyword+"_G"])
	rsigma=1./np.sqrt(pix[depth_keyword+"_R"])
	zsigma=1./np.sqrt(pix[depth_keyword+"_Z"])
	efficiency =np.zeros(npix)
	#efficiency2 =np.zeros(npix)
	#efficiency3 =np.zeros(npix)
	for j in range(0,len(pix)):# loop on sky pixels using funciton

		if j%5000==0 : print("{}/{}, {:3.2f}%".format(j,ii.size,float(j)/ii.size*100.))
	
		gsig = gsigma[j]#*10**(0.4*R_G*pix["EBV"][j])
		rsig = rsigma[j]#*10**(0.4*R_R*pix["EBV"][j])
		zsig = zsigma[j]#*10**(0.4*R_Z*pix["EBV"][j])
		ebv = pix["EBV"][j]
		#do for three redshift ranges and for overall
		eff = ELGeffcalcExt(gsig,rsig,zsig,ebv,south=True)    
		efficiency[j]=eff
	outf = os.getenv('SCRATCH')+'/ELGMCeffHSCHPext.fits'
	
	collist = []
	collist.append(fits.Column(name='HPXPIXEL',format='K',array=pix['HPXPIXEL']))
	collist.append(fits.Column(name='EFFALL',format='D',array=efficiency))
	hdu  = fits.BinTableHDU.from_columns(fits.ColDefs(collist))
	hdu.writeto(outf,overwrite=True)
		
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