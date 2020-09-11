import sys,os
import numpy as np
import matplotlib.pyplot as plt
from desitarget import cuts
import fitsio
import astropy.io.fits as fits
import healpy as hp

colorcuts_function = cuts.isELG_colors

#deep DECaLS imaging, with photozs from HSC
truthf = '/project/projectdirs/desi/users/ajross/MCdata/desi_mcsyst_truth.dr7.34ra38.-7dec-3.fits' 
truth  = fitsio.read(truthf,1)
gmag  = truth["g"]
w = gmag < 24.5
truth = truth[w]
gmag  = truth["g"]
rmag  = truth["r"]
zmag  = truth["z"]
photz = truth['hsc_mizuki_photoz_best']

pixfn      = '/project/projectdirs/desi/target/catalogs/dr8/0.31.1/pixweight/pixweight-dr8-0.31.1.fits' #update this to be more recent


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

#set up correlation matrix for fluxes

ml = np.zeros(3)
cv = np.ones((3,3))*.5 #just given them all correlation of 0.5 for now
cv[0][0] = 1.
cv[1][1] = 1.
cv[2][2] = 1.

cg = np.random.default_rng().multivariate_normal(ml,cv,len(gflux))
cg = cg.transpose()





def ELGeffcalcExt(gsig,rsig,zsig,wtg,wtr,wtz,south=True,snrc=True,zmin=-1,zmax=20,corr=True,gf=1.,rf=1.,zf=1.,rsel=False):
    '''
    calculate the ELG efficiency for given g,r,z flux uncertainties and a given region's selection
    gsig, rsig, zsig are 1sigma flux uncertainties for g,r,z
    wtg,wtr,wtz are Milky Way transmission coefficients (i.e. Galactic extinction < 1 multiplied by flux to account for loss)
    South toggles whether north or south target selection cuts get used (truth data is DECaLS, so maybe should always be south until that is updated)
    zmin,zmax control redshift range of photozs from truth
    corr toggles whether or not correlation is assumed between flux measurements
    gf,rf,zf allow one to test what happens if the flux is multiplied by these factors
    rsel toggles whether the selection or the efficiency is returned
    '''
    wz = (photz > zmin) & (photz <= zmax)
    if corr:
        mgflux = gflux[wz]*wtg*gf + cg[0][wz]*gsig
        mrflux = rflux[wz]*wtr*rf + cg[1][wz]*rsig
        mzflux = zflux[wz]*wtz*zf + cg[2][wz]*zsig
    else:
        mgflux = gflux[wz]*wtg*gf + grand[wz]*gsig
        mrflux = rflux[wz]*wtr*rf + rrand[wz]*rsig
        mzflux = zflux[wz]*wtz*zf + zrand[wz]*zsig
   
    selection = colorcuts_function(gflux=mgflux/wtg, rflux=mrflux/wtr, zflux=mzflux/wtz, w1flux=w1flux, w2flux=w2flux, south=south) 
    selection_snr   = np.zeros_like(mgflux, dtype=bool)
    snrg = mgflux/gsig
    selection_snr = selection_snr | (snrg > 6.)
    snrr = mrflux/rsig
    selection_snr = selection_snr | (snrr > 6.)
    snrz = mzflux/zsig
    selection_snr = selection_snr | (snrz > 6.)

    flatmap = mgflux/(gsig)**2+mrflux/(rsig)**2+mzflux/(zsig)**2
    fdiv = 1./(gsig)**2+1./rsig**2+1./(zsig)**2
    flatmap   /= np.maximum(1.e-16, fdiv)
    combined_snr = flatmap * np.sqrt(fdiv) #combined signal to noise matching Dustin's vode for flat sed

    selection_snr = selection_snr | (combined_snr > 6)
    redmap = mgflux/(gsig)**2/2.5+mrflux/rsig**2+mzflux/(zsig)**2/0.4
    sediv = 1./(gsig*2.5)**2+1./rsig**2+1./(zsig*0.4)**2
    redmap   /= np.maximum(1.e-16, sediv)
    combined_snrred = redmap * np.sqrt(sediv) #combined signal to noise; red sed
    selection_snr = selection_snr | (combined_snrred>6.)
    
    if snrc:
        selection *= selection_snr
    
    if rsel:
        return selection #just return the selection if rsel is True    
    
    efficiency=np.mean(selection.astype(float))/true_mean
    
    return efficiency

selmed = ELGeffcalcExt(0.023,0.041,.06,1.,1.,1.,rsel=True) #slightly worse than median, no extinction, one could improve this

def getrelnz(sigg,sigr,sigz,wtg,wtr,wtz,bs=0.05,zmin=0.6,zmax=1.4,south=True):
    nb = int((zmax*1.001-zmin)/bs)
    print(nb)
    zbe = np.linspace(zmin,zmax,nb+1) #bin edges for histogram, so nb +1
    zl = np.arange(zmin+bs/2.,zmax,bs) #bin centers
    nzmed = np.histogram(photz[selmed],range=(zmin,zmax),bins=zbe)
    sel = ELGeffcalcExt(sigg,sigr,sigz,wtg,wtr,wtz,south=south,rsel=True)
    nztest = np.histogram(photz[sel],range=(zmin,zmax),bins=zbe)
    nzrel = nztest[0]/nzmed[0]
    
    return zl,nzrel


def thphi2radec(theta,phi):
    return 180./np.pi*phi,-(180./np.pi*theta-90)


def mkeffmap(south=True,outf='/ELGMCeffHSCHPext.fits'):
    #read in DR8 properties map
    
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
    gsigma=1./np.sqrt(pix[depth_keyword+"_G"])
    rsigma=1./np.sqrt(pix[depth_keyword+"_R"])
    zsigma=1./np.sqrt(pix[depth_keyword+"_Z"])
    efficiency =np.zeros(npix)
    wtgp = 10.**(-0.4*R_G*pix['EBV'])
    wtrp = 10.**(-0.4*R_R*pix['EBV'])
    wtzp = 10.**(-0.4*R_Z*pix['EBV'])
    for j in range(0,len(pix)):# loop on sky pixels using funciton

        if j%5000==0 : print("{}/{}, {:3.2f}%".format(j,ii.size,float(j)/ii.size*100.))
    
        gsig = gsigma[j]
        rsig = rsigma[j]
        zsig = zsigma[j]
        ebv = pix["EBV"][j]
        #do for three redshift ranges and for overall
        eff = ELGeffcalcExt(gsig,rsig,zsig,wtgp[j],wtrp[j],wtzp[j],south=True)    
        efficiency[j]=eff
    outf = os.getenv('SCRATCH')+outf
    
    collist = []
    collist.append(fits.Column(name='HPXPIXEL',format='K',array=pix['HPXPIXEL']))
    collist.append(fits.Column(name='EFFALL',format='D',array=efficiency))
    hdu  = fits.BinTableHDU.from_columns(fits.ColDefs(collist))
    hdu.writeto(outf,overwrite=True)
        
if __name__ == '__main__':
    '''
    For now, nothing
    '''
    
    print('nothing in main')