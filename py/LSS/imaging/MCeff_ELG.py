import sys,os
import numpy as np
import matplotlib.pyplot as plt
from desitarget import cuts
import fitsio
import astropy.io.fits as fits
import healpy as hp
from scipy.special import erf
from astropy.table import Table

colorcuts_function = cuts.isELG_colors

#deep DECaLS imaging, with photozs from HSC
truthf = '/project/projectdirs/desi/users/ajross/MCdata/desi_mcsyst_truth.dr7.34ra38.-7dec-3.fits' 
truth  = fitsio.read(truthf,1)
gmag  = truth["g"]
w = gmag < 24.5
#truth = truth[w]
gmag  = truth["g"]
rmag  = truth["r"]
zmag  = truth["z"]
photz = truth['hsc_mizuki_photoz_best']

pixfn = '/global/cfs/cdirs/desi/survey/catalogs/pixweight_maps_all/pixweight-1-dark.fits' 


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


def perturb_flux(ina,outf='test.fits'):
	'''
	ina should be input array containing necessary columns
	the idea here is that input photometry + flux errors and their cov given by cv an output distribution consistent with Obiwan could be produced
	'''

	vv = np.zeros(3)
	cc = np.ones((3,3))
	cc[0][0] = 1.86
	cc[1][1] = 1.75
	cc[2][2] = 1.64
	cc[0][1] = 0.643
	cc[1][0] = cc[0][1]
	cc[0][2] = 0.321
	cc[2][0] = 0.321
	cc[1][2] = 0.341
	cc[2][1] = cc[1][2]
	pg = np.random.default_rng().multivariate_normal(vv,cc,len(ina)) #this provides correlated vectors for perturbing fluxes
	gflux = ina['input_flux_g'] #column name from Obiwan file
	rflux = ina['input_flux_r'] #column name from Obiwan file
	zflux = ina['input_flux_z'] #column name from Obiwan file
	wtg = ina['input_mw_transmission_g']
	wtr = ina['input_mw_transmission_r']
	wtz = ina['input_mw_transmission_z']
	gsig = (1.35/ina['galdepth_g'])**.5 #factors are based on ivar/galdepth from obiwan output
	rsig = (1.44/ina['galdepth_r'])**.5
	zsig = (1.66/ina['galdepth_z'])**.5

	mgflux = gflux*wtg + pg[0]*gsig
	mrflux = rflux*wtr + pg[1]*rsig
	mzflux = zflux*wtz + pg[2]*zsig

	snrg = mgflux/gsig    
	snrr = mrflux/rsig    
	snrz = mzflux/zsig

	flatmap = mgflux/(gsig)**2+mrflux/(rsig)**2+mzflux/(zsig)**2
	fdiv = 1./(gsig)**2+1./rsig**2+1./(zsig)**2
	flatmap   /= np.maximum(1.e-16, fdiv)
	combined_snr2 = flatmap**2.*fdiv

	redmap = mgflux/(gsig)**2/2.5+mrflux/rsig**2+mzflux/(zsig)**2/0.4
	sediv = 1./(gsig*2.5)**2+1./rsig**2+1./(zsig*0.4)**2
	redmap   /= np.maximum(1.e-16, sediv)
	combined_snrred2 = redmap**2. * (sediv) 

	to = Table([gflux,rflux,zflux,wtg,wtr,wtz,ina['galdepth_g'],ina['galdepth_r'],ina['galdepth_z'],snrg,snrr,snrz,combined_snr2,combined_snrred2],\
			   names=('input_flux_g','input_flux_r','input_flux_z','input_mw_transmission_g','input_mw_transmission_r','input_mw_transmission_z','galdepth_g','galdepth_r','galdepth_z','snr_g','snr_r','snr_z','combined_snr2','combined_snrred2'))
	to.write(outf,format='fits',overwrite=True)  
	return True         
    
    


def ELGeffcalcExt(gsig,rsig,zsig,wtg,wtr,wtz,south=True,snrc=True,zmin=-1,zmax=20,corr=True,gf=1.,rf=1.,zf=1.,dg=0,dr=0,dz=0,sg=0,gfluxcut=None,rsel=False,vis=False,gefac=0):
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
        mgflux = gflux[wz]*wtg*(gf+(1.-gf)*erf(gflux[wz]*gefac)) + cg[0][wz]*gsig+dg
        mrflux = rflux[wz]*wtr*rf + cg[1][wz]*rsig+dr
        mzflux = zflux[wz]*wtz*zf + cg[2][wz]*zsig+dz
    else:
        mgflux = gflux[wz]*wtg*gf + grand[wz]*gsig+dg
        mrflux = rflux[wz]*wtr*rf + rrand[wz]*rsig+dr
        mzflux = zflux[wz]*wtz*zf + zrand[wz]*zsig+dz
   
    selection = colorcuts_function(gflux=mgflux/wtg, rflux=mrflux/wtr, zflux=mzflux/wtz, w1flux=w1flux, w2flux=w2flux, south=south) 
    selection_snr   = np.zeros_like(mgflux, dtype=bool)
    snrg = mgflux/gsig    
    snrr = mrflux/rsig    
    snrz = mzflux/zsig
    selection_snr = selection_snr | (snrr > 6.)
    selection_snr = selection_snr | (snrg > 6.)
    selection_snr = selection_snr | (snrz > 6.)

    flatmap = mgflux/(gsig)**2+mrflux/(rsig)**2+mzflux/(zsig)**2
    fdiv = 1./(gsig)**2+1./rsig**2+1./(zsig)**2
    flatmap   /= np.maximum(1.e-16, fdiv)
    #combined_snr = flatmap * np.sqrt(fdiv) #combined signal to noise matching Dustin's vode for flat sed
    combined_snr2 = flatmap**2.*fdiv #faster to remove sqrt?
    #selection_snr = selection_snr | (combined_snr > 6)
    #selection_snr = selection_snr | (combined_snr2 > 36)
    redmap = mgflux/(gsig)**2/2.5+mrflux/rsig**2+mzflux/(zsig)**2/0.4
    sediv = 1./(gsig*2.5)**2+1./rsig**2+1./(zsig*0.4)**2
    redmap   /= np.maximum(1.e-16, sediv)
    #combined_snrred = redmap * np.sqrt(sediv) #combined signal to noise; red sed
    combined_snrred2 = redmap**2. * (sediv) #faster to remove sqrt?
    #selection_snr = selection_snr | (combined_snrred>6.)
    #selection_snr = selection_snr | (combined_snrred2>36.)
    selection_snr = selection_snr & ((snrg>0) & (snrr>0) & (snrz > 0))
    if snrc:
        selection *= selection_snr

    if gfluxcut:
        selg = mgflux/wtg > gfluxcut
        selection *= selg        
        
    if rsel:
        return selection #just return the selection if rsel is True    
    
    efficiency=np.mean(selection.astype(float))/true_mean
    
    if vis:
        plt.hist(snrg[selection],bins=100,range=(0,15),label='g',histtype='step')
        plt.xlabel('S/N')
        plt.hist(snrr[selection],bins=100,range=(0,15),label='r',histtype='step')
        #plt.xlabel('S/N')
        plt.hist(snrz[selection],bins=100,range=(0,15),label='z',histtype='step')
        #plt.xlabel('S/N')
        plt.legend()
        plt.show()
        plt.hist(mgflux[selection],bins=100,range=(0,2),label='g',histtype='step')
        plt.xlabel('flux')
        plt.hist(mrflux[selection],bins=100,range=(0,2),label='r',histtype='step')
        #plt.xlabel('S/N')
        plt.hist(mzflux[selection],bins=100,range=(0,2),label='z',histtype='step')
        #plt.xlabel('S/N')
        plt.legend()
        plt.show()
  
    return efficiency

def ELGeffcalcExt_dect(gsig,rsig,zsig,wtg,wtr,wtz,south=True,zmin=-1,zmax=20,gf=1.,rf=1.,zf=1.,rsel=False):
    '''
    calculate the ELG efficiency for given g,r,z flux uncertainties and a given region's selection
    only consider effect of needing 6sigma detection
    gsig, rsig, zsig are 1sigma flux uncertainties for g,r,z
    wtg,wtr,wtz are Milky Way transmission coefficients (i.e. Galactic extinction < 1 multiplied by flux to account for loss)
    South toggles whether north or south target selection cuts get used (truth data is DECaLS, so maybe should always be south until that is updated)
    zmin,zmax control redshift range of photozs from truth
    corr toggles whether or not correlation is assumed between flux measurements
    gf,rf,zf allow one to test what happens if the flux is multiplied by these factors
    rsel toggles whether the selection or the efficiency is returned
    '''
    
    wz = (photz > zmin) & (photz <= zmax)
    mgflux = gflux[wz]*wtg*gf 
    mrflux = rflux[wz]*wtr*rf 
    mzflux = zflux[wz]*wtz*zf 
   
    selection = colorcuts_function(gflux=mgflux/wtg, rflux=mrflux/wtr, zflux=mzflux/wtz, w1flux=w1flux, w2flux=w2flux, south=south) 
    selection_snr   = np.zeros_like(mgflux, dtype=bool)
    snrg = mgflux/gsig    
    snrr = mrflux/rsig    
    snrz = mzflux/zsig
    selection_snr = selection_snr | (snrr > 6.)
    selection_snr = selection_snr | (snrg > 6.)
    selection_snr = selection_snr | (snrz > 6.)

    flatmap = mgflux/(gsig)**2+mrflux/(rsig)**2+mzflux/(zsig)**2
    fdiv = 1./(gsig)**2+1./rsig**2+1./(zsig)**2
    flatmap   /= np.maximum(1.e-16, fdiv)
    #combined_snr = flatmap * np.sqrt(fdiv) #combined signal to noise matching Dustin's vode for flat sed
    combined_snr2 = flatmap**2.*fdiv #faster to remove sqrt?
    #selection_snr = selection_snr | (combined_snr > 6)
    selection_snr = selection_snr | (combined_snr2 > 36)
    redmap = mgflux/(gsig)**2/2.5+mrflux/rsig**2+mzflux/(zsig)**2/0.4
    sediv = 1./(gsig*2.5)**2+1./rsig**2+1./(zsig*0.4)**2
    redmap   /= np.maximum(1.e-16, sediv)
    #combined_snrred = redmap * np.sqrt(sediv) #combined signal to noise; red sed
    combined_snrred2 = redmap**2. * (sediv) #faster to remove sqrt?
    #selection_snr = selection_snr | (combined_snrred>6.)
    selection_snr = selection_snr | (combined_snrred2>36.)
    selection_snr = selection_snr & ((snrg>0) & (snrr>0) & (snrz > 0))
    
    selection *= selection_snr
    
    if rsel:
        return selection #just return the selection if rsel is True    
    
    efficiency=np.mean(selection.astype(float))/true_mean
    
    return efficiency

def getELGdist(gsig,rsig,zsig,ebv,south=True,zmin=-1,zmax=20,corr=True,gf=1.,rf=1.,zf=1.):
    '''
    get truth and perturbed fluxes for given g,r,z flux uncertainties and a given region's selection
    gsig, rsig, zsig are 1sigma flux uncertainties for g,r,z
    ebv is Milky Way E(B-V) dust extinction
    South toggles whether north or south target selection cuts get used (truth data is DECaLS, so maybe should always be south until that is updated)
    zmin,zmax control redshift range of photozs from truth
    corr toggles whether or not correlation is assumed between flux measurements
    gf,rf,zf allow one to test what happens if the flux is multiplied by these factors
    rsel toggles whether the selection or the efficiency is returned
    '''
    wtg = 10.**(-0.4*R_G*ebv)
    wtr = 10.**(-0.4*R_R*ebv)
    wtz = 10.**(-0.4*R_Z*ebv)
    wz = (photz > zmin) & (photz <= zmax)

    if south == False:
        #shifting the true flux to north from south, using negative exponent compared to https://github.com/desihub/desitarget/blob/master/py/desitarget/cuts.py#L72
        gfluxc = gflux * 10**(0.4*0.004) * (gflux/rflux)**(0.059)
        rfluxc = rflux * 10**(-0.4*0.003) * (rflux/zflux)**(0.024)
        zfluxc = zflux * 10**(-0.4*0.013) * (rflux/zflux)**(-0.015)
    else:
        gfluxc = gflux
        rfluxc = rflux
        zfluxc = zflux

    
    if corr:
        mgflux = gfluxc[wz]*wtg*gf + cg[0][wz]*gsig
        mrflux = rfluxc[wz]*wtr*rf + cg[1][wz]*rsig
        mzflux = zfluxc[wz]*wtz*zf + cg[2][wz]*zsig
    else:
        mgflux = gfluxc[wz]*wtg*gf + grand[wz]*gsig
        mrflux = rfluxc[wz]*wtr*rf + rrand[wz]*rsig
        mzflux = zfluxc[wz]*wtz*zf + zrand[wz]*zsig
   

    selection = colorcuts_function(gflux=mgflux/wtg, rflux=mrflux/wtr, zflux=mzflux/wtz, w1flux=w1flux, w2flux=w2flux, south=south)
    ebvs = np.ones(len(mgflux))*ebv
    gsigs = np.ones(len(mgflux))*gsig
    rsigs = np.ones(len(mgflux))*rsig
    zsigs = np.ones(len(mgflux))*zsig
    arrtot = np.array([gfluxc,rfluxc,zfluxc,mgflux,mrflux,mzflux,ebvs,gsigs,rsigs,zsigs])
    dt = [('True_g_flux', float), ('True_r_flux', float), ('True_z_flux', float),('g_flux', float), ('r_flux', float), ('z_flux', float),('EBV', float),('sigma_g_flux', float), ('sigma_r_flux', float), ('sigma_z_flux', float)]
    arrtot = np.rec.fromarrays(arrtot,dtype=dt) 
    arrtot = arrtot[selection]
    return arrtot

def cutSN(inl):
    selection_snr   = np.zeros_like(inl['g_flux'], dtype=bool)
    mgflux = inl['g_flux']
    mrflux = inl['r_flux']
    mzflux = inl['z_flux']
    gsig = inl['sigma_g_flux']
    rsig = inl['sigma_r_flux']
    zsig = inl['sigma_z_flux']
    snrg = mgflux/gsig   
    snrr = mrflux/rsig       
    snrz = mzflux/zsig  
    selection_snr = selection_snr | (snrr > 6.)
    selection_snr = selection_snr | (snrg > 6.)
    selection_snr = selection_snr | (snrz > 6.)

    flatmap = mgflux/(gsig)**2+mrflux/(rsig)**2+mzflux/(zsig)**2
    fdiv = 1./(gsig)**2+1./rsig**2+1./(zsig)**2
    flatmap   /= np.maximum(1.e-16, fdiv)
    #combined_snr = flatmap * np.sqrt(fdiv) #combined signal to noise matching Dustin's vode for flat sed
    combined_snr2 = flatmap**2.*fdiv #faster to remove sqrt?
    #selection_snr = selection_snr | (combined_snr > 6)
    selection_snr = selection_snr | (combined_snr2 > 36)
    redmap = mgflux/(gsig)**2/2.5+mrflux/rsig**2+mzflux/(zsig)**2/0.4
    sediv = 1./(gsig*2.5)**2+1./rsig**2+1./(zsig*0.4)**2
    redmap   /= np.maximum(1.e-16, sediv)
    #combined_snrred = redmap * np.sqrt(sediv) #combined signal to noise; red sed
    combined_snrred2 = redmap**2. * (sediv) #faster to remove sqrt?
    #selection_snr = selection_snr | (combined_snrred>6.)
    selection_snr = selection_snr | (combined_snrred2>36.)
    selection_snr = selection_snr & ((snrg>0) & (snrr>0) & (snrz > 0))
    
    #selection_snr = selection_snr & (snrg > 6.)
    return inl[selection_snr]


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


def mkeffmap(south=True,outf='/DR9mELGMCeffHSCHPextnocorr.fits',corr=False):
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
        #ebv = pix["EBV"][j]
        #do for three redshift ranges and for overall
        eff = ELGeffcalcExt(gsig,rsig,zsig,wtgp[j],wtrp[j],wtzp[j],south=True,corr=corr)    
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