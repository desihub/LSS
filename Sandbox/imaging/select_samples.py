import fitsio
from astropy.table import Table,join,vstack
import numpy as np
import glob
import os
import sys
from desitarget import cuts
from desitarget import targetmask
import astropy.io.fits as fits

dirsweeps = '/global/project/projectdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/'
dirsweepn = '/global/project/projectdirs/cosmo/data/legacysurvey/dr9/north/sweep/9.0/'
targroot = '/project/projectdirs/desi/target/catalogs/dr9m/0.44.0/targets/main/resolve/'

sfs = glob.glob(dirsweeps+'sweep*')
sfn = glob.glob(dirsweepn+'sweep*')

outdir = '/project/projectdirs/desi/users/ajross/dr9/'

def gather_targets(type,fo='targetDR9m44.fits',prog='dark'):
	#just concatenate all of the targets for a given type, keeping only the columns quoted below
	print(targroot+prog)
	fns = glob.glob(targroot+prog+'/*.fits')
	ncat     = len(fns)
	print('data is split into '+str(ncat)+' healpix files')
	keys = ['RA', 'DEC', 'BRICKID', 'BRICKNAME','MORPHTYPE','DCHISQ','FLUX_G', 'FLUX_R', 'FLUX_Z','MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z','FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z','NOBS_G', 'NOBS_R', 'NOBS_Z','PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'GALDEPTH_G', 'GALDEPTH_R',\
	   'GALDEPTH_Z','FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z',\
	   'MASKBITS', 'EBV', 'PHOTSYS','TARGETID','DESI_TARGET']
	#check to make sure those were copied correctly
	f = fitsio.read(fns[0])
	for key in keys:
	   try:
		   d = f[key]
	   except:
		   print(key+' not in target file!')
	bs = targetmask.desi_mask[type]  
	
	outf = outdir+type +fo   
	print('file will be written to '+outf)  
	
	data = fitsio.read(fns[0],columns=keys)
	data = data[(data['DESI_TARGET'] & bs)>0]
	for i in range(1,ncat):
	    print(i)
	    datan = fitsio.read(fns[i],columns=keys)
	    datan = datan[(datan['DESI_TARGET'] & bs)>0]
	    data = np.hstack((data,datan))
	    print(len(data))
	
	
	fitsio.write(outf,data,clobber=True)
	print('wrote to '+outf)
	del data
# 	put information together, takes a couple of minutes
# 	
# 	mydict   = {}
# 	for key in keys:
# 		mydict[key] = []
# 	
# 	for i in range(0,ncat):
# 		data = fitsio.read(fns[i],columns=keys)
# 		data = data[(data['DESI_TARGET'] & bs)>0]
# 		for key in keys:
# 			mydict[key] += data[key].tolist()
# 		print(i)
# 		
# 	print(str(len(mydict['RA']))+' '+type)
# 	
# 	 
# 	collist = []
# 	for key in keys:
# 		fmt = fits.open(fns[0])[1].columns[key].format
# 		collist.append(fits.Column(name=key,format=fmt,array=mydict[key]))
# 		print(key)
# 	hdu  = fits.BinTableHDU.from_columns(fits.ColDefs(collist))
# 	hdu.writeto(outf,overwrite=True)
# 	print('wrote to '+outf)
# 	del collist
# 	del mydict
	
def starsel_sweep(f,gfluxmin):
    w = f['TYPE'] == 'PSF '
    gflux = f['FLUX_G']/f['MW_TRANSMISSION_G']
    w &= gflux > gfluxmin
    return w

def typesel(f,type,south=True,ebvfac=1.,Rv=3.1):
	if ebvfac == 1. and Rv == 3.1:
		gflux = f['FLUX_G']/f['MW_TRANSMISSION_G'] 
		rflux = f['FLUX_R']/f['MW_TRANSMISSION_R']   
		zflux = f['FLUX_Z']/f['MW_TRANSMISSION_Z'] 
		w1flux = f['FLUX_W1']/f['MW_TRANSMISSION_W1']
		zfiberflux = f['FIBERFLUX_Z']/f['MW_TRANSMISSION_Z'] 
	else:
		Rg = 3.214*ebvfac
		Rr = 2.165*ebvfac
		Rz = 1.211*ebvfac
		Rw1 = 0.184*ebvfac
		if Rv < 3.1:
			#linear interpolation from Schlafly 2011 table
			Rg = (3.739-3.273)*(3.1-Rv)*ebvfac+Rg
			Rr = (2.113-2.176)*(3.1-Rv)*ebvfac+Rr	 
			Rz = (1.175-1.217)*(3.1-Rv)*ebvfac+Rz
			Rw1 = (-.1)*(Rv-3.1)*ebvfac+Rw1
		if Rv > 3.1:
			#linear interpolation from Schlafly 2011 table
			Rg = (3.006-3.273)*(Rv-3.1)*ebvfac+Rg
			Rr = (2.205-2.176)*(Rv-3.1)*ebvfac+Rr	 
			Rz = (1.236-1.217)*(Rv-3.1)*ebvfac+Rz
			Rw1 = (-.05)*(Rv-3.1)*ebvfac+Rw1
		print('ebvfac,Rv,Rg,Rr,Rz,Rw1')
		print(ebvfac,Rv,Rg,Rr,Rz,Rw1)
		wtg = 10**(-0.4*Rg*f['EBV'])
		wtr = 10**(-0.4*Rr*f['EBV'])
		wtz = 10**(-0.4*Rz*f['EBV'])
		wtw = 10**(-0.4*Rw1*f['EBV'])
		gflux = f['FLUX_G']/wtg
		rflux = f['FLUX_R']/wtr 
		zflux = f['FLUX_Z']/wtz
		w1flux = f['FLUX_W1']/wtw
		zfiberflux = f['FIBERFLUX_Z']/wtz
	if type == 'LRG':
		w = cuts.isLRG_colors(gflux, rflux, zflux, w1flux,zfiberflux, south=south)
	if type == 'ELG':
		w = cuts.isELG_colors(gflux, rflux, zflux, w1flux,zfiberflux, south=south)
	return w
  


    
def putstar_me(gfluxmin,south=True):
    if south == True:
        fls = sfs
    else:
        fls = sfn
    f = fitsio.read(fls[0])
    w0 = starsel_sweep(f,gfluxmin)
    fw = f[w0]
    #dt = fw.dtype
    dt = []
    cols = ['BRICKID','OBJID','RA','DEC','DCHISQ','EBV','FLUX_G','FLUX_R','FLUX_Z','MW_TRANSMISSION_G','MW_TRANSMISSION_R','MW_TRANSMISSION_Z',\
    'PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GAIA_PHOT_G_MEAN_MAG','GAIA_ASTROMETRIC_EXCESS_NOISE','MASKBITS']
    for colname in cols:
    	dt.append((colname,fw.dtype[colname]))
    new = np.empty(len(fw),dtype=dt)
    for colname in cols:
        new[colname][...] = fw[colname][...]
    tm = new       
    #for i in range(1,10):
    for i in range(1,len(fls)):
        #try:
        f = fitsio.read(fls[i])
        w = starsel_sweep(f,gfluxmin)
        fw = f[w]
        new = np.empty(len(fw),dtype=dt)
        for colname in cols:
            new[colname][...] = fw[colname][...]
        tmn = np.concatenate((tm, new),axis=0)
        print(i,len(tmn))
        tm = tmn
        #except:
        #    print(i)
    NS = 'south'
    if south != True:
        NS = 'north'
    
    s = 0
    es = ''
    while s == 0:
        outf = outdir+'mysweeps/stars_gfluxg'+str(gfluxmin)+'_'+NS+es+'.fits'
        try:
            fitsio.read(outf)
            es += 'n'
            print(es)
            if len(es) > 10:
            	return 'es too long, probably a bug'
        except:
            s = 1
        
    #tm = Table(tm,names=fi.dtype.names)
    fits = fitsio.FITS(outf,'rw')
    fits.write(tm, names=cols,overwrite=True)    
 
def puttype(type,south=True,ebvfac=1.,Rv=3.1):
    if south == True:
        fls = sfs
    else:
        fls = sfn
    f = fitsio.read(fls[0])
    w0 = typesel(f,type,south=south,ebvfac=ebvfac,Rv=Rv)
    fw = f[w0]
    dt = fw.dtype
    new = np.empty(len(fw),dtype=dt)
    #cols = ['BRICKID','OBJID','RA','DEC','DCHISQ','EBV','FLUX_G','FLUX_R','FLUX_Z','MW_TRANSMISSION_G','MW_TRANSMISSION_R','MW_TRANSMISSION_Z',\
    #'PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GAIA_PHOT_G_MEAN_MAG','GAIA_ASTROMETRIC_EXCESS_NOISE','MASKBITS']
    #for colname in fw.dtype.names:
    cols = fw.dtype.names
    for colname in cols:
        new[colname][...] = fw[colname][...]
    tm = new       
    for i in range(1,len(fls)):
        #try:
        f = fitsio.read(fls[i])
        w = typesel(f,type,south=south,ebvfac=ebvfac,Rv=Rv)
        fw = f[w]
        new = np.empty(len(fw),dtype=dt)
        for colname in cols:
            new[colname][...] = fw[colname][...]
        tmn = np.concatenate((tm, new),axis=0)
        print(i,len(tmn))
        tm = tmn
        #except:
        #    print(i)
    NS = 'south'
    if south != True:
        NS = 'north'
    
    s = 0
    es = 'ebvfac'+str(ebvfac)+'Rv'+str(Rv)
    outf = outdir+'mysweeps/'+type+'dr8_'+NS+es+'.fits'
#     while s == 0:
#         
#         try:
#             fitsio.read(outf)
#             es += 'n'
#             print(es)
#             if len(es) > 10:
#             	return 'es too long, probably a bug'
#         except:
#             s = 1
        
    #tm = Table(tm,names=fi.dtype.names)
    fits = fitsio.FITS(outf,'rw')
    fits.write(tm, names=dt.names,overwrite=True)    

  


  
    
if __name__ == '__main__':
	#get target catalogs for the main dark time programs
	gather_targets('ELG')
	gather_targets('LRG')
	gather_targets('QSO')
	#fluxlim = float(str(sys.argv[1]))
	#putstar(fluxlim)
	#this was meant to test dependence with changes to extinction parameters explicitly
	#ebf = 1
	#rv = 2.6
	#puttype('LRG',ebvfac=ebf,Rv=rv)
	#print('LRG south done') 
	#puttype('LRG',south=False,ebvfac=ebf,Rv=rv)  
	#print('LRG north done')
	#puttype('ELG',ebvfac=ebf,Rv=rv)
	#print('ELG south done') 
	#puttype('ELG',south=False,ebvfac=ebf,Rv=rv)  
	#print('ELG north done')	