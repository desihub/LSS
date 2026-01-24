import fitsio
from astropy.table import Table,join,vstack
import numpy as np
import glob
import os
import sys
from desitarget import cuts
from desitarget import targetmask
from desitarget.sv1 import sv1_targetmask
from desitarget.sv3 import sv3_targetmask
import astropy.io.fits as fits


def mask(dd,mb=[1,5,6,7,11,12,13]):
    keep = (dd['NOBS_G']>0) & (dd['NOBS_R']>0) & (dd['NOBS_Z']>0)
    print(len(dd[keep]))
    
    keepelg = keep
    for bit in mb:
        keepelg &= ((dd['MASKBITS'] & 2**bit)==0)
    print(len(dd[keepelg]))
    dd = dd[keepelg] 
    return dd       


def gather_targets(type,targroot,outdir,tarver,survey,prog='dark',keys=[]):
	#just concatenate all of the targets for a given type, keeping only the columns quoted below
	print(targroot+prog)
	fns = glob.glob(targroot+prog+'/*.fits')
	ncat     = len(fns)
	print('data is split into '+str(ncat)+' healpix files')
		#check to make sure those were copied correctly
	f = fitsio.read(fns[0])
	for key in keys:
	   try:
		   d = f[key]
	   except:
		   print(key+' not in target file!')
	if survey == 'main':
	    if type[:3] == 'BGS':
	        bs = targetmask.bgs_mask[type]
	        tp = 'BGS_TARGET'
	    else:
	        bs = targetmask.desi_mask[type]
	        tp = 'DESI_TARGET'  
	    ws = '' 
	if survey == 'sv1':
	    bs = sv1_targetmask.desi_mask[type]
	    tp = 'SV1_DESI_TARGET'
	    ws = 'sv1'
	if survey == 'sv3':
	    bs = sv3_targetmask.desi_mask[type]
	    tp = 'SV3_DESI_TARGET'
	    ws = 'sv3'
	print(type+' selection bit is '+str(bs))
	
	outf = outdir+type+ws +'targetsDR9v'+tarver.strip('.')+'.fits'   
	print('file will be written to '+outf)  
	
	data = fitsio.read(fns[0],columns=keys)
	data = data[(data[tp] & bs)>0]
	for i in range(1,ncat):
	    print(i)
	    datan = fitsio.read(fns[i],columns=keys)
	    datan = datan[(datan[tp] & bs)>0]
	    data = np.hstack((data,datan))
	    print(len(data))
	
	
	fitsio.write(outf,data,clobber=True)
	print('wrote to '+outf)
	del data
	#return outf

def gather_targets(type,targroot,outf,tarver,survey,prog='dark',keys=None):
	#just concatenate all of the targets for a given type, keeping only the columns quoted below
	print(targroot+prog)
	fns = glob.glob(targroot+prog+'/*.fits')
	ncat     = len(fns)
	print('data is split into '+str(ncat)+' healpix files')
		#check to make sure those were copied correctly
	f = fitsio.read(fns[0])
	if keys is not None:
		for key in keys:
		   try:
			   d = f[key]
		   except:
			   print(key+' not in target file!')
	else:
	    keys = list(f.dtype.names)
	if survey == 'main':
	    if type[:3] == 'BGS' and 'extra' not in targroot:
	        bs = targetmask.bgs_mask[type]
	        tp = 'BGS_TARGET'
	    else:
			print(type,targroot)
	        bs = targetmask.desi_mask[type]
	        tp = 'DESI_TARGET'  
	    ws = '' 
	if survey == 'sv1':
	    bs = sv1_targetmask.desi_mask[type]
	    tp = 'SV1_DESI_TARGET'
	    ws = 'sv1'
	if survey == 'sv3':
	    bs = sv3_targetmask.desi_mask[type]
	    tp = 'SV3_DESI_TARGET'
	    ws = 'sv3'
	print(type+' selection bit is '+str(bs))
	
	#outf = outdir+type+ws +'targetsDR9v'+tarver.strip('.')+'.fits'   
	print('file will be written to '+outf)  
	
	data = fitsio.read(fns[0],columns=keys)
	data = data[(data[tp] & bs)>0]
	for i in range(1,ncat):
	    print(i)
	    datan = fitsio.read(fns[i],columns=keys)
	    datan = datan[(datan[tp] & bs)>0]
	    data = np.hstack((data,datan))
	    print(len(data))
	
	
	fitsio.write(outf,data,clobber=True)
	print('wrote to '+outf)
	del data

	
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
  


    
def putstar_me(gfluxmin,sfs,sfn,south=True):
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
 
def puttype(type,sfs,sfn,south=True,ebvfac=1.,Rv=3.1):
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

  


  
    
