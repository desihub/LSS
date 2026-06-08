#!/usr/bin/env python

import os
import sys
import numpy as np
from glob import glob
from astropy.io import fits
import fitsio
from time import time
from fiberassign.utils import Logger
import healpy as hp
from astropy.coordinates import SkyCoord
import astropy.units as units
from argparse import ArgumentParser

# AR reading arguments
parser = ArgumentParser()
parser.add_argument('--outdir',  help='output directory',type=str,default=os.getenv('CSCRATCH')+'/dr9m/elg/',required=False,metavar='OUTDIR')
#
args     = parser.parse_args()
log      = Logger.get()

# input directories/files
maintargdir= '/global/cfs/cdirs/desi/target/catalogs/dr9m/0.44.0/targets/main/resolve/dark/'
sv1targdir = '/global/cfs/cdirs/desi/target/catalogs/dr9m/0.44.0/targets/sv1/resolve/dark/'
randdir    = '/global/cfs/cdirs/desi/target/catalogs/dr9m/0.44.0/randoms/resolve/'
sv1pixfn   = '/global/cfs/cdirs/desi/target/catalogs/dr9m/0.44.0/pixweight/sv1/resolve/dark/sv1pixweight-dark.fits'

start    = time()
log.info('{:.1f}s\tstart'.format(time()-start))

def get_isdes(ra,dec):
    hdu = fits.open('/global/cscratch1/sd/raichoor/desits/des_hpmask.fits')
    nside,nest = hdu[1].header['HPXNSIDE'],hdu[1].header['HPXNEST']
    hppix     = hp.ang2pix(nside,(90.-dec)*np.pi/180.,ra*np.pi/180.,nest=nest)
    isdes     = np.zeros(len(ra),dtype=bool)
    isdes[np.in1d(hppix,hdu[1].data['hppix'])] = True
    return isdes

def get_subs(d,gkey):
	if (gkey=='gtot'):
		g = 22.5-2.5*np.log10(d['FLUX_G']/d['MW_TRANSMISSION_G'])
		gfaint_n,gfaint_s = 23.6,23.5
	elif (gkey=='gfib'):
		g = 22.5-2.5*np.log10(d['FIBERFLUX_G']/d['MW_TRANSMISSION_G'])
		gfaint_n,gfaint_s = 24.2,24.1
	else:
		sys.exit('wrong gkey; exiting')
	gr   = -2.5*np.log10(d['FLUX_G']/d['MW_TRANSMISSION_G']/d['FLUX_R']*d['MW_TRANSMISSION_R'])
	rz   = -2.5*np.log10(d['FLUX_R']/d['MW_TRANSMISSION_R']/d['FLUX_Z']*d['MW_TRANSMISSION_Z'])
	#
	mydict= {key:np.zeros(len(d),dtype=bool) for key in ['fdr','faint','blue','red']}
	for cap,lowzcut_zp,gfaint in zip(
		[(d['PHOTSYS']=='N'),(d['PHOTSYS']=='S')],
		[-0.35,-0.15],
		[gfaint_n,gfaint_s]):
		coii             = (gr+1.2*rz<1.6-7.2*(g-gfaint))
		fdrbox           = (rz>0.3) & (rz<1.6) & (gr<1.15*rz+lowzcut_zp) & (gr<-1.2*rz+1.6)
		mydict['fdr']   |= (cap) & (coii) & (fdrbox) & (g<gfaint)
		mydict['faint'] |= (cap) & (coii) & (fdrbox) & (g>gfaint)
		mydict['blue']  |= (cap) & (coii) & (rz<0.3)
		mydict['red']   |= (cap) & (coii) & (gr>-1.2*rz+1.6)
	return mydict


# initialising pixweight with ADM pixweight
h        = fits.open(sv1pixfn)
nside,nest = h[1].header['HPXNSIDE'],h[1].header['HPXNEST']
npix     = hp.nside2npix(nside)
pixarea  = hp.nside2pixarea(nside,degrees=True)
collist  = []
keys = ['HPXPIXEL', 'FRACAREA', 'STARDENS', 'EBV',
		'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z',
		'GALDEPTH_G', 'GALDEPTH_R', 'GALDEPTH_Z', 'PSFDEPTH_W1', 'PSFDEPTH_W2',
		'PSFSIZE_G', 'PSFSIZE_R', 'PSFSIZE_Z',
		'FRACAREA_13312', 'FRACAREA_13314', 'FRACAREA_14562',
		'ELG', 'ELG_SV_GTOT', 'ELG_SV_GFIB', 'ELG_FDR_GTOT', 'ELG_FDR_GFIB']
for key in keys:
	collist += [h[1].columns[key]]
log.info('{:.1f}s\tadm pixweight read'.format(time()-start))
# fracarea_elg (bright,galaxy,cluster + nobs_grz>0)
dens     = 0.
fns      = np.sort(glob(randdir+'randoms-1-*.fits'))
ras,decs = [],[]
for ifn,fn in enumerate(fns):
	dens += float(fits.getheader(fn,1)['DENSITY'])
	d     = fitsio.read(fn,columns=['RA','DEC','NOBS_G','NOBS_R','NOBS_Z','MASKBITS'])
	keep  = (d['NOBS_G']>0) & (d['NOBS_R']>0) & (d['NOBS_Z']>0)
	for b in [1,12,13]:
		keep &= ((d['MASKBITS'] & 2**b)==0)
	log.info('{:.1f}s\t{}\t({:.0f}/{:.0f})\t-> {:.0f}\trandoms'.format(time()-start,fn.split('/')[-1],ifn,len(fns)-1,keep.sum()))
	ras  += d['RA'] [keep].tolist()
	decs += d['DEC'][keep].tolist()
ras,decs    = np.array(ras),np.array(decs)
pixs        = hp.ang2pix(nside,(90.-decs)*np.pi/180.,ras*np.pi/180.,nest=nest)
i,c         = np.unique(pixs,return_counts=True)
fracarea    = np.zeros(npix)
fracarea[i] = c / dens / pixarea
collist    += [fits.Column(name='FRACAREA_ELG',format='E',array=fracarea)]
log.info('{:.1f}s\tfracarea_elg done'.format(time()-start))
# north,south,des
theta,phi= hp.pix2ang(nside,np.arange(npix),nest=nest)
ras,decs = 180./np.pi*phi,90.-180./np.pi*theta
c        = SkyCoord(ras*units.degree,decs*units.degree, frame='icrs')
isnorth  = (fracarea>0) & (c.galactic.b.value>0) & (decs>32.375)
issouth  = (fracarea>0) & ((c.galactic.b.value<0) | ((c.galactic.b.value>0) & (decs<32.375)))
isdes    = (fracarea>0) & (get_isdes(ras,decs))
collist += [fits.Column(name='RA',     format='E',array=ras)]
collist += [fits.Column(name='DEC',    format='E',array=decs)]
collist += [fits.Column(name='ISNORTH',format='L',array=isnorth)]
collist += [fits.Column(name='ISSOUTH',format='L',array=issouth)]
collist += [fits.Column(name='ISDES',  format='L',array=isdes)]
# looping through targets
dens   = {key:np.zeros(npix) for key in ['gtot','gfib']}
keys   = ['blue','red','fdr','faint']
for gkey in ['gtot','gfib']:
	dens[gkey] = {key:np.zeros(npix) for key in keys}
fns    = np.sort(glob(sv1targdir+'sv1targets-dark-hp-*.fits'))
for ifn,fn in enumerate(fns):
	d  = fitsio.read(fn,columns=[
				'RA','DEC','PHOTSYS',
				'FLUX_G','FLUX_R','FLUX_Z','FIBERFLUX_G',
				'MW_TRANSMISSION_G','MW_TRANSMISSION_R','MW_TRANSMISSION_Z',
				'SV1_DESI_TARGET'])
	d  = d[(d['SV1_DESI_TARGET'] & 2**1)>0]
	log.info('{:.1f}s\t{}\t({:.0f}/{:.0f})\t-> {:.0f}\telg targets'.format(time()-start,fn.split('/')[-1],ifn,len(fns)-1,len(d)))
	for gkey in ['gtot','gfib']:
		subs = get_subs(d,gkey)
		for key in keys:
			dk                  = d[subs[key]]
			pixs                = hp.ang2pix(nside,(90.-dk['DEC'])*np.pi/180.,dk['RA']*np.pi/180.,nest=nest)
			i,c                 = np.unique(pixs,return_counts=True)
			tmp                 = (fracarea[i]>0)
			ii,cc               = i[tmp],c[tmp]
			dens[gkey][key][ii]+= cc / fracarea[ii] / pixarea 
for gkey in ['gtot','gfib']:
	dens[gkey]['sv'] = np.zeros(npix)
	for key in keys:
		dens[gkey]['sv'] += dens[gkey][key]
#
for gkey in ['gtot','gfib']:
	for key in keys+['sv']:
		collist += [fits.Column(name=key.upper()+'_'+gkey.upper(),format='E',array=dens[gkey][key])]
hdu  = fits.BinTableHDU.from_columns(fits.ColDefs(collist))
hdu.writeto(args.outdir+'/dr9m-sv1pixweight-dark-elg.fits',overwrite=True)
log.info('{:.1f}s\tmypixweight done'.format(time()-start))
