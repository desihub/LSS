import os
import logging
import datetime
import argparse
import re

import numpy as np
import yaml

import fitsio
from astropy.table import Table,vstack

import pyrecon
from pyrecon import  utils,IterativeFFTParticleReconstruction,MultiGridReconstruction,setup_logging

from LSS.tabulated_cosmo import TabulatedDESI
cosmo = TabulatedDESI()
comoving_distance = cosmo.comoving_radial_distance

parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='everest')
parser.add_argument("--nthreads",default=1)
parser.add_argument("--rectype",help="IFT or MG supported so far",default='MG')
parser.add_argument("--convention",help="recsym or reciso supported so far",default='reciso')
parser.add_argument("--nran",help="how many of the random files to concatenate",default=5)


args = parser.parse_args()
print(args)

basedir = args.basedir
version = args.version
specrel = args.verspec

nran = int(args.nran)


maindir = basedir +'/main/LSS/'

ldirspec = maindir+specrel+'/'
    
dirout = ldirspec+'LSScats/'+version+'/'


regl = ['_DN','_N']
position_columns = ['RA','DEC','Z']
zmin = 0.3
zmax = 1.2
bias = 1.8
beta = 0.4
ff = beta*bias

setup_logging()

if args.rectype == 'MG':
    recfunc = MultiGridReconstruction
if args.rectype == 'IFT':
    recfunc = IterativeFFTParticleReconstruction

def getrdz_fromxyz(cat):
    distance, ra, dec = utils.cartesian_to_sky(cat)
    # DESI tabulated values go to z = 10. Do we want to go further?
    distance_to_redshift = utils.DistanceToRedshift(comoving_distance, zmax=10)
    z = distance_to_redshift(distance)
    return ra,dec,z 
    
for reg in regl:
    fb = dirout+'LRGzdone'+reg

    
    fcd = fb+'_clustering.dat.fits'    
    dat_cat = fitsio.read(fcd)
    seld = dat_cat['Z'] > zmin
    seld &= dat_cat['Z'] < zmax
    dat_cat = dat_cat[seld]
    print('using '+str(len(dat_cat))+' data entries')
    
    randoms_fn = [fb+ '_{:d}_clustering.ran.fits'.format( iran) for iran in range(nran)]
    ran_cat = vstack([Table.read(fn) for fn in randoms_fn])
    selr = ran_cat['Z'] > zmin
    selr &= ran_cat['Z'] < zmax
    ran_cat = ran_cat[selr]
    print('using '+str(len(ran_cat))+' random entries')

    dat_dis = comoving_distance(dat_cat[position_columns[2]])
    pos_dat = utils.sky_to_cartesian(dat_dis,dat_cat[position_columns[0]],dat_cat[position_columns[1]])
    ran_dis = comoving_distance(ran_cat[position_columns[2]])
    pos_ran = utils.sky_to_cartesian(ran_dis,ran_cat[position_columns[0]],ran_cat[position_columns[1]])
    recon = recfunc(f=ff, bias=bias, cellsize=7, los='local', positions=pos_ran, nthreads=int(args.nthreads), fft_engine='fftw', fft_plan='estimate')
    print('grid set up',flush=True)
    recon.assign_data(pos_dat,dat_cat['WEIGHT'])
    print('data assigned',flush=True)
    recon.assign_randoms(pos_ran,ran_cat['WEIGHT'])
    print('randoms assigned',flush=True)
    recon.set_density_contrast()
    print('density constrast calculated, now doing recon',flush=True)
    recon.run()
    print('recon has been run',flush=True)
    
    positions_rec = {}
    if args.rectype == 'IFT':
        positions_rec['data'] = pos_dat - recon.read_shifts('data')#, field='disp+rsd')
    else:
        positions_rec['data'] = pos_dat - recon.read_shifts(pos_dat)#, field='disp+rsd')
    
    positions_rec['randoms'] = pos_ran - recon.read_shifts(pos_ran, field='disp+rsd' if args.convention == 'recsym' else 'disp')
    
    fcro = fb+'_clustering_'+args.rectype+args.convention+'.ran.fits'
    fcdo = fb+'_clustering_'+args.rectype+args.convention+'.dat.fits'
    
    datt = Table(dat_cat)
    ra,dec,z = getrdz_fromxyz(positions_rec['data'])
    datt['RA'] = ra
    datt['DEC'] = dec
    datt['Z'] = z
    datt.write(fcdo,format='fits',overwrite=True)
    print('wrote data to '+fcdo)
    
    rant = Table(ran_cat)
    ra,dec,z = getrdz_fromxyz(positions_rec['randoms'])
    rant['RA'] = ra
    rant['DEC'] = dec
    rant['Z'] = z
    rant.write(fcro,format='fits',overwrite=True)
    print('wrote data to '+fcro)
    del datt
    del rant