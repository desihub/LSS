import os
import argparse
import numpy as np
import healpy as hp
import fitsio
from mockfactory import Catalog
import desimodel.footprint

def radec2thphi(ra,dec):
    return (-dec+90.)*np.pi/180.,ra*np.pi/180.


def map_ntile(cat, nside = 1024):

    theta,phi = radec2thphi(cat['RA'],cat['DEC'])
    pix = hp.ang2pix(nside,theta,phi,nest=True)
    map_c = np.zeros(12*nside**2)
    map_ntile = np.zeros(12*nside**2)
    for ii in range(0,len(pix)):
        pi = pix[ii]
        ntile = cat['NTILE'][ii]
        map_c[pi] += 1.
        map_ntile[pi] += ntile
    sel = map_c > 0
    map_ntile[sel] /= map_c[sel]
    print(np.min(map_ntile),np.max(map_ntile))
    return map_ntile



parser = argparse.ArgumentParser()
parser.add_argument('--prog', help='tracer to be selected', type=str, choices=['dark', 'bright'],default='dark')
parser.add_argument('--survey', help='e.g., SV3 or main', type=str, choices=['SV3', 'DA02', 'main','Y1'], default='Y1')
parser.add_argument('--verspec', help='version for redshifts', type=str, default='iron')
parser.add_argument('--version', help='catalog version', type=str, default='test')
parser.add_argument('--nside', help='nside for healpix map', type=int, default=1024)
parser.add_argument('--nran', help='number of random files to combine together (1-18 available)', type=int, default=18)
args = parser.parse_args()


indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/'+args.verspec+'/LSScats/'+args.version+'/'


args = parser.parse_args()

print('test')

tracer = 'QSO'
if args.prog == 'bright':
    tracer = 'BGS_BRIGHT'

if args.nran == 1:
    ran_fn = f'{indir}/{tracer}_1_full.ran.fits'
    ran_cat = fitsio.read(ran_fn)
else:
    ran_fns = [f'{indir}/{tracer}_{idx}_full.ran.fits' for idx in range(args.nran)]
    ran_cat = Catalog.read(ran_fns)


save_fn = indir+f'/healpix_map_ntile_{args.prog}.fits'

nside = args.nside
map_nt = map_ntile(ran_cat, nside=nside)

hp.write_map(save_fn, map_nt, overwrite=True)
