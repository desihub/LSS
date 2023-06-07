import os
import argparse
import numpy as np
import healpy as hp
import fitsio
from mockfactory import Catalog
import desimodel.footprint

def radec2thphi(ra,dec):
    return (-dec+90.)*np.pi/180.,ra*np.pi/180.

def cut_randoms(cat_fns, tls_fn):

    tls = fitsio.read(tls_fn)
    cat = Catalog.read(cat_fns)
    inds = desimodel.footprint.is_point_in_desi(tls,cat['RA'], cat['DEC'])
    cat_cut = cat[inds]

    return cat_cut

def map_count(cat, nside = 1024):

    theta,phi = radec2thphi(cat['RA'],cat['DEC'])
    pix = hp.ang2pix(nside,theta,phi,nest=True)
    pix_non_zero,pix_count = np.unique(pix,return_counts=True)
    map_c = np.zeros(12*nside**2)
    map_c[pix_non_zero] = pix_count

    return map_c



parser = argparse.ArgumentParser()
parser.add_argument('--tracer', help='tracer to be selected', type=str, default='QSO')
parser.add_argument('--basedir', help='where to find catalogs', type=str, default='/global/cfs/cdirs/desi/survey/catalogs')
parser.add_argument('--survey', help='e.g., SV3 or main', type=str, choices=['SV3', 'DA02', 'main','Y1'], default='Y1')
parser.add_argument('--verspec', help='version for redshifts', type=str, default='iron')
parser.add_argument('--version', help='catalog version', type=str, default='test')
parser.add_argument('--nside', help='nside for healpix map', type=int, default=1024)
parser.add_argument('--nran', help='number of random files to combine together (1-18 available)', type=int, default=18)
parser.add_argument('--outdir', help='base directory for output', type=str, default=None)

args = parser.parse_args()


target_fns = [f'/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-{idx}.fits' for idx in range(args.nran)]
fibered_fns = [f'//global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/{args.verspec}/LSScats/{args.version}/{args.tracer}_{idx}_full.ran.fits' for idx in range(args.nran)]

if args.tracer=='BGS_BRIGHT':
    tile_fn = f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/tiles-BRIGHT.fits'
else:
    tile_fn = f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/tiles-DARK.fits'

if args.outdir is None:
    outdir = f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/{args.verspec}/LSScats/{args.version}/'
else:
    outdir = args.outdir
print(outdir)
save_fn = os.path.join(outdir,f'healpix_map_ran_comp_{args.tracer}.fits')

target = cut_randoms(target_fns, tile_fn)
fibered = cut_randoms(fibered_fns, tile_fn)
nside = args.nside
map_fib = map_count(fibered, nside=nside)
map_target = map_count(target, nside=nside)

map_weight = np.divide(map_fib, map_target, out=np.zeros_like(map_fib), where=map_fib!=0)

hp.write_map(save_fn, map_weight, overwrite=True)
