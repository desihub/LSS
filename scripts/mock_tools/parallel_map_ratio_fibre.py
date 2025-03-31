#In order to run this script, you need the cosmodesi environment
#source /global/common/software/desi/desi_environment.sh main
#source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
#Then, request two nodes interactively (you need 2 for memory, the code works with ~400GB of data)
#salloc -N 2 --qos interactive --time 04:00:00 --constraint cpu --account desi 
#then, run
# srun -n 18 python /global/cfs/cdirs/desi/users/schiarenza/LSS/scripts/mock_tools/parallel_map_ratio_fibre.py 

import os
import argparse
import numpy as np
import healpy as hp
import fitsio
#from mockfactory import Catalog -> This causes issues with parallelization and is not needed.
import desimodel.footprint
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def radec2thphi(ra, dec):
    return (-dec + 90.) * np.pi / 180., ra * np.pi / 180.

def cut_randoms(cat_fn, tls_fn):
    print(f"Rank {rank}: Reading catalog {cat_fn}")  # Debugging
    tls = fitsio.read(tls_fn)
    cat = fitsio.read(cat_fn)
    inds = desimodel.footprint.is_point_in_desi(tls, cat['RA'], cat['DEC'])
    cat_cut = cat[inds]
    print(f"Rank {rank}: Finished cutting randoms")  # Debugging
    return cat_cut

def map_count(cat, nside=1024):
    print(f"Rank {rank}: Creating Healpix map")  # Debugging
    theta, phi = radec2thphi(cat['RA'], cat['DEC'])
    pix = hp.ang2pix(nside, theta, phi, nest=True)
    pix_non_zero, pix_count = np.unique(pix, return_counts=True)
    map_c = np.zeros(12 * nside ** 2)
    map_c[pix_non_zero] = pix_count
    print(f"Rank {rank}: Healpix map created")  # Debugging
    return map_c

# Argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('--tracer', type=str, default='QSO')
parser.add_argument('--basedir', type=str, default='/global/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument('--survey', type=str, choices=['SV3', 'DA02', 'main', 'Y1'], default='Y1')
parser.add_argument('--verspec', type=str, default='iron')
parser.add_argument('--version', type=str, default='v1.5')
parser.add_argument('--nside', type=int, default=1024)
parser.add_argument('--nran', type=int, default=18)
parser.add_argument('--outdir', type=str, default=None)

args = parser.parse_args()

if rank == 0:
    print(f"Starting MPI with {size} processes")

# Ensure number of MPI ranks does not exceed the number of random files
if size > args.nran:
    if rank == 0:
        print(f"Warning: More MPI ranks ({size}) than random files ({args.nran}). Some ranks will be idle.")
    size = args.nran  # Limit MPI size to number of files

# Assign idx based on rank
if rank < args.nran:  
    idx = rank  # Each rank gets a unique index
else:
    idx = None  # Extra ranks do nothing

if args.tracer == 'BGS_BRIGHT':
    tile_fn = f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/tiles-BRIGHT.fits'
else:
    tile_fn = f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/tiles-DARK.fits'

if args.outdir is None:
    outdir = f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/{args.verspec}/LSScats/{args.version}/'
else:
    outdir = args.outdir

save_fn = os.path.join(outdir, f'healpix_map_ran_comp_{args.tracer}_nside_{args.nside}.fits')

if idx is not None:
    print(f"Rank {rank}: Processing file index {idx}")  # Debugging
    target_fn = f'/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-{idx}.fits'
    fibered_fn = f'/global/cfs/cdirs/desi/survey/catalogs/{args.survey}/LSS/{args.verspec}/LSScats/{args.version}/{args.tracer}_{idx}_full.ran.fits'

    target = cut_randoms(target_fn, tile_fn)
    fibered = cut_randoms(fibered_fn, tile_fn)
    nside = args.nside
    map_fib = map_count(fibered, nside=nside)
    map_target = map_count(target, nside=nside)
    print(f"Rank {rank}: Finished processing index {idx}")  # Debugging
else:
    map_fib = None
    map_target = None

# Gather results at rank 0
if rank == 0:
    print("Rank 0: Gathering results from all ranks")  # Debugging
map_fib_all = comm.gather(map_fib, root=0)
map_target_all = comm.gather(map_target, root=0)

if rank == 0:
    print("Rank 0: Received all results. Summing maps.")  # Debugging
    map_fib_all = [m for m in map_fib_all if m is not None]
    map_target_all = [m for m in map_target_all if m is not None]

    # Sum all individual maps
    map_fib_final = np.sum(map_fib_all, axis=0)
    map_target_final = np.sum(map_target_all, axis=0)
    
    print("Rank 0: Computing final weight map")  # Debugging
    map_weight = np.divide(map_fib_final, map_target_final, out=np.zeros_like(map_fib_final), where=map_fib_final != 0)

    print(f"Rank 0: Saving final healpix map to {save_fn}")  # Debugging
    hp.write_map(save_fn, map_weight, overwrite=True)

if rank == 0:
    print("Processing complete. Exiting.")  # Final message