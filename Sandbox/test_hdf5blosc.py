#!/usr/bin/env python

"""
Experiment with reformatting a random file

#- Create environment with hdf5plugin
module load python
conda create -n hdf5blosc fitsio h5py hdf5plugin ipython astropy
conda activate hdf5blosc
"""

import os, sys, glob
import numpy as np
from astropy.table import Table

import fitsio
import numpy as np
import h5py
import hdf5plugin

import argparse
p = argparse.ArgumentParser()
p.add_argument('-i', '--infile', required=True, help="input altmtl random FITS filename")
p.add_argument('-o', '--outdir', default=os.getcwd(), help="output directory (default %(default)s)")
# p.add_argument('-n', '--night', type=int, help="night to process")
# p.add_argument('--debug', action="store_true", help="...")
args = p.parse_args()

# for easier debugging cutting and pasting 
infile = args.infile
outdir = args.outdir

"""
# cd /global/cfs/cdirs/desi/users/sjbailey/dev/altmtl_random_format
infile = 'LRG_NGC_1_clustering.ran.fits'
outdir = os.getcwd()
"""

assert infile.endswith('.fits')

r1 = Table.read(infile)

def print_column_types(table):
    for col in table.colnames:
        print(f'| {col} | {table[col].dtype} |')

def reduce_column_precision(table):
    """
    Returns output table with same columns but reduced precision optimized for BLOSC compression
    """
    result = Table()
    result.meta.update(table.meta)
    for col in table.colnames:
        if col in ('RA', 'DEC'):
            #- float64 -> float32 -> float64
            result[col] = table[col].astype(np.float32).astype(np.float64)
        elif col in ('FRAC_TLOBS_TILES', 'NX') or col.startswith('WEIGHT'):
            #- float64 -> float16 -> float64
            result[col] = table[col].astype(np.float16).astype(np.float64)
        elif col in ('NTILE',):
            #- intXX -> int8
            result[col] = table[col].astype(np.int8)
        else:
            #- default unchanged type
            result[col] = table[col]

    assert result.colnames == table.colnames
    return result

#- Optimize column dtypes
r2 = reduce_column_precision(r1)

print('Input table column types')
print_column_types(r1)

print('\nOptimized column types (some were cast via lower precision then restored)')
print_column_types(r2)

#- Write as FITS
outfile = os.path.basename(infile).replace('.fits', '.reformat.fits')
outfile = f'{outdir}/{outfile}'
r2.write(outfile, overwrite=True)

#- Write both tables as blosc-compressed hdf5
def write_hdf5_blosc(filename, table):
    """Write table to filename using hdf5 blosc compression; code adapted from Joe DeRose"""
    if os.path.exists(filename):
        print(f'Replacing {filename}')
        os.remove(filename)

    for k in table.dtype.names:
        data = table[k]
        dt = table.dtype[k]
        if dt == '<U1':
            dt = 'S1'  
            data = np.array(data, dtype=dt)
            
        with h5py.File(filename, 'a') as fn:
            # Using Blosc with default settings
            fn.create_dataset(k, data=data, dtype=dt,
                                compression=hdf5plugin.Blosc(cname='zstd', clevel=5))

outfile = os.path.join(outdir, os.path.basename(infile).replace('.fits', '.h5'))
write_hdf5_blosc(outfile, r1)

outfile = os.path.join(outdir, os.path.basename(infile).replace('.fits', '.reformat.h5'))
write_hdf5_blosc(outfile, r2)

