"""
Functions for calculating PIP weights
"""
import os
import h5py
import numpy as np
from astropy.table import Table
from fiberassign.tiles import load_tiles
from fiberassign.targets import (Targets, TargetsAvailable, TargetTree,
                                 LocationsAvailable, load_target_table)
from fiberassign.assign import Assignment


def get_tiles(tilefile, footprint):
    """
    Load tile information
    """
    # Read tiles we are using
    tileselect = None
    if tilefile is not None:
        tileselect = list()
        with open(tilefile, "r") as f:
            for line in f:
                # Try to convert the first column to an integer
                try:
                    tileselect.append(int(line.split()[0]))
                except ValueError:
                    pass
    tiles = load_tiles(tiles_file=footprint,
                       select=tileselect)

    return tiles

def get_targets(mtlfile, skyfile):
    """
    Load target and randoms information
    """
    # Read mtl file
    mtl = Table.read(mtlfile)
    if 'SUBPRIORITY' not in mtl:
        mtl['SUBPRIORITY'] = np.ones(len(mtl))
    if 'OBSCONDITIONS' not in mtl:
        mtl['OBSCONDITIONS'] = np.ones(len(mtl), dtype=int)
    if 'DESI_TARGET' not in mtl:
        mtl['DESI_TARGET'] = np.ones(len(mtl), dtype=int)

    # Load science targets
    tgs = Targets()
    load_target_table(tgs, mtl)

    # Load sky targets
    sky = None
    if skyfile:
        sky = Table.read(skyfile)
        load_target_table(tgs, sky)

    return mtl, tgs, sky

def setup_fba(mtl, sky, tiles, hw):
    """
    Set up tiles, targets, etc. for fiber assignment
    """
    # Load targets and target tree
    tgs = Targets()
    load_target_table(tgs, mtl)
    if sky:
        load_target_table(tgs, sky)
    tree = TargetTree(tgs)

    # Compute available targets / locations
    tgsavail = TargetsAvailable(hw, tgs, tiles, tree)
    favail = LocationsAvailable(tgsavail)
    asgn = Assignment(tgs, tgsavail, favail)

    return asgn

def update_bitweights(realization, asgn, tiles, tg_science, bitweights):
    """
    Update bit weights for assigned science targets
    """
    assigned = []
    for tile_id in tiles.id:
        try: # Find which targets were assigned
            adata = asgn.tile_location_target(tile_id)
            for loc, tgid in adata.items():
                assigned.append(tgid)
        except:
            pass

    # Update bit weights
    idas = np.isin(tg_science, assigned)
    min_targ = realization * len(tg_science)
    max_targ = (realization + 1) * len(tg_science)
    bitweights[min_targ:max_targ] = idas

    return idas, bitweights

def pack_bitweights(array):
    """
    Creates an array of bitwise weights stored as 64-bit signed integers
    Input: a 2D boolean array of shape (Ngal, Nreal), where Ngal is the total number 
           of target galaxies, and Nreal is the number of fibre assignment realizations.
    Output: returns a 2D array of 64-bit signed integers. 
    """
    Nbits=64
    dtype=np.int64
    Ngal, Nreal = array.shape           # total number of realizations and number of target galaxies
    Nout = (Nreal + Nbits - 1) // Nbits # number of output columns
    # intermediate arrays
    bitw8 = np.zeros((Ngal, 8), dtype="i")   # array of individual bits of 8 realizations
    bitweights = np.zeros(Ngal, dtype=dtype) # array of 64-bit integers
    # array to store final output
    output_array = np.zeros((Ngal, Nout), dtype=dtype)
    idx_out = 0 # initial column in output_array
    # loop through realizations to build bitwise weights
    for i in range(Nreal):
        bitw8[array[:,i], i%8] = 1
        arr = np.array(np.packbits(bitw8[:,::-1]), dtype=dtype)
        bitweights = np.bitwise_or(bitweights, np.left_shift(arr, 8*((i%Nbits)//8)))
        if (i+1)%Nbits == 0 or i+1 == Nreal:
            output_array[:,idx_out] = bitweights
            bitweights[:] = 0
            idx_out += 1
        if (i+1)%8 == 0:
            bitw8[:] = 0
    return output_array

def write_output(outtype, outdir, fileformat, targets, idas, bitvectors):
    """
    Write output file containing bit weights
    """
    try:
        os.makedirs(outdir)
    except:
        pass

    if fileformat == 'fits':
        outfile = os.path.join(outdir, '{}.fits'.format(outtype))
        output = Table()
        if outtype == 'targeted' or outtype == 'parent':
            output['TARGETID'] = targets['TARGETID']
        output['RA'] = targets['RA']
        output['DEC'] = targets['DEC']
        output['Z'] = targets['Z']
        for i in range(bitvectors.shape[1]):
            if outtype == 'targeted':
                output['BITWEIGHT{}'.format(i)] = bitvectors[:,i][idas]
            elif outtype == 'parent':
                output['BITWEIGHT{}'.format(i)] = -np.ones(len(bitvectors[:,i]), dtype=int)
            elif outtype == 'randoms':
                output['BITWEIGHT{}'.format(i)] = -np.ones(len(targets), dtype=int)
        output.write(outfile)

    elif fileformat == 'hdf5':
        outfile = os.path.join(outdir, '{}.hdf5'.format(outtype))
        outfile = h5py.File(outfile, 'w')
        if outtype == 'targeted' or outtype == 'parent':
            outfile.create_dataset('TARGETID', data=targets['TARGETID'])
        outfile.create_dataset('RA', data=targets['RA'])
        outfile.create_dataset('DEC', data=targets['DEC'])
        outfile.create_dataset('Z', data=targets['Z'])
        for i in range(bitvectors.shape[1]):
            if outtype == 'targeted':
                outfile.create_dataset('BITWEIGHT{}'.format(i), data=bitvectors[:,i][idas])
            elif outtype == 'parent':
                outfile.create_dataset('BITWEIGHT{}'.format(i), data=-np.ones(len(bitvectors[:,i]), dtype=int))
            elif outtype == 'randoms':
                outfile.create_dataset('BITWEIGHT{}'.format(i), data=-np.ones(len(targets), dtype=int))
        outfile.close()

    else:
        raise TypeError("Must specify either fits or hdf5 as output file format.")

    return

