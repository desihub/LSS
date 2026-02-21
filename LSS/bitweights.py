from desiutil.iers import freeze_iers
freeze_iers()

"""
Functions for calculating PIP weights
"""
import os
import h5py
import numpy as np
from astropy.table import Table
try:
    from fiberassign.targets import (Targets, TargetsAvailable, TargetTree,
                                 LocationsAvailable, load_target_table)
except:
    from fiberassign.targets import (Targets, TargetsAvailable,
                                LocationsAvailable, load_target_table)

from fiberassign.assign import Assignment


def get_targets(mtlfile, skyfile,faver='5.0.0'):
    """
    Load target and information
    """
    # Read mtl file
    mtl = Table.read(mtlfile)
    if 'SUBPRIORITY' not in mtl.dtype.names:
        mtl['SUBPRIORITY'] = np.ones(len(mtl))
    if 'OBSCONDITIONS' not in mtl.dtype.names:
        mtl['OBSCONDITIONS'] = np.ones(len(mtl), dtype=int)
    if 'DESI_TARGET' not in mtl.dtype.names:
        mtl['DESI_TARGET'] = np.ones(len(mtl), dtype=int)

    # Load targets
    tgs = Targets()
    # Load sky targets
    sky = None
    if skyfile:
        sky = Table.read(skyfile)

    mver = int(faver[:1])
    if mver < 5:
        load_target_table(tgs, mtl)
        if sky is not None:
            load_target_table(tgs, sky)
    else:
        from fiberassign.targets import TargetTagalong,create_tagalong 
        tagalong = create_tagalong()
        load_target_table(tgs,tagalong, mtl) 
        if sky is not None:
            load_target_table(tgs, tagalong,sky)
            
    return mtl, tgs, sky

def setup_fba(mtl, sky, tiles, hw,faver='5.0.0'):
    """
    Set up tiles, targets, etc. for fiber assignment
    """
    # Load targets and target tree
    tgs = Targets()
    mver = int(faver[:1])
    if mver < 5:
        load_target_table(tgs, mtl)
        if sky:
            load_target_table(tgs, sky)

    else:
        from fiberassign.targets import TargetTagalong,create_tagalong 
        tagalong = create_tagalong()
        load_target_table(tgs,tagalong, mtl)  
        if sky:
            load_target_table(tgs, tagalong,sky)
        
    if mver < 3:
        from fiberassign.targets import (TargetTree)
        tree = TargetTree(tgs)

    # Compute available targets / locations
    if faver == '2.3.0':
        tgsavail = TargetsAvailable(hw, tgs, tiles, tree)
        favail = LocationsAvailable(tgsavail)
        asgn = Assignment(tgs, tgsavail, favail)
    if faver == '2.4.0' or faver == '2.5.0' or faver == '2.5.1':
        tgsavail = TargetsAvailable(hw, tgs, tiles, tree)
        favail = LocationsAvailable(tgsavail)
        asgn = Assignment(tgs, tgsavail, favail,{}) #this is needed for fiberassign 2.4 and higher(?)
    if mver >= 3 and mver < 5:
        from fiberassign.targets import targets_in_tiles
        tile_targetids, tile_x, tile_y = targets_in_tiles(hw, tgs, tiles)
        tgsavail = TargetsAvailable(hw, tiles, tile_targetids, tile_x, tile_y)
        favail = LocationsAvailable(tgsavail)
        asgn = Assignment(tgs, tgsavail, favail,{}) #this is needed for fiberassign 2.4 and higher(?)
    if mver >= 5:
        from fiberassign.targets import targets_in_tiles
        tile_targetids, tile_x, tile_y, tile_xy_cs5 = targets_in_tiles(hw, tgs, tiles,tagalong)
        tgsavail = TargetsAvailable(hw, tiles, tile_targetids, tile_x, tile_y)
        favail = LocationsAvailable(tgsavail)
        asgn = Assignment(tgs, tgsavail, favail,{}) #this is needed for fiberassign 2.4 and higher(?) to fake stucksky

    return asgn

def update_bitweights(realization, asgn, tileids, tg_ids, tg_ids2idx, bitweights):
    """
    Update bit weights for assigned science targets
    """
    for tileid in tileids:
        try: # Find which targets were assigned
            adata = asgn.tile_location_target(tileid)
            for loc, tgid in adata.items():
                idx = tg_ids2idx[tgid]
                bitweights[realization * len(tg_ids) + idx] = True
        except:
            pass

    return bitweights

def pack_bitweights(array):
    """
    Creates an array of bitwise weights stored as 64-bit signed integers
    Input: a 2D boolean array of shape (Ngal, Nreal), where Ngal is the total number 
           of target galaxies, and Nreal is the number of fibre assignment realizations.
    Output: returns a 2D array of 64-bit signed integers. 
    """
    Nbits = 64
    dtype = np.int64
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

def unpack_bitweights(we):
    Nwe= 1
    Nbits = 64
    Ngal = np.shape(we)[0]
    Nreal = Nbits*Nwe
    print('Nbits, Nwe = ',Nbits,Nwe)
    print('Nreal = ',Nreal)
    print('Ngal = ',Ngal)
    true8=[np.uint8(255) for n in range(0, Ngal)]
    array_bool = np.zeros((Ngal,Nreal), dtype=bool)
    for j in range(Nwe):
        lg = np.zeros((Ngal, Nbits), dtype=bool)
        for i in range(Nbits//8):
            chunk8 = np.uint8(np.bitwise_and(np.right_shift(we,8*i), true8))
            lg[:,Nbits-8*(i+1):Nbits-i*8] = np.reshape(np.unpackbits(chunk8), (Ngal, 8))
        array_bool[:,j*Nbits:(j+1)*Nbits] = lg[:,::-1]
    return array_bool

def write_output(outdir, outfilename, overwrite, fileformat, targets, bitvectors, desi_target_key=None):
    """
    Write output file containing bit weights
    """
    try:
        os.makedirs(outdir)
    except:
        pass

    # Output fits files
    if fileformat == 'fits':
        outfile = os.path.join(outdir, outfilename)
        output = Table()
        output['TARGETID'] = targets['TARGETID']
        if desi_target_key:
            output['{}'.format(desi_target_key)] = targets['{}'.format(desi_target_key)]
        output['RA'] = targets['RA']
        output['DEC'] = targets['DEC']
        output['Z'] = targets['RSDZ']
        for i in range(bitvectors.shape[1]):
            output['BITWEIGHT{}'.format(i)] = bitvectors[:,i]
        output.write(outfile, overwrite=overwrite)

    # Output hdf5 files
    elif fileformat == 'hdf5':
        outfile = os.path.join(outdir, outfilename)
        outfile = h5py.File(outfile, 'w')
        outfile.create_dataset('TARGETID', data=targets['TARGETID'])
        if desi_target_key:
            outfile.create_dataset('{}'.format(desi_target_key), data=targets['{}'.format(desi_target_key)])
        outfile.create_dataset('RA', data=targets['RA'])
        outfile.create_dataset('DEC', data=targets['DEC'])
        outfile.create_dataset('Z', data=targets['Z'])
        for i in range(bitvectors.shape[1]):
            outfile.create_dataset('BITWEIGHT{}'.format(i), data=bitvectors[:,i])
        outfile.close()

    else:
        raise TypeError("Must specify either fits or hdf5 as output file format.")

    return

