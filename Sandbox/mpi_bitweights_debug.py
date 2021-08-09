#!/usr/bin/env python
"""
Compute bit weights in parallel to be used for PIP or IIP weights

Example run at NERSC:

srun -N 16 -c 8 -C haswell -A desi --qos=interactive -t 0:30:0 mpi_bitweights
     --mtl LRG_oneperztrue_clus.dat.fits --format hdf5 --outdir output --realizations 128

Notes:
- Runs on master version of all DESI repos except for fiberassign (use 2.3.0 for now)
- Must add path_to_LSS_repo/LSS/bin to $PATH and path_to_LSS/LSS/py to $PYTHONPATH
- bit weights are stored as 64 bit integers, so the number of realizations must
  be a multiple of 64 (i.e. 64, 128, 192, etc. realizations correspond to output
  arrays BITWEIGHT0, BITWEIGHT1, BITWEIGHT2, etc.)
- n nodes (-N) x n cpus (-c) must be a factor of the number of realizations

Required Arguments:
  mtl     (str) : target fits file containing TARGETID, RA, DEC, Z

Optional Arguments:
  sky          (str) : fits file containing sky targets
  tiles        (str) : fits file containing tile information
  format       (str) : output file format (fits or hdf5)
  outdir       (str) : output directory
  realizations (int) : number of realizations

Outputs:
  targeted.fits : fits (hdf5) file containing targeted science objects
"""
import argparse

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mtl", type=str, required=True,
                        help="fits file with mtl, should contain only targets available to fibers")
    parser.add_argument("--sky", type=str, required=False, help="fits file with sky targets")
    parser.add_argument("--tiles", type=str, required=False, default=None,
                        help="Optional fits file containing tile information to use in the footprint,"
                        "default uses desimodel footprint.")
    parser.add_argument("--format", type=str, required=False, default="fits",
                        help="File format for outputs (either fits or hdf5)")
    parser.add_argument("--outdir", type=str, required=False, default=".", help="Output directory")
    parser.add_argument("--realizations", type=int, required=False, default=1, help="Number of realizations")
    parser.add_argument("--survey", type=str, required=False, default="sky", help="e.g., sv3 or main")

    args = parser.parse_args()
    return args

def main(args):
    print('started')
    # Initialize mpi
    #from mpi4py import MPI
    #comm = MPI.COMM_WORLD
    #mpi_procs = comm.size
    #mpi_rank = comm.rank

    import numpy as np
    from astropy.table import Table, vstack
    from LSS.bitweights import (get_targets, setup_fba, update_bitweights,
                                pack_bitweights, write_output)
    from fiberassign.hardware import load_hardware
    from fiberassign.tiles import load_tiles
    from fiberassign._internal import Tiles
    from fiberassign.assign import run
    from fiberassign.targets import (Targets, load_target_table)
    
    from matplotlib import pyplot as plt

    # Load tile data
    tiledata = None
    if args.tiles:
        #tiledata = Table.read(args.tiles)
         #tiles = load_tiles(args.tiles)
         tiles = load_tiles(tiles_file=args.tiles, select=None)#,obstime=args.obsdate, obstheta=args.fieldrot, obsha=args.ha)
    else:
         print('YOU SHOULD REALLY USE A TILES FILE (unless this is for the original e2e mock)!')
         tiles = load_tiles()

    # Get target and information
    mtl, tgs, sky = get_targets(args.mtl, args.sky)
    
    tg_ids = tgs.ids()
    tg_ids2idx = {y: x for x, y in enumerate(tg_ids)}
    n_target = len(tg_ids)
    print('there are targets '+str(n_target))

    # Divide realizations among processes
    n_realization = args.realizations
    realizations = np.arange(n_realization, dtype=np.int32)
    #my_realizations = np.array_split(realizations, mpi_procs)[mpi_rank]

    # Bit weight array for all targets and realizations
    bitweights = np.zeros(n_realization * n_target, dtype=bool)

    # Load hardward for fiber assignment
    hw = load_hardware(rundate='2021-04-06T00:39:37') #rundate for first SV3 tiles

    for realization in realizations:#my_realizations:
        # Set the seed based on the realization, so that the result is reproducible
        # regardless of which process is working on the realization
        np.random.seed(realization)

        # Randomize subpriorities
        mtl['SUBPRIORITY'] = np.random.random_sample(size=len(mtl))

        # Loop through tiles individually to load appropriate tile information
        if tiledata:
            for tdata in tiledata:
                # Load tile
                tileid = [tdata['TILEID']]
                tilera = [tdata['RA']]
                tileha = [tdata['FA_HA']]
                tiledec = [tdata['DEC']]
                rundate = [tdata['RUNDATE']]
                obscond = tdata['OBSCONDITIONS']
                theta = [tdata['FIELDROT']]
                tile = Tiles(tileid, tilera, tiledec, obscond, rundate, theta, tileha)
    
                # Load appropriate hardware state
                hw = load_hardware(rundate=rundate[0])
    
                # Set up and run fiber assignment
                asgn = setup_fba(mtl, sky, tile, hw)
                run(asgn)
    
                # Update bit weights for assigned science targets
                bitweights = update_bitweights(realization, asgn, tileid, tg_ids, tg_ids2idx, bitweights)

        # Run all tiles simultaneously if we don't have tile information
        else:
            # Set up and run fiber assignment
            asgn = setup_fba(mtl, sky, tiles, hw)
            run(asgn)

            # Update bit weights for assigned science targets
            bitweights = update_bitweights(realization, asgn, tiles.id, tg_ids, tg_ids2idx, bitweights)

    # Gather weights from all processes
    #gather_weights = None
    #if mpi_rank == 0:
    #gather_weights = np.empty(len(bitweights), dtype=bool)
    #min_targ = np.min(my_realizations) * n_target
    #max_targ = (np.max(my_realizations) + 1) * n_target
    #comm.Gather(bitweights[min_targ:max_targ], gather_weights, root=0)

    print(len(bitweights),sum(bitweights))
    #if mpi_rank == 0:
    s = 0
    for tileid in tiles.id:
        adata = asgn.tile_location_target(tileid)
        for loc, tgid in adata.items():
            if s == 0:
                tids = tgid
            else:
                tids = np.concatenate([tids,tgid])
    tidsu = np.unique(tids)
    print(len(tidsu))
    w = np.isin(mtl['TARGETID'],tidsu)

    #plt.plot(mtl[bitweights]['RA'],mtl[bitweights]['DEC'],',k')
    plt.plot(mtl[w]['RA'],mtl[w]['DEC'],',k')
    plt.xlim(178,188)
    plt.ylim(-5,5)
    plt.show()
    if n_realization == False:
        # Pack weights into 64 bit integers
        #weights = np.array_split(gather_weights, n_realization)
        #bweights = np.array(weights).T
        #bitvectors = pack_bitweights(bweights)
        bitvectors = pack_bitweights(bitweights)

        # Figure out targeted sample
        idas = np.zeros(n_target, dtype=bool)
        for real in range(n_realization):
            was_assigned = weights[real] == True
            idas[was_assigned] = True

        print('number of targets that got assigned '+str(np.sum(idas)))
        # Find which target bits to output, if SV3_DESI_TARGET is there, use that, DESI_TARGET otherwise
        desi_target_key = None
        for key in mtl.keys():
            if 'SV3_DESI_TARGET' == key:
                desi_target_key = key
        if desi_target_key == None:
            for key in mtl.keys():
                if 'DESI_TARGET' == key:
                    desi_target_key = key
        

        # Combine sky and mtl
        if sky:
            outmtl = Table()
            outmtl['TARGETID'] = mtl['TARGETID']
            if desi_target_key:
                outmtl['{}'.format(desi_target_key)] = mtl['{}'.format(desi_target_key)]
            outmtl['RA'] = mtl['RA']
            outmtl['DEC'] = mtl['DEC']
            outmtl['Z'] = mtl['Z']
            outsky = Table()
            outsky['TARGETID'] = sky['TARGETID']
            if desi_target_key:
                outsky['{}'.format(desi_target_key)] = sky['{}'.format(desi_target_key)]
            outsky['RA'] = sky['RA']
            outsky['DEC'] = sky['DEC']
            outsky['Z'] = np.empty(len(sky))
            targets = vstack([outmtl, outsky])
        else:
            targets = mtl

        # Write output files
        write_output(args.outdir, args.format, targets[idas], bitvectors[idas], desi_target_key=desi_target_key)

if __name__ == "__main__":
    args = parse()
    print('started')
    print(args)
    main(args)
