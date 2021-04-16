#!/usr/bin/env python

"""
Compute PIP weights in parallel
"""

from mpi4py import MPI
import os, sys
import h5py
import argparse
import numpy as np
from astropy.table import Table
from fiberassign.hardware import load_hardware
from fiberassign.tiles import load_tiles
from fiberassign.targets import (Targets,
                                 TargetsAvailable,
                                 TargetTree,
                                 LocationsAvailable,
                                 load_target_table)
from fiberassign.assign import Assignment, run


def main():
    comm = MPI.COMM_WORLD
    mpi_procs = comm.size
    mpi_rank = comm.rank

    parser = argparse.ArgumentParser()
    parser.add_argument("--survey_log", type=str, required=False,
                        help="Eventually we would pass in a file containing the log"
                        " of when each fiber assignment was run and for which tiles, "
                        "along with the options that were used.")
    parser.add_argument("--sky", type=str, required=False,
                        help="Input file with sky or supp_sky targets.  "
                        "These target files are assumed to be constant and not "
                        "tracked by the MTL ledger.")
    parser.add_argument("--mtl", type=str, required=True,
                        help="The MTL ledger.  This is still a work in progress and"
                        " I am not sure what the interface will be, but given the "
                        "fiber assignment dates in the survey log, we should be able"
                        " to get the MTL state at that time.  For now, this option"
                        " is just one or more target files.")
    parser.add_argument("--randoms", type=str, required=True,
                        help="FITS file containing random targets.")
    parser.add_argument("--footprint", type=str, required=False, default=None,
                        help="Optional FITS file defining the footprint.  If"
                        " not specified, the default footprint from desimodel"
                        " is used.")
    parser.add_argument("--tiles", type=str, required=False, default=None,
                        help="Optional text file containing a subset of the"
                        " tile IDs to use in the footprint, one ID per line."
                        " Default uses all tiles in the footprint.")
    parser.add_argument("--format", type=str, required=False, default="fits",
                        help="File format for outputs (either fits or hdf5).")
    parser.add_argument("--outdir", type=str, required=False, default=".",
                        help="Output directory.")
    parser.add_argument("--realizations", type=int, required=False, default=128,
                        help="Number of realizations.")
    args = parser.parse_args()

    # Read tiles we are using
    tileselect = None
    if args.tiles is not None:
        tileselect = list()
        with open(args.tiles, "r") as f:
            for line in f:
                # Try to convert the first column to an integer.
                try:
                    tileselect.append(int(line.split()[0]))
                except ValueError:
                    pass
    tiles = load_tiles(tiles_file=args.footprint,
                       select=tileselect)

    # Load mtl
    mtl = Table.read(args.mtl)
    mtl['SUBPRIORITY'] = np.ones(len(mtl))
    mtl['OBSCONDITIONS'] = np.ones(len(mtl), dtype=int)
    mtl['DESI_TARGET'] = np.ones(len(mtl), dtype=int)

    # Load science targets and get IDs
    tgs = Targets()
    load_target_table(tgs, mtl)
    tg_science = tgs.ids()
    n_target = len(tg_science)

    # Load sky targets
    if args.sky:
        sky = Table.read(args.sky)
        load_target_table(tgs, sky)

    # Divide up realizations among the processes.
    n_realization = args.realizations
    realizations = np.arange(n_realization, dtype=np.int32)
    my_realizations = np.array_split(realizations, mpi_procs)[mpi_rank]

    # Bitarray for all targets and realizations
    bitweights = np.zeros(n_realization * n_target, dtype=bool)

    hw = load_hardware()

    for realization in my_realizations:
        # Set the seed based on the realization, so that the result is reproducible
        # regardless of which process is working on the realization.
        np.random.seed(realization)

        # Randomize science target subpriority for this realization
        mtl['SUBPRIORITY'] = np.random.random_sample(size=n_target)

        # Update targets given new subpriorities
        tgs = Targets()
        load_target_table(tgs, mtl)
        tree = TargetTree(tgs)

        # Compute available targets / locations
        tgsavail = TargetsAvailable(hw, tgs, tiles, tree)
        favail = LocationsAvailable(tgsavail)
        asgn = Assignment(tgs, tgsavail, favail)

        # Run fiber assignment for this realization
        run(asgn)

        # Update bit arrays for assigned science targets
        assigned = []
        for tile_id in tiles.id:
            try:
                adata = asgn.tile_location_target(tile_id)
                for loc, tgid in adata.items():
                    assigned.append(tgid)
            except:
                pass
        idas = np.isin(tg_science, assigned)
        bitweights[realization*n_target:(realization+1)*n_target] = idas

    # Gather weights from all processes
    gather_weights = None
    if mpi_rank == 0:
        gather_weights = np.empty(len(bitweights), dtype=bool)
    min_targ = np.min(my_realizations) * n_target
    max_targ = (np.max(my_realizations) + 1) * n_target
    comm.Gather(bitweights[min_targ:max_targ], gather_weights, root=0)

    if mpi_rank == 0:
        # Pack weights into 64 bit integers
        weights = np.array_split(gather_weights, n_realization)
        bweights = np.array(weights).T
        n_vector = n_realization // 64
        bitvectors = [[] for _ in range(n_vector)]
        for w in bweights:
            for v in range(n_vector):
                bitvectors[v].append(np.packbits(list(w)).view(np.int)[v])

        # Load randoms
        randoms = Table.read(args.randoms)

        # Write output
        try:
            os.makedirs(args.outdir)
        except:
            pass

        if args.format == 'fits':
            # Output targeted file
            tfile = os.path.join(args.outdir, 'targeted.fits')
            toutput = Table()
            toutput['TARGETID'] = mtl['TARGETID'][idas]
            toutput['RA'] = mtl['RA'][idas]
            toutput['DEC'] = mtl['DEC'][idas]
            toutput['Z'] = mtl['Z'][idas]
            for t,vec in enumerate(bitvectors):
                toutput['BITWEIGHT{}'.format(t)] = np.array(vec)[idas]
            toutput.write(tfile)

            # Output parent file
            pfile = os.path.join(args.outdir, 'parent.fits')
            poutput = Table()
            poutput['TARGETID'] = mtl['TARGETID']
            poutput['RA'] = mtl['RA']
            poutput['DEC'] = mtl['DEC']
            poutput['Z'] = mtl['Z']
            for p,vec in enumerate(bitvectors):
                poutput['BITWEIGHT{}'.format(p)] = -np.ones(len(vec), dtype=int)
            poutput.write(pfile)

            # Output randoms file
            rfile = os.path.join(args.outdir, 'randoms.fits')
            routput = Table()
            routput['RA'] = randoms['RA']
            routput['DEC'] = randoms['DEC']
            routput['Z'] = randoms['Z']
            for r in range(len(bitvectors)):
                routput['BITWEIGHT{}'.format(r)] = -np.ones(len(randoms), dtype=int)
            routput.write(rfile)

        elif args.format == 'hdf5':
            # Output targeted file
            targetfile = os.path.join(args.outdir, 'targeted.hdf5')
            tfile = h5py.File(targetfile, 'w')
            tfile.create_dataset('TARGETID', data=mtl['TARGETID'][idas])
            tfile.create_dataset('RA', data=mtl['RA'][idas])
            tfile.create_dataset('DEC', data=mtl['DEC'][idas])
            tfile.create_dataset('Z', data=mtl['Z'][idas])
            for t,vec in enumerate(bitvectors):
                tfile.create_dataset('BITWEIGHT{}'.format(t), data=np.array(vec)[idas])
            tfile.close()

            # Output parent file
            parentfile = os.path.join(args.outdir, 'parent.hdf5')
            pfile = h5py.File(parentfile, 'w')
            pfile.create_dataset('TARGETID', data=mtl['TARGETID'])
            pfile.create_dataset('RA', data=mtl['RA'])
            pfile.create_dataset('DEC', data=mtl['DEC'])
            pfile.create_dataset('Z', data=mtl['Z'])
            for p,vec in enumerate(bitvectors):
                pfile.create_dataset('BITWEIGHT{}'.format(p), data=-np.ones(len(vec), dtype=int))
            pfile.close()

            # Output randoms file
            randomfile = os.path.join(args.outdir, 'randoms.hdf5')
            rfile = h5py.File(randomfile, 'w')
            rfile.create_dataset('RA', data=randoms['RA'])
            rfile.create_dataset('DEC', data=randoms['DEC'])
            rfile.create_dataset('Z', data=randoms['Z'])
            for r in range(len(bitvectors)):
                rfile.create_dataset('BITWEIGHT{}'.format(r), data=-np.ones(len(randoms), dtype=int))
            rfile.close()

        else:
            raise TypeError("Must specify either fits or hdf5 as output file format.")

if __name__ == "__main__":
    main()
