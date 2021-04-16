#!/usr/bin/env python

"""
Example of how one could compute PIP weights in parallel.
"""

from mpi4py import MPI
import os, sys
import argparse
import numpy as np
from astropy.table import Table
from fiberassign.utils import Logger
from fiberassign.hardware import load_hardware
from fiberassign.tiles import load_tiles
from fiberassign.targets import (Targets,
                                 TargetsAvailable,
                                 TargetTree,
                                 LocationsAvailable,
                                 load_target_table)
from fiberassign.assign import Assignment, run


def main():
    log = Logger.get()

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
    parser.add_argument("--truth", type=str, required=False,
                        help="Truth information used to access and output redshift.")
    parser.add_argument("--footprint", type=str, required=False, default=None,
                        help="Optional FITS file defining the footprint.  If"
                        " not specified, the default footprint from desimodel"
                        " is used.")
    parser.add_argument("--tiles", type=str, required=False, default=None,
                        help="Optional text file containing a subset of the"
                        " tile IDs to use in the footprint, one ID per line."
                        " Default uses all tiles in the footprint.")
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

    # Load mtl and truth, make sure targets are the same
    mtl = Table.read(args.mtl)
    if not 'SUBPRIORITY' in mtl.keys():
        mtl['SUBPRIORITY'] = np.ones(len(mtl))
        mtl['OBSCONDITIONS'] = np.ones(len(mtl), dtype=int)
        mtl['DESI_TARGET'] = np.ones(len(mtl), dtype=int)
    if args.truth:
        truth = Table.read(args.truth)
        assert mtl['TARGETID'].all() == truth['TARGETID'].all(), 'MTL and truth targets are different'

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
    bitweights = np.zeros((n_realization, n_target),dtype=bool)

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
        bitweights[realization][idas] = True

    # Gather weights from all processes
    # Number of mpi processes must equal number of realizations (TODO: fix this)
    gather_weights = None
    if mpi_rank == 0:
        gather_weights = np.empty((n_realization, n_target), dtype=bool)
    comm.Gather(bitweights[mpi_rank], gather_weights, root=0)

    if mpi_rank == 0:
        # Pack bits into 64 bit integers
        bweights = np.array(gather_weights).T
        n_vector = n_realization // 64
        bitvectors = [[] for _ in range(n_vector)]
        for w in bweights:
            for v in range(n_vector):
                bitvectors[v].append(np.packbits(list(w)).view(np.int)[v])

        # Get redshifts and spectral type
        if args.truth:
            z = truth['TRUEZ']
            templatetype = truth['TEMPLATETYPE']
            templatetype = np.array([t.strip() for t in templatetype], dtype=str)
        else:
            z = mtl['Z']

        # Write output
        try:
            os.makedirs(args.outdir)
        except:
            pass

        outfile = os.path.join(args.outdir,'bitweight_vectors.fits')
        output = Table()
        output['TARGETID'] = mtl['TARGETID']
        output['RA'] = mtl['RA']
        output['DEC'] = mtl['DEC']
        output['Z'] = z
        output['ASSIGNEDID'] = idas
        if args.truth:
            output['TEMPLATETYPE'] = templatetype
        for i,vec in enumerate(bitvectors):
            output['BITWEIGHT{}'.format(i)] = vec
        output.write(outfile)

if __name__ == "__main__":
    main()
