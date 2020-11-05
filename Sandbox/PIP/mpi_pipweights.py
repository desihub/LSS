#!/usr/bin/env python

"""
Example of how one could compute PIP weights in parallel.
"""

from mpi4py import MPI

import os
import sys

import argparse

from bitarray import bitarray

import fitsio

import fiberassign

from fiberassign.utils import Logger, distribute_discrete

from fiberassign.hardware import load_hardware

from fiberassign.tiles import load_tiles

from fiberassign.targets import (
    Targets,
    TargetsAvailable,
    TargetTree,
    LocationsAvailable,
    load_target_file
)

from fiberassign.assign import Assignment, run



def main():
    log = Logger.get()

    mpi_procs = MPI.COMM_WORLD.size
    mpi_rank = MPI.COMM_WORLD.rank

    parser = argparse.ArgumentParser()

    parser.add_argument("--survey_log", type=str, required=False,
                        help="Eventually we would pass in a file containing the log"
                        " of when each fiber assignment was run and for which tiles, "
                        "along with the options that were used.")

    parser.add_argument("--sky", type=str, required=False, nargs="+",
                        help="Input file with sky or supp_sky targets.  "
                        "These target files are assumed to be constant and not "
                        "tracked by the MTL ledger.")

    parser.add_argument("--mtl", type=str, required=True, nargs="+",
                        help="The MTL ledger.  This is still a work in progress and"
                        " I am not sure what the interface will be, but given the "
                        "fiber assignment dates in the survey log, we should be able"
                        " to get the MTL state at that time.  For now, this option"
                        " is just one or more target files.")

    parser.add_argument("--footprint", type=str, required=False, default=None,
                        help="Optional FITS file defining the footprint.  If"
                        " not specified, the default footprint from desimodel"
                        " is used.")

    parser.add_argument("--tiles", type=str, required=False, default=None,
                        help="Optional text file containing a subset of the"
                        " tile IDs to use in the footprint, one ID per line."
                        " Default uses all tiles in the footprint.")

    parser.add_argument("--out", type=str, required=False, default=None,
                        help="Output directory.")

    parser.add_argument("--realizations", type=int, required=False, default=10,
                        help="Number of realizations.")

    args = parser.parse_args()

    if args.sky is None:
        args.sky = list()

    # Set output directory
    if args.out is None:
        args.out = "."

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
    tiles = load_tiles(
        tiles_file=args.footprint,
        select=tileselect,
    )

    # Create empty target list
    tgs = Targets()

    # Append each input target file.  These target files must all be of the
    # same survey type, and will set the Targets object to be of that survey.

    for tgfile in args.targets:
        load_target_file(tgs, tgfile)

    # Just the science target IDs
    tg_science = tgs.ids()
    tg_science2indx = {y: x for x, y in enumerate(tg_science)}

    # Now load the sky target files.
    survey = tgs.survey()
    for tgfile in args.sky:
        load_target_file(tgs, tgfile)

    # Divide up realizations among the processes.

    n_realization = args.realizations
    realizations = np.arange(n_realization, dtype=np.int32)
    my_realizations = np.array_split(realizations, mpi_procs)[mpi_rank]

    # Bitarray for all targets and realizations
    tgarray = bitarray(len(tg_science) * n_realization)
    tgarray.setall(False)

    # Target tree
    tree = TargetTree(tgs)

    for realization in my_realizations:
        # Set the seed based on the realization, so that the result is reproducible
        # regardless of which process is working on the realization.
        np.random.seed(realization)

        # Comment out the next block to avoid randomizing subpriority
        # ----
        # Randomize science target subpriority for this realization
        new_subpriority = np.random.random_sample(size=len(tg_science))
        for indx, tgid in enumerate(tg_science):
            tg = tgs.get(tgid)
            tg.subpriority = new_subpriority[indx]

        # Comment out the next block to avoid dithering tiles
        # ----
        # Dither tiles centers by the same



        # Compute available targets / locations

        tgsavail = TargetsAvailable(hw, tgs, tiles, tree)

        favail = LocationsAvailable(tgsavail)

        asgn = Assignment(tgs, tgsavail, favail)

        # Replay the survey log for each time fiber assignment was run.  For now, this
        # is just doing the tiles all at once.
        for assign_event in range(1):
            # In the future, load MTL updates to the obs remaining for each target here

            # Read hardware properties- in the future, pass in the assignment run date
            # to this function.
            hw = load_hardware()

            # Run assignment for this event.
            run(asgn)

            # Update bit arrays for assigned science targets
            for tile_id in tiles.ids():
                adata = asgn.tile_location_target(tile_id)
                for loc, tgid in adata.items():
                    try:
                        idx = tg_science2indx[tgid]
                        tgarray[idx * n_realization + realization] = True
                    except KeyError:
                        # Not a science target
                        pass

    # Reduce bitarrays to root process.  The bitarray type conforms to the
    # buffer protocol.

    tgall = None
    if mpi_rank == 0:
        tgall = bitarray(tgarray)
        tgall.setall(False)

    MPI.COMM_WORLD.Reduce(tgarray, tgall, op=MPI.BOR, root=0)

    # Write it out

    if mpi_rank == 0:
        pass


if __name__ == "__main__":
    main()
