#!/global/common/software/desi/perlmutter/desiconda/20230111-2.1.0/conda/bin/python -u
from desiutil.iers import freeze_iers
freeze_iers()

from multiprocessing import Pool
from LSS.SV3 import altmtltools as amt
from astropy.table import Table, vstack, join
#import altmtltools as amt
from desiutil.log import get_logger
import dill
from sys import argv
import os
import numpy as np
import multiprocessing as mp
import logging
import atexit
import glob
import cProfile, pstats, io
from pstats import SortKey
import argparse

#Base directory for the alternate MTLs created in the InitializeAltMTLs script


parser = argparse.ArgumentParser(
                    prog = 'RunAltMTLParallel',
                    description = 'Progresses alternate MTLs through the MTL update loop in parallel. More documentation available on the DESI wiki. ')
parser.add_argument('-a', '--altMTLBaseDir', dest='altMTLBaseDir', required=True, type=str, help = 'the path to the location where alt MTLs are stored, up to, but not including survey and obscon information.')

parser.add_argument('-obscon', '--obscon', dest='obscon', default='DARK', help = 'observation conditions, either BRIGHT or DARK.', required = False, type = str)
parser.add_argument('-s', '--survey', dest='survey', default='sv3', help = 'DESI survey to create Alt MTLs for. Either sv3 or main.', required = False, type = str)
parser.add_argument('-sec', '--secondary', dest = 'secondary', default=False, action='store_true', help = 'set flag to incorporate secondary targets.')
parser.add_argument('-mock', '--mock', dest = 'mock', default=False, action='store_true', help = 'set flag if running pipeline on mocks.')
parser.add_argument('-tf', '--targfile', dest='targfile', required=False, default = None, type=str, help = 'Location for target file for mocks or data. Only required if mocks are being processed.')
parser.add_argument('-v', '--verbose', dest = 'verbose', default=False, action='store_true', help = 'set flag to enter verbose mode')
parser.add_argument('-qr', '--quickRestart', dest = 'quickRestart', default=False, action='store_true', help = 'set flag to remove any AMTL updates that have already been performed. Useful for rapidfire debugging of steps in this part of the pipeline.')
parser.add_argument('-prof', '--profile', dest = 'profile', default=False, action='store_true', help = 'set flag to profile code time usage. This flag may not profile all components of any particular stage of the AMTL pipeline. ')
parser.add_argument('-d', '--debug', dest = 'debug', default=False, action='store_true', help = 'set flag to enter debug mode.')
parser.add_argument('-nfl', '--NumObsNotFromLedger', dest = 'numobs_from_ledger', default=True, action='store_false', help = 'If True (flag is NOT set) then inherit the number of observations so far from the ledger rather than expecting it to have a reasonable value in the zcat.')

parser.add_argument('-redoFA', '--redoFA', dest = 'redoFA', default=False, action='store_true', help = 'pass this flag to regenerate already existing fiber assignment files.')

parser.add_argument('-getosubp', '--getosubp', action='store_true', dest='getosubp', default=False, help = 'WARNING: THIS FLAG SHOULD ONLY BE USED FOR DEBUGGING AND NEVER FOR MOCKS. Pass this flag to grab subpriorities directly from the real survey MTLs for fiberassignment.', required = False)
parser.add_argument('-ppn', '--ProcPerNode', dest='ProcPerNode', default=None, help = 'Number of processes to spawn per requested node. If not specified, determined automatically from NERSC_HOST.', required = False, type = int)
parser.add_argument('-rmbd', '--realMTLBaseDir', dest='mtldir', default='/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/', help = 'Location of the real (or mock) MTLs that serve as the basis for the alternate MTLs. Defaults to location of data MTLs. Do NOT include survey or obscon information here. ', required = False, type = str)
parser.add_argument('-zcd', '--zCatDir', dest='zcatdir', default='/global/cfs/cdirs/desi/spectro/redux/daily/', help = 'Location of the real redshift catalogs for use in alt MTL loop.  Defaults to location of survey zcatalogs.', required = False, type = str)

print(argv)

args = parser.parse_args()
log = get_logger()

if args.mock:
    assert(not (args.targfile is None))
    print('args.getosubp')
    print(args.getosubp)
    assert(not (args.getosubp))

# Leave confirmation file in output directory if using original subpriorities
if args.getosubp:
    from pathlib import Path
    Path(args.altMTLBaseDir + '/GETOSUBPTRUE').touch()

#Get information about environment for multiprocessing
NodeID = int(os.getenv('SLURM_NODEID'))
SlurmNProcs = int(os.getenv('SLURM_NPROCS'))
try:
    NNodes = int(os.getenv('SLURM_JOB_NUM_NODES'))
except:
    log.warning('no SLURM_JOB_NUM_NODES env set. You may not be on a compute node.')
    NNodes = 1

if args.ProcPerNode is None:
    if 'cori' in os.getenv['NERSC_HOST'].lower():
        args.ProcPerNode = 32
    elif 'perlmutter' in os.getenv['NERSC_HOST'].lower():
        args.ProcPerNode = 128
    else:
        raise ValueError('Code is only supported on NERSC Cori and NERSC perlmutter.')

NProc = int(NNodes*args.ProcPerNode)
log.info('NProc = {0:d}'.format(NProc))
log.info('NNodes = {0:d}'.format(NNodes))




# These should be constant/default. 
# If this changes, add these to argparse

ndirs = None
multiproc = True
singleDate = True


def procFunc(nproc):
    if args.verbose:
        log.debug('calling procFunc')
    if not(args.targfile is None):
        targets = Table.read(args.targfile)
        print('targets.dtype')
        print(targets.dtype)
        print('targets[0:5]')
        print(targets[0:5])
        print('targets TARGETID,RA,DEC')
        print(targets['TARGETID'][0:5])
        print(targets['RA'][0:5])
        print(targets['DEC'][0:5])
    else:
        targets = None
    retval = amt.loop_alt_ledger(args.obscon, survey = args.survey, mtldir = args.mtldir, zcatdir = args.zcatdir, altmtlbasedir = args.altMTLBaseDir, ndirs = ndirs, numobs_from_ledger = args.numobs_from_ledger,secondary = args.secondary, getosubp = args.getosubp, quickRestart = args.quickRestart, multiproc = multiproc, nproc = nproc, singleDate = singleDate, redoFA = args.redoFA, mock = args.mock, targets = targets, debug = args.debug, verbose = args.verbose)
    if args.verbose:
        log.debug('finished with one iteration of procFunc')
    if type(retval) == int:
        if args.verbose:
            log.debug('retval')
            log.debug(retval)
        if retval == 151:
            raise ValueError('No more data. Ending script.')
        return retval
    elif args.verbose:
        print('retval')
        print(retval)
    
    return 42

inds = []
start = int(NodeID*NProc/SlurmNProcs)
end = int((NodeID + 1)*NProc/SlurmNProcs)
log.info('NodeID = {0:d}'.format(NodeID))
log.info('StartProc = {0:d}'.format(start))
log.info('EndProc = {0:d}'.format(end))


for i in range(start, end):
    log.info('Process i = {0}'.format(i))
    files = glob.glob(args.altMTLBaseDir + "Univ{0:03d}/*".format(i))
    if len(files):
        pass
    else:
        log.info('no files in dir number {0}, not processing that directory.'.format(i))
        continue
    inds.append(i)
    
assert(len(inds))
p = Pool(NProc)
atexit.register(p.close)
result = p.map(procFunc,inds)
