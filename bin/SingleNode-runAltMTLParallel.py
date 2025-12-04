from desiutil.iers import freeze_iers
freeze_iers()

from multiprocessing import Pool
from LSS.SV3 import altmtltools as amt
#from LSS.main import mockaltmtltools as mockamt
from astropy.table import Table, vstack, join
#import altmtltools as amt
from desiutil.log import get_logger
#import dill
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
parser.add_argument('-rep', '--reproducing', action='store_true', dest='reproducing', default=False, help = 'WARNING: THIS FLAG SHOULD ONLY BE USED FOR DEBUGGING. Pass this flag to confirm to the alt mtl code that you are trying to reproduce real MTLs. This option should (must?) be used in conjunction with --shuffleSubpriorities.', required = False)
parser.add_argument('-prof', '--profile', dest = 'profile', default=False, action='store_true', help = 'set flag to profile code time usage. This flag may not profile all components of any particular stage of the AMTL pipeline. ')
parser.add_argument('-d', '--debug', dest = 'debug', default=False, action='store_true', help = 'set flag to enter debug mode.')
parser.add_argument('-nfl', '--NumObsNotFromLedger', dest = 'numobs_from_ledger', default=True, action='store_false', help = 'If True (flag is NOT set) then inherit the number of observations so far from the ledger rather than expecting it to have a reasonable value in the zcat.')

parser.add_argument('-redoFA', '--redoFA', dest = 'redoFA', default=False, action='store_true', help = 'pass this flag to regenerate already existing fiber assignment files.')

parser.add_argument('-getosubp', '--getosubp', action='store_true', dest='getosubp', default=False, help = 'WARNING: THIS FLAG SHOULD ONLY BE USED FOR DEBUGGING AND NEVER FOR MOCKS. Pass this flag to grab subpriorities directly from the real survey MTLs for fiberassignment.', required = False)
parser.add_argument('-md', '--multiDate', action='store_true', dest='multiDate', default=False, help = 'Currently this flag is being debugged. In the future, it will switch between interactive submission of each date as a separate job (True) and of all nights to be looped through until a single job`s time runs out. ', required = False)
parser.add_argument('-ppn', '--ProcPerNode', dest='ProcPerNode', default=None, help = 'Number of processes to spawn per requested node. If not specified, determined automatically from NERSC_HOST.', required = False, type = int)
parser.add_argument('-rmbd', '--realMTLBaseDir', dest='mtldir', default='/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/', help = 'Location of the real (or mock) MTLs that serve as the basis for the alternate MTLs. Defaults to location of data MTLs. Do NOT include survey or obscon information here. ', required = False, type = str)
parser.add_argument('-zcd', '--zCatDir', dest='zcatdir', default='/global/cfs/cdirs/desi/spectro/redux/daily/', help = 'Location of the real redshift catalogs for use in alt MTL loop.  Defaults to location of survey zcatalogs.', required = False, type = str)
parser.add_argument('-npc', '--nproc',  dest = 'nproc',required = False, default = 128, help = 'Number of realizations in AltMTL Directory.', type = int)

print(argv)

args = parser.parse_args()

args.ProcPerNode = 128

print(args)



log = get_logger()
'''
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
''';



# These should be constant/default. 
# If this changes, add these to argparse

ndirs = None
multiproc = True
singleDate = not(args.multiDate)


def procFunc(nproc):
    return
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
#    if args.mock:
    print(args.obscon, args.survey, args.mtldir, args.zcatdir, args.altMTLBaseDir, ndirs, args.numobs_from_ledger,args.secondary,args.getosubp, args.quickRestart, multiproc, nproc, singleDate, args.redoFA, args.mock, targets, args.debug, args.verbose, args.reproducing)
    retval = amt.loop_alt_ledger(args.obscon, survey = args.survey, mtldir = args.mtldir, zcatdir = args.zcatdir, altmtlbasedir = args.altMTLBaseDir, ndirs = ndirs, numobs_from_ledger = args.numobs_from_ledger,secondary = args.secondary, getosubp = args.getosubp, quickRestart = args.quickRestart, multiproc = multiproc, nproc = nproc, singleDate = singleDate, redoFA = args.redoFA, mock = args.mock, targets = targets, debug = args.debug, verbose = args.verbose, reproducing = args.reproducing)
        #retval = mockamt.loop_alt_ledger(args.obscon, survey = args.survey, mtldir = args.mtldir, zcatdir = args.zcatdir, altmtlbasedir = args.altMTLBaseDir, ndirs = ndirs, numobs_from_ledger = args.numobs_from_ledger,secondary = args.secondary, getosubp = args.getosubp, quickRestart = args.quickRestart, multiproc = multiproc, nproc = nproc, singleDate = singleDate, redoFA = args.redoFA, mock = args.mock, targets = targets, debug = args.debug, verbose = args.verbose, reproducing = args.reproducing)
    #else:
    #    retval = amt.loop_alt_ledger(args.obscon, survey = args.survey, mtldir = args.mtldir, zcatdir = args.zcatdir, altmtlbasedir = args.altMTLBaseDir, ndirs = ndirs, numobs_from_ledger = args.numobs_from_ledger,secondary = args.secondary, getosubp = args.getosubp, quickRestart = args.quickRestart, multiproc = multiproc, nproc = nproc, singleDate = singleDate, redoFA = args.redoFA, mock = args.mock, targets = targets, debug = args.debug, verbose = args.verbose, reproducing = args.reproducing)
    if args.verbose:
        log.debug('finished with one iteration of procFunc')
    if type(retval) == int:
        if args.verbose:
            log.debug('retval')
            log.debug(retval)
        if retval == 151:
            log.info('No more data. Ending script.')
            return 151
        return retval
    elif args.verbose:
        print('retval')
        print(retval)
    
    return 42

#indices of univ to update
if args.nproc == 1:
    inds = [0]
else:
    inds = np.arange(args.nproc)

print(inds)
    
#assert(len(inds))
#p = Pool()
#atexit.register(p.close)
#result = p.map(procFunc,inds)


#run as 
#python3 python3 /global/homes/d/desica/LSScode/LSS/bin/SingleNode-runAltMTLParallel.py --altMTLBaseDir=/pscratch/sd/d/desica/AltMTL -obscon=DARK --survey=main --realMTLBaseDir='/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/' --multiDate --nproc 128