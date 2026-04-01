##!/global/common/software/desi/perlmutter/desiconda/20240425-2.2.0/conda/bin/python -u
from desiutil.iers import freeze_iers
freeze_iers()

#keeping all common imports
import multiprocessing as mp
from multiprocessing import Pool
import atexit

import desitarget
#from desitarget import io

#import desitarget.io as io
import glob
from LSS.SV3.altmtltools import initializeAlternateMTLs
import numpy as np
import os 
from sys import argv
from desiutil.log import get_logger
import cProfile, pstats, io
from pstats import SortKey
import argparse
log = get_logger()

#Keeping all arguments from previous version
parser = argparse.ArgumentParser(
                    prog = 'InitializeAltMTLParallel',
                    description = 'Using a set of initial MTLs as an input, create a set of alternate MTLs. More documentation available on the DESI wiki. ')
parser.add_argument('-o', '--outputMTLDirBase', dest='outputMTLDirBase', required=True, type=str, help = 'the path to the location where MTLs will be written initially. If dontusetmp is set, this is also the final location. If it is not set, then finalDir is the final Location.')
parser.add_argument('-obscon', '--obscon', dest='obscon', default='DARK', help = 'observation conditions, either BRIGHT or DARK.', required = False, type = str)
parser.add_argument('-s', '--survey', dest='survey', default='sv3', help = 'DESI survey to create Alt MTLs for. Either sv3 or main.', required = False, type = str)
parser.add_argument('-fd', '--finalDir', dest='finalDir', default=None, help = 'Location to copy final MTLs if not using tmp dirs. if dontusetmp is set, this will be set to $CSCRATCH/$PSCRATCH', required = False, type = str)
parser.add_argument('-sd', '--startDate', dest='startDate', default=None, help = 'Date at which to start the MTL Loop. Format should be either YYYYMMDD or the rundate format of YYYY-MM-DDTHH:MM:SS+00:00. WARNING: Currently does NOT retain MTL updates from before startDate.', required = False, type = str)
parser.add_argument('-ed', '--endDate', dest='endDate', default=None, help = 'Date at which to end the MTL Loop. Format should be either YYYYMMDD or the rundate format of YYYY-MM-DDTHH:MM:SS+00:00.', required = False, type = str)
parser.add_argument('-seed', '--seed', dest='seed', default=31415, help = 'Random seed to ensure reproducability.', required = False, type = int)
parser.add_argument('-n', '--ndir', dest='ndir', default=128, help = 'Random seed to ensure reproducability.', required = False, type = int)
parser.add_argument('-v', '--verbose', dest = 'verbose', default=False, action='store_true', help = 'set flag to enter verbose mode')
parser.add_argument('-prof', '--profile', dest = 'profile', default=False, action='store_true', help = 'set flag to profile code time usage.')
parser.add_argument('-d', '--debug', dest = 'debug', default=False, action='store_true', help = 'set flag to enter debug mode.')
parser.add_argument('-nt', '--dontUseTemp', dest = 'usetmp', default=True, action='store_false', help = 'pass flag to NOT use local tmp directories on compute nodes for initial writing of alt MTLs.')
parser.add_argument('-ow', '--overwrite', dest = 'overwrite', default=False, action='store_true', help = 'pass this flag to regenerate already existing alt MTLs.')
parser.add_argument('-sb', '--saveBackup', dest = 'saveBackup', default=False, action='store_true', help = 'Save the initial MTLs post shuffling in a backup directory underneath the main alt MTL directories.')
parser.add_argument('-hpfn', '--HPListFile', dest='HPListFile', default=None, help = 'Name of a text file consisting only of one line of comma separated healpixel numbers for which the code will generate alt MTLs. If not specified, it will be automatically determined from the survey name.', required = False, type = str)

parser.add_argument('--version', action='version', version='There are no version numbers.', help = 'There are no version numbers.')
parser.add_argument('-sbp', '--shuffleBrightPriorities', action='store_true', dest='shuffleBrightPriorities', default=False, help = 'Pass this flag to shuffle top level PRIORITIES in BRIGHT survey Alt MTLs. If this argument is true, should also be setting PromoteFracBGSFaint.', required = False)
parser.add_argument('-sep', '--shuffleELGPriorities', action='store_true', dest='shuffleELGPriorities', default=False, help = 'Pass this flag to shuffle top level ELG PRIORITIES in DARK survey Alt MTLs. If this argument is true, should also be setting PromoteFracELG.', required = False)
parser.add_argument('-dssp', '--dontShuffleSubpriorities', action='store_false', dest='shuffleSubpriorities', default=True, help = 'WARNING: THIS FLAG SHOULD ONLY BE USED FOR DEBUGGING. Pass this flag to NOT shuffle subpriorities. This option must be used in conjunction with --reproducing.', required = False)
parser.add_argument('-rep', '--reproducing', action='store_true', dest='reproducing', default=False, help = 'WARNING: THIS FLAG SHOULD ONLY BE USED FOR DEBUGGING. Pass this flag to confirm to the alt mtl code that you are trying to reproduce real MTLs. This option should (must?) be used in conjunction with --shuffleSubpriorities.', required = False)
parser.add_argument('-pfbf', '--PromoteFracBGSFaint', dest='PromoteFracBGSFaint', default=0.2, help = 'What (decimal) percentage of BGS_FAINT targets to promote to BGS_FAINT_HIP. This argument is only used if shuffleBrightPriorities is passed and if obscon = BRIGHT.', required = False, type = float)
parser.add_argument('-pfe', '--PromoteFracELG', dest='PromoteFracELG', default=0.1, help = 'What (decimal) percentage of ELG_LOP targets to promote to ELG_HIP. This argument is only used if shuffleELGPriorities is passed and if obscon = DARK.', required = False, type = float)
parser.add_argument('-elb', '--exampleLedgerBase', dest='exampleLedgerBase', default='/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/', help = 'Location of the real (or mock) MTLs that serve as the basis for the alternate MTLs. Defaults to location of data MTLs. Do NOT include survey or obscon information here. ', required = False, type = str)
parser.add_argument('-p2L', '--path2LSS', dest='path2LSS', default='~/.local/desicode/LSS/', help = 'location where LSS repository is cloned.', required = False, type = str)

parser.add_argument('-ppn', '--ProcPerNode', dest='ProcPerNode', default=None, help = 'Number of processes to spawn per requested node. If not specified, determined automatically from NERSC_HOST.', required = False, type = int)
args = parser.parse_args()

#keeping overwriting check
if ('trunk' in args.outputMTLDirBase.lower()) or  ('ops' in args.outputMTLDirBase.lower()):
    raise ValueError("In order to prevent accidental overwriting of the real MTLs, please remove \'ops\' and \'trunk\' from your MTL output directory")

#simplifying to just perlmutter
if args.ProcPerNode is None:
    if 'perlmutter' in os.getenv['NERSC_HOST'].lower():
        args.ProcPerNode = 128
    else:
        raise ValueError('Code is only supported on NERSC Cori and NERSC perlmutter.')

# If folder doesn't exist, then create it.
if not os.path.isdir(args.outputMTLDirBase):
    os.makedirs(args.outputMTLDirBase)

if args.usetmp:
    testdir = args.finalDir
else:
    testdir = args.outputMTLDirBase

if not os.path.isdir(testdir):
    os.makedirs(testdir)
if os.path.exists(testdir + 'SeedFile'):
    temp = open(testdir + 'SeedFile', 'r')
    tempseed = int(temp.readlines()[0].split()[0])
    log.info('Seed file already exists with seed {0:d}'.format(tempseed))

    if int(tempseed) == int(args.seed):
        log.info('random seeds are saved, continuing')
    else:
        raise RuntimeError('different random seed {0:d} provided than for initial generation {1:d}'.format(int(seed), int(tempseed)))
else:
    with open(args.outputMTLDirBase + 'SeedFile', 'w') as f:
        f.write(str(args.seed))

HPList = np.array(open(args.HPListFile,'r').readlines()[0].split(',')).astype(int)
if args.verbose or args.debug:
    log.debug('HPList')
    log.debug(HPList)
    log.debug('First healpixel: {0:d}'.format(HPList[0]))
    log.debug('Last healpixel: {0:d}'.format(HPList[-1]))
    log.debug('Number of healpixels: {0:d}'.format(int(len(HPList))))

#Dropping these
#NodeID = int(os.getenv('SLURM_NODEID'))
#SlurmNProcs = int(os.getenv('SLURM_NPROCS'))

NNodes = 1
NProc = int(NNodes*args.ProcPerNode)

outputMTLDir = args.outputMTLDirBase + "Univ{0:03d}/"

def procFunc(nproc):
    if 'sv' in args.survey.lower():
        log.info('sv survey')
        mtlprestr = args.survey.lower()
    else:
        log.info('non sv survey')
        mtlprestr = ''

    if args.usetmp:
        finalDir = args.finalDir.format(nproc)
    else:
        #TEST THIS
        finalDir = outputMTLDir.format(nproc)

    if os.path.exists(finalDir + '/{0}/{2}/{3}mtl-{2}-hp-{1}.ecsv'.format(args.survey.lower(),HPList[-1], args.obscon.lower(), mtlprestr)):
        log.info('Alt MTL for last HP in list exists. Exiting script')
        log.info(outputMTLDir + '/{0}/{2}/{3}mtl-{2}-hp-{1}.ecsv'.format(args.survey.lower(),HPList[-1], args.obscon.lower(), mtlprestr))
        return 42
    for hpnum in HPList:
        if os.path.exists(finalDir + '/{0}/{2}/{3}mtl-{2}-hp-{1}.ecsv'.format(args.survey.lower(),hpnum, args.obscon.lower(), mtlprestr)) and (not args.overwrite):
            log.info('Alt MTL for HP {0:d} already exists. Set -ow or --overwrite to force regeneration. '.format(hpnum))
            continue
        log.info('hpnum = {0}'.format(hpnum))
        exampleLedger = args.exampleLedgerBase + '/{0}/{2}/{3}mtl-{2}-hp-{1}.ecsv'.format(args.survey.lower(),hpnum, args.obscon.lower(), mtlprestr)
        log.info('exampleLedger = {0}'.format(exampleLedger))
        log.info('exampleLedgerBase = {0}'.format(args.exampleLedgerBase))
        if args.usetmp and (args.debug or args.verbose or args.profile):
            log.info('outputMTLDir, nproc {0}'.format(nproc))
            log.info(outputMTLDir)
            log.info('finalDir, nproc{0}'.format(nproc))
            log.info(args.finalDir)
        log.info('shuffleSubpriorities: {0}'.format(args.shuffleSubpriorities))
        log.info('reproducing: {0}'.format(args.reproducing))
        initializeAlternateMTLs(exampleLedger, outputMTLDir, genSubset = nproc, seed = args.seed, obscon = args.obscon, survey = args.survey, saveBackup = args.saveBackup, hpnum = hpnum, overwrite = args.overwrite, reproducing = args.reproducing, shuffleSubpriorities = args.shuffleSubpriorities, startDate=args.startDate, endDate=args.endDate, profile = args.profile, usetmp=args.usetmp, finalDir=args.finalDir, debug = args.debug, verbose = args.verbose,
            shuffleBrightPriorities = args.shuffleBrightPriorities, shuffleELGPriorities = args.shuffleELGPriorities, PromoteFracBGSFaint = args.PromoteFracBGSFaint
            , PromoteFracELG = args.PromoteFracELG)
    return 0

inds = []
start = 0
end = args.ndir
if args.verbose or args.debug:
    log.debug('start')
    log.debug(start)
    log.debug('end')
    log.debug(end)
if args.ndir < start:
    raise ValueError('ndir is too low for the number of nodes requested. Either request more realizations (ndir) or fewer nodes')
for i in range(start, end):
    if i >= args.ndir: 
        break
    if args.debug or args.verbose or args.profile:
        log.info('i')
        log.info(i)
    inds.append(i)
    
NProc = len(inds)
assert(len(inds))

log.info('running on NProc = {0} processes'.format(NProc))
p = Pool(NProc)
atexit.register(p.close)
log.info('running procFunc now on inds:')
log.info(inds)

result = p.map(procFunc,inds)













































