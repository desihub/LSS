#!/global/common/software/desi/cori/desiconda/20211217-2.0.0/conda/bin/python -u
import multiprocessing as mp
from multiprocessing import Pool
import logging
import atexit
import desitarget.io as io
import glob
from LSS.SV3.altmtltools import initializeAlternateMTLs
import numpy as np
import os 
from sys import argv
from desiutil.log import get_logger
import cProfile, pstats, io
from pstats import SortKey
profile = False
if profile:
    pr = cProfile.Profile()
    pr.enable()

log = get_logger()
debug = False #These will eventually be CLAs
verbose = False  #These will eventually be CLAs
ProcPerNode = 32 #These will eventually be CLAs
#startDate = '2021-11-19T20:43:24+00:00'
startDate = None #These will eventually be CLAs
if debug or verbose or profile:
    log.info('CLAs for script')
    log.info(argv)
seed = int(argv[1])
ndir = int(argv[2])
try:
    overwrite = bool(int(argv[3]))
except:
    raise ValueError('Invalid non-integer value of overwrite: {0}'.format(argv[3]))
obscon = argv[4]
survey = argv[5]
outputMTLDirBase = argv[6]
if ('trunk' in outputMTLDirBase.lower()) or  ('ops' in outputMTLDirBase.lower()):
    raise ValueError("In order to prevent accidental overwriting of the real MTLs, please remove \'ops\' and \'trunk\' from your MTL output directory")
log.info('output directory for alternate MTLs: {0}'.format(outputMTLDirBase))
HPListFile = argv[7]
try:
    shuffleBrightPriorities = bool(int(argv[8]))
except:
    if obscon.lower() == 'dark':
        log.info('Ignoring invalid noninteger value of shuffleBrightPriorities: {0} because dark time MTLs are being initialized'.format(argv[8]))
        shuffleBrightPriorities = False
    else:
        raise ValueError('Invalid non-integer value of shuffleBrightPriorities: {0}'.format(argv[8]))

try:
    PromoteFracBGSFaint = float(argv[9])
except:
    if obscon.lower() == 'dark':
        log.info('Ignoring invalid nonfloat value of PromoteFracBGSFaint: {0} because dark time MTLs are being initialized'.format(argv[9]))
        PromoteFracBGSFaint = 0.2
    else:
        raise ValueError('Invalid non-float value of PromoteFracBGSFaint: {0}'.format(argv[9]))
        PromoteFracBGSFaint = 0.2

exampleledgerbase = argv[10] #"/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/{0}/{2}/{3}mtl-{2}-hp-{1}.ecsv"
NNodes = int(argv[11])
# If folder doesn't exist, then create it.
if not os.path.isdir(outputMTLDirBase):
    os.makedirs(outputMTLDirBase)
if os.path.exists(outputMTLDirBase + 'SeedFile'):
    temp = open(outputMTLDirBase + 'SeedFile', 'r')
    tempseed = int(temp.readlines()[0].split()[0])
    log.info('Seed file already exists with seed {0:d}'.format(tempseed))

    if int(tempseed) == int(seed):
        log.info('random seeds are saved, continuing')
    else:
        raise RuntimeError('different random seed {0:d} provided than for initial generation {1:d}'.format(int(seed), int(tempseed)))
else:
    with open(outputMTLDirBase + 'SeedFile', 'w') as f:
        f.write(str(seed))

HPList = np.array(open(HPListFile,'r').readlines()[0].split(',')).astype(int)
log.info('First healpixel: {0:d}'.format(HPList[0]))
log.info('Last healpixel: {0:d}'.format(HPList[-1]))
log.info('Number of healpixels: {0:d}'.format(int(len(HPList))))

NodeID = int(os.getenv('SLURM_NODEID'))
SlurmNProcs = int(os.getenv('SLURM_NPROCS'))

NProc = int(NNodes*ProcPerNode)

log.info('requested number of nodes: {0:d}'.format(NNodes))
log.info('requested number of directories/realizations: {0:d}'.format(ndir))
log.info('requested number of processes: {0:d}'.format(NProc))

outputMTLDir = outputMTLDirBase + "Univ{0:03d}/"

HPList = np.array(open(HPListFile,'r').readlines()[0].split(',')).astype(int)

usetmp = argv[12] #outputMTLDir.startswith('/dev/shm/') | outputMTLDir.startswith('/tmp/')
if usetmp:
    finalDir = argv[13]
else:
    finalDir = None

try:
    shuffleSubpriorities = bool(int(argv[14]))
except:
    raise ValueError('Invalid non-integer value of shuffleSubpriorities: {0}'.format(argv[14]))


try:
    reproducing = bool(int(argv[15]))
except:
    raise ValueError('Invalid non-integer value of reproducing: {0}'.format(argv[15]))

def procFunc(nproc):
    if 'sv' in survey.lower():
        log.info('sv survey')
        mtlprestr = survey.lower()
    else:
        log.info('non sv survey')
        mtlprestr = ''

    if os.path.exists(outputMTLDir + '/{0}/{2}/{3}mtl-{2}-hp-{1}.ecsv'.format(survey.lower(),HPList[-1], obscon.lower(), mtlprestr)):
        log.info('Alt MTL for last HP in list exists. Exiting script')
        log.info(outputMTLDir + '/{0}/{2}/{3}mtl-{2}-hp-{1}.ecsv'.format(survey.lower(),HPList[-1], obscon.lower(), mtlprestr))
        return 42
    for hpnum in HPList:
        exampleledger = exampleledgerbase + '/{0}/{2}/{3}mtl-{2}-hp-{1}.ecsv'.format(survey.lower(),hpnum, obscon.lower(), mtlprestr)
        if usetmp and (debug or verbose or profile):
            log.info('outputMTLDir, nproc {0}'.format(nproc))
            log.info(outputMTLDir)
            log.info('finalDir, nproc{0}'.format(nproc))
            log.info(finalDir)
        initializeAlternateMTLs(exampleledger, outputMTLDir, genSubset = nproc, seed = seed, obscon = obscon, survey = survey, saveBackup = True, hpnum = hpnum, overwrite = overwrite, reproducing = reproducing, shuffleSubpriorities = shuffleSubpriorities, startDate = startDate, profile = profile, usetmp=usetmp, finalDir=finalDir, debug = debug, verbose = verbose)
    return 0
inds = []
start = int(NodeID*NProc/SlurmNProcs)
end = int((NodeID + 1)*NProc/SlurmNProcs)
print('start')
print(start)
print('end')
print(end)
if ndir < start:
    raise ValueError('ndir is too low for the number of nodes requested. Either request more realizations (ndir) or fewer nodes')
for i in range(start, end):
    if i >= ndir: 
        break
    if debug or verbose or profile:
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
if profile:
    pr.disable()
    s = io.StringIO()
    sortby = SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    ps.dump_stats(outputMTLDirBase + '/InitializeAltMTLParallel.prof')
    print(s.getvalue())