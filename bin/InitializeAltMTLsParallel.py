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
profile = True
pr = cProfile.Profile()

log = get_logger()

print('in python script REMOVE BEFORE PUSHING')
log.info('in python script REMOVE BEFORE PUSHING')
log.info(argv)
seed = int(argv[1])
ndir = int(argv[2])
try:
    overwrite = bool(int(argv[3]))
except:
    raise ValueError('Invalid non-integer value of overwrite: {0}'.format(overwrite))
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

NProc = int(NNodes*64)

log.info('requested number of nodes: {0:d}'.format(NNodes))
log.info('requested number of processes: {0:d}'.format(ndir))



#outputMTLDirBase = "/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_mainTest_{0}dirs/".format(ndir)
outputMTLDir = outputMTLDirBase + "Univ{0:03d}/"


HPList = np.array(open(HPListFile,'r').readlines()[0].split(',')).astype(int)
print(HPList)


def procFunc(nproc):
    log.info('starting fxn call')
    if 'sv' in survey.lower():
        log.info('sv survey')
        mtlprestr = survey.lower()
    else:
        log.info('non sv survey')
        mtlprestr = ''

    if os.path.exists(outputMTLDir + '/{0}/{2}/{3}mtl-{2}-hp-{1}.ecsv'.format(survey.lower(),HPList[-1], obscon.lower(), mtlprestr)):
        log.info('pathname')
        log.info(outputMTLDir + '/{0}/{2}/{3}mtl-{2}-hp-{1}.ecsv'.format(survey.lower(),HPList[-1], obscon.lower(), mtlprestr))
        return 42
    log.info('still going')
    for hpnum in HPList:
        exampleledger = exampleledgerbase + '/{0}/{2}/{3}mtl-{2}-hp-{1}.ecsv'.format(survey.lower(),hpnum, obscon.lower(), mtlprestr)
        initializeAlternateMTLs(exampleledger, outputMTLDir, genSubset = nproc, seed = seed, obscon = obscon, survey = survey, saveBackup = True, hpnum = hpnum, overwrite = overwrite, profile = profile)
    print('ending function call')
    return 42           

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
    print('i')
    print(i)
    inds.append(i)
    

NProc = len(inds)
assert(len(inds))
    
print('b')
print(inds)
print('running on NProc = {0} processes'.format(NProc))
p = Pool(NProc)
atexit.register(p.close)
log.info('running procFunc now on inds:')
log.info(inds)
pr.enable()
result = p.map(procFunc,inds)
pr.disable()
print('c')
s = io.StringIO()
sortby = SortKey.CUMULATIVE
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
ps.dump_stats(outputMTLDirBase + '/InitializeAltMTLParallel.prof')
print(s.getvalue())


