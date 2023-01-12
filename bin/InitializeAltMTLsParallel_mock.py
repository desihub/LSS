#!/global/common/software/desi/cori/desiconda/20211217-2.0.0/conda/bin/python -u
from astropy.table import Table
import astropy
import multiprocessing as mp
from multiprocessing import Pool
import logging
import atexit
import desitarget.io as io
import glob
from LSS.SV3.mockaltmtltools import initializeAlternateMTLs
import numpy as np
import os 
from sys import argv
from desiutil.log import get_logger
log = get_logger()

#List of healpixels for SV3 (dark and bright are same)

print(argv)
seed = int(argv[1])
ndir = int(argv[2])
try:
    overwrite = bool(int(argv[3]))
except:
    raise ValueError('Invalid non-integer value of overwrite: {0}'.format(overwrite))
obscon = argv[4]
survey = argv[5]
outputMTLDirBase = argv[6]
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

with open(outputMTLDirBase + 'SeedFile', 'w') as f:
    f.write(str(seed))

HPList = np.array(open(HPListFile,'r').readlines()[0].split(',')).astype(int)
print(HPList)
NodeID = int(os.getenv('SLURM_NODEID'))
SlurmNProcs = int(os.getenv('SLURM_NPROCS'))

NProc = int(NNodes*64)





#outputMTLDirBase = "/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/altmtl/debug_jl/alt_mtls_mainTest_{0}dirs/".format(ndir)
outputMTLDir = outputMTLDirBase + "Univ{0:03d}/"


HPList = np.array(open(HPListFile,'r').readlines()[0].split(',')).astype(int)
print(HPList)

# If folder doesn't exist, then create it.
if not os.path.isdir(outputMTLDirBase):
    os.makedirs(outputMTLDirBase)

with open(outputMTLDirBase + 'SeedFile', 'w') as f:
    f.write(str(seed))

def procFunc(nproc):
    print('starting fxn call')
    for hpnum in HPList:
        if 'sv' in survey.lower():
            mtlprestr = survey.lower()
        else:
            mtlprestr = ''
        exampleledger = exampleledgerbase + '/{0}/{2}/{3}mtl-{2}-hp-{1}.ecsv'.format(survey.lower(),hpnum, obscon.lower(), mtlprestr)
        if os.path.isfile(exampleledger):
            initializeAlternateMTLs(exampleledger, outputMTLDir, genSubset = nproc, seed = seed, obscon = obscon, survey = survey, saveBackup = True, hpnum = hpnum, overwrite = overwrite)
        else:
            print(hpnum, 'not present')
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
p = Pool(NProc)
atexit.register(p.close)
result = p.map(procFunc,inds)
print('c')








