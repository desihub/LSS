#!/global/common/software/desi/cori/desiconda/20211217-2.0.0/conda/bin/python -u
from multiprocessing import Pool
from LSS.SV3 import mockaltmtltools as amt
#import altmtltools as amt
import dill
from sys import argv
import os
from astropy.table import Table
import numpy as np
import astropy
import multiprocessing as mp
import logging
import atexit
import glob

print('argv')
print(argv)

#Base directory for the alternate MTLs created in the InitializeAltMTLs script
altmtlbasedir=argv[3]
#Include secondary targets?
secondary = bool(int(argv[4]))
#Observing conditions to process the observations
obscon = argv[5].lower()
#Survey whose observations you are processing
survey = argv[6]
numobs_from_ledger = bool(int(argv[7]))
#Force redo fiber assignment if it has already been done. 
redoFA = bool(int(argv[8]))

#Get information about environment for multiprocessing
NodeID = int(os.getenv('SLURM_NODEID'))
SlurmNProcs = int(os.getenv('SLURM_NPROCS'))

NNodes = int(argv[1])
NProc = int(NNodes*64)
print('NProc')
print(NProc)
print('NNodes')
print(NNodes)

try:
    try:
        print('argv[2]')
        print(argv[2])
        QRArg = int(argv[2])
    except:
        print('argv[2]v2')
        print(argv[2])
        QRArg = str(argv[2])
    print('QR Arg')
    print(QRArg)
    quickRestart = bool(QRArg)
except:

    print('No Quick Restart Argument given (or incorrect argument). Defaulting to false.')
    quickRestart = False

print('quick Restart')
print(quickRestart)
mtldir = '/global/cscratch1/sd/acarnero/mtl_test/init_mock000'
### AUREmtldir = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/'
zcatdir = '/global/cfs/cdirs/desi/spectro/redux/daily/'
ndirs = None
getosubp = False
multiproc = True
singleDate = True


def procFunc(nproc):
    print('starting fxn call')
    amt.loop_alt_ledger(obscon, survey = survey, mtldir = mtldir, zcatdir = zcatdir, altmtlbasedir = altmtlbasedir, ndirs = ndirs, numobs_from_ledger = numobs_from_ledger,secondary = secondary, getosubp = getosubp, quickRestart = quickRestart, multiproc = multiproc, nproc = nproc, singleDate = singleDate, redoFA = redoFA)
    print('ending function call')
    return 42           
#amt.quickRestartFxn(ndirs = 1, altmtlbasedir = altmtlbasedir, survey = 'sv3', obscon = 'dark')

#amt.loop_alt_ledger('dark', survey = survey, mtldir = mtldir, zcatdir = zcatdir, altmtlbasedir = altmtlbasedir, ndirs = 1, numobs_from_ledger = True,secondary = False, getosubp = False, quickRestart = True, multiproc = True, nproc = 0, redoFA = True, singleDate = '20210406')
#for d in zdates:
#    print('zdate')#
#    print(d)
#    amt.loop_alt_ledger('dark', survey = survey, mtldir = mtldir, zcatdir = zcatdir, altmtlbasedir = altmtlbasedir, ndirs = 1, numobs_from_ledger = True,secondary = False, getosubp = False, quickRestart = False, multiproc = True, nproc = 0, redoFA = True, singleDate = d)

#procFunc(0)

inds = []
start = int(NodeID*NProc/SlurmNProcs)
end = int((NodeID + 1)*NProc/SlurmNProcs)
print("NodeID")
print(NodeID)
print('start')
print(start)
print('end')
print(end)

for i in range(start, end):
    print('i')
    print(i)
    files = glob.glob(altmtlbasedir + "Univ{0:03d}/*".format(i))
    if len(files):
        pass
    else:
        print('no files in dir number {0}'.format(i))
        continue
    inds.append(i)
    
assert(len(inds))
    
print('b')
print(inds)
p = Pool(NProc)
atexit.register(p.close)
result = p.map(procFunc,inds)
print('c')
