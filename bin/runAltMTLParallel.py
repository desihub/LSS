#!/global/common/software/desi/cori/desiconda/20211217-2.0.0/conda/bin/python -u
from multiprocessing import Pool
from LSS.SV3 import altmtltools as amt
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
pr = cProfile.Profile()
pr.enable()
log = get_logger()
print('argv')
print(argv)

#Base directory for the alternate MTLs created in the InitializeAltMTLs script
altmtlbasedir=argv[3]
#Include secondary targets?
assert((str(argv[4]) == '0') | (str(argv[4]) == '1'))
secondary = bool(int(argv[4]))
#Observing conditions to process the observations
obscon = argv[5].lower()
#Survey whose observations you are processing
survey = argv[6]
numobs_from_ledger = bool(int(argv[7]))
#Force redo fiber assignment if it has already been done. 
assert((str(argv[8]) == '0') | (str(argv[8]) == '1'))
redoFA = bool(int(argv[8]))
#getosubp: grab subpriorities from the original (exampleledgerbase) MTLs
#This should only be turned on for testing/debugging purposes
assert((str(argv[9]) == '0') | (str(argv[9]) == '1'))

getosubp = bool(int(argv[9]))
log.info('getosubp value = {0}'.format(getosubp))
#Get information about environment for multiprocessing
NodeID = int(os.getenv('SLURM_NODEID'))
SlurmNProcs = int(os.getenv('SLURM_NPROCS'))
NProcPerNode=32
NNodes = int(argv[1])
NProc = int(NNodes*NProcPerNode)
log.info('NProc = {0:d}'.format(NProc))
log.info('NNodes = {0:d}'.format(NNodes))


try:
    try:
        QRArg = int(argv[2])
    except:
        QRArg = str(argv[2])
    quickRestart = bool(QRArg)
except:

    log.info('No Quick Restart Argument given (or incorrect, non-boolean, non-integer, argument). Defaulting to false.')
    quickRestart = False

print('quick Restart')
print(quickRestart)
mtldir = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/'
zcatdir = '/global/cfs/cdirs/desi/spectro/redux/daily/'
ndirs = None
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



inds = []
start = int(NodeID*NProc/SlurmNProcs)
end = int((NodeID + 1)*NProc/SlurmNProcs)
log.info('NodeID = {0:d}'.format(NodeID))
log.info('StartProc = {0:d}'.format(start))
log.info('EndProc = {0:d}'.format(end))


for i in range(start, end):
    log.info('Process i = {0}'.format(i))
    files = glob.glob(altmtlbasedir + "Univ{0:03d}/*".format(i))
    if len(files):
        pass
    else:
        log.info('no files in dir number {0}, not processing that directory.'.format(i))
        continue
    inds.append(i)
    
assert(len(inds))
    
print('b')
print(inds)
p = Pool(NProc)
atexit.register(p.close)
result = p.map(procFunc,inds)
print('c')

pr.disable()
s = io.StringIO()
sortby = SortKey.CUMULATIVE
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
ps.dump_stats(altmtlbasedir + '/runAltMTLParallel.prof')
print(s.getvalue())