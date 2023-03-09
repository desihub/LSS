#!/global/common/software/desi/cori/desiconda/20211217-2.0.0/conda/bin/python -u
from LSS.SV3.altmtltools import writeBitweights
from desiutil.log import get_logger
from LSS.bitweights import pack_bitweights
from sys import argv
import numpy as np
import argparse
import os
import multiprocessing as mp
from multiprocessing import Pool
import logging
import atexit
log = get_logger()
parser = argparse.ArgumentParser(
                    prog = 'MakeBitweights',
                    description = 'Convert a set of {ndir} realizations of alternate MTLs into bitweights and (optionally) probobs')
parser.add_argument('-o', '--outdir', dest='outdir', required=True, type=str, help = 'base output directory.')
parser.add_argument('-obscon', '--obscon', dest='obscon', default='DARK', help = 'observation conditions, either BRIGHT or DARK.', required = False, type = str)
parser.add_argument('-s', '--survey', dest='survey', default='sv3', help = 'DESI survey to create Alt MTLs for. Either sv3 or main.', required = False, type = str)
parser.add_argument('-n', '--ndir', dest='ndir', default=128, help = 'Random seed to ensure reproducability.', required = False, type = int)
parser.add_argument('-ppn', '--ProcPerNode', dest='ProcPerNode', default=None, help = 'Number of processes to spawn per requested node. If not specified, determined automatically from NERSC_HOST.', required = False, type = int)
parser.add_argument('-hpfn', '--HPListFile', dest='HPListFile', default=None, help = 'Name of a text file consisting only of one line of comma separated healpixel numbers for which the code will generate alt MTLs. If not specified, it will be automatically determined from the survey name.', required = False, type = str)
parser.add_argument('-ow', '--overwrite', dest = 'overwrite', default=False, action='store_true', help = 'pass this flag to regenerate already existing alt MTLs.')
parser.add_argument('-v', '--verbose', dest = 'verbose', default=False, action='store_true', help = 'set flag to enter verbose mode')
parser.add_argument('-d', '--debug', dest = 'debug', default=False, action='store_true', help = 'set flag to enter debug mode.')
'''
survey = argv[1]
obscon = argv[2]
ndirs = int(argv[3])
ProcPerNode = int(argv[4])
HPListFile = argv[5]
outdir = argv[6]
overwrite = bool(int(argv[7]))
'''
args = parser.parse_args()

mtlBaseDir = args.outdir + '/Univ{0:03d}/'
HPList = np.array(open(args.HPListFile,'r').readlines()[0].split(',')).astype(int)
print(HPList)

#mtlBaseDir = '/global/cscratch1/sd/jlasker/TestGeneralizedAltMTLScripts/alt_mtls_64dirs/Univ{0:03d}/'
#outdir = '/global/cscratch1/sd/jlasker/TestGeneralizedAltMTLScripts/alt_mtls_64dirs/'
#bw = makeBitweights(mtlBaseDir, ndirs = 64, hplist = hplist, debug = False)
#writeBitweights(mtlBaseDir, ndirs = 128, hplist = sv3dark, debug = False, outdir = outdir, survey = 'sv3', obscon = 'dark', allFiles = True)
#writeBitweights(mtlBaseDir, ndirs = 128, hplist = sv3dark, debug = False, outdir = outdir, survey = 'sv3', obscon = 'bright', allFiles = True)
#writeBitweights(mtlBaseDir, ndirs = None, hplist = None, debug = False, outdir = None, obscon = "dark", survey = 'sv3', overwrite = False, allFiles = False, splitByReal = False, splitNChunks = None)
def procFunc(nproc):

    thisHPList = np.array_split(HPList, args.ProcPerNode)[nproc]

    for hp in thisHPList:
        writeBitweights(mtlBaseDir, ndirs = args.ndir, hplist = [hp], debug = args.debug, verbose = args.verbose, outdir = args.outdir, survey = args.survey, obscon = args.obscon.lower(), allFiles = False, overwrite = args.overwrite)
#if survey.lower() == 'main':
#    for hp in HPList:
        
#else:
#    writeBitweights(mtlBaseDir, ndirs = args.ndir, hplist = HPList, debug = args.debug, verbose = args.verbose, outdir = outdir, survey = survey, obscon = obscon.lower(), allFiles = True, overwrite = overwrite)

#writeBitweights(mtlBaseDir, ndirs = ndir, hplist = HPList, debug = False, outdir = outdir, survey = survey, obscon = obscon.lower(), allFiles = True, overwrite = overwrite, splitByReal = splitByReal, splitNChunks = splitNChunks)
try:
    NNodes = int(os.getenv('SLURM_JOB_NUM_NODES'))
except:
    log.warning('no SLURM_JOB_NUM_NODES env set. You may not be on a compute node.')
    NNodes = 1
NodeID = int(os.getenv('SLURM_NODEID'))
SlurmNProcs = int(os.getenv('SLURM_NPROCS'))
NProc = int(NNodes*args.ProcPerNode)
if args.verbose or args.debug:
    log.debug('requested number of nodes: {0:d}'.format(NNodes))
    log.debug('requested number of directories/realizations: {0:d}'.format(args.ndir))
    log.debug('requested number of processes: {0:d}'.format(NProc))
inds = []
start = int(NodeID*NProc/SlurmNProcs)
end = int((NodeID + 1)*NProc/SlurmNProcs)
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

