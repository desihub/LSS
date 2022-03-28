#!/global/common/software/desi/cori/desiconda/20200801-1.4.0-spec/conda/bin/python -u 
from desiutil.log import get_logger
from LSS.SV3.altmtltools import writeBitweights
from LSS.bitweights import pack_bitweights
from sys import argv
import numpy as np
log = get_logger()

survey = argv[1]
obscon = argv[2]
ndirs = int(argv[3])
splitByReal = bool(int(argv[4]))
splitNChunks = int(argv[5])
HPListFile = argv[6]
outdir = argv[7]
overwrite = argv[8]
mtlBaseDir = outdir + '/Univ{0:03d}/'
HPList = np.array(open(HPListFile,'r').readlines()[0].split(',')).astype(int)
print(HPList)

#mtlBaseDir = '/global/cscratch1/sd/jlasker/TestGeneralizedAltMTLScripts/alt_mtls_64dirs/Univ{0:03d}/'
#outdir = '/global/cscratch1/sd/jlasker/TestGeneralizedAltMTLScripts/alt_mtls_64dirs/'
#bw = makeBitweights(mtlBaseDir, ndirs = 64, hplist = hplist, debug = False)
#writeBitweights(mtlBaseDir, ndirs = 128, hplist = sv3dark, debug = False, outdir = outdir, survey = 'sv3', obscon = 'dark', allFiles = True)
#writeBitweights(mtlBaseDir, ndirs = 128, hplist = sv3dark, debug = False, outdir = outdir, survey = 'sv3', obscon = 'bright', allFiles = True)
#writeBitweights(mtlBaseDir, ndirs = None, hplist = None, debug = False, outdir = None, obscon = "dark", survey = 'sv3', overwrite = False, allFiles = False, splitByReal = False, splitNChunks = None)
writeBitweights(mtlBaseDir, ndirs = ndirs, hplist = HPList, debug = False, outdir = outdir, survey = survey, obscon = obscon.lower(), allFiles = True, overwrite = overwrite)
