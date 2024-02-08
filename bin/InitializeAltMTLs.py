#!/global/common/software/desi/cori/desiconda/20211217-2.0.0/conda/bin/python -u
from desiutil.iers import freeze_iers
freeze_iers()

from astropy.table import Table
import desitarget.io as io
import glob
from LSS.SV3.altmtltools import initializeAlternateMTLs
import numpy as np
import os 
from sys import argv
from desiutil.log import get_logger
log = get_logger()

#List of healpixels for SV3 (dark and bright are same)

   
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

# If folder doesn't exist, then create it.
if not os.path.isdir(outputMTLDirBase):
    os.makedirs(outputMTLDirBase)

with open(outputMTLDirBase + 'SeedFile', 'w') as f:
    f.write(str(seed))

HPList = np.array(open(HPListFile,'r').readlines()[0].split(',')).astype(int)
print(HPList)

for i in HPList:
    print(i)
    if 'sv' in survey.lower():
        mtlprestr = survey.lower()
    else:
        mtlprestr = ''
    exampleledger = exampleledgerbase + '/{0}/{2}/{3}mtl-{2}-hp-{1}.ecsv'.format(survey.lower(),i, obscon.lower(), mtlprestr)
    #print(exampleledger)
    outputMTLDir = outputMTLDirBase + "Univ{0:03d}/"
    initializeAlternateMTLs(exampleledger, outputMTLDir, nAlt = ndir, seed = seed, obscon = obscon, survey = survey, saveBackup = True, hpnum = i,shuffleBrightPriorities = shuffleBrightPriorities, PromoteFracBGSFaint = PromoteFracBGSFaint, overwrite = overwrite)

