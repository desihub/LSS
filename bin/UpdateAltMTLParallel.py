#!/global/common/software/desi/perlmutter/desiconda/20240425-2.2.0/conda/bin/python -u
from desiutil.iers import freeze_iers
freeze_iers()

import os
import numpy as np
from astropy.table import Table
import argparse
from multiprocessing import Pool
from LSS.SV3 import altmtltools as amt

parser = argparse.ArgumentParser(prog = 'UpdateAltMTLParallel', description = 'Updates existing alt ledgers by generating a new tiletracker entry and running dateloop')
parser.add_argument('-a', '--altMTLBaseDir', dest='altMTLBaseDir', required=True, type=str, help = 'the path to the location where alt MTLs are stored, up to, but not including survey and obscon information.')
parser.add_argument('-d', '--endDate', required=True, type=int, help = 'The date up to which the alt mtl ledgers will be updated')
parser.add_argument('-o', '--obscon', dest='obscon', default='DARK', help = 'observation conditions, either BRIGHT or DARK.', required = False, type = str)
parser.add_argument('-s', '--survey', dest='survey', default='sv3', help = 'DESI survey to create Alt MTLs for. Either sv3 or main.', required = False, type = str)
parser.add_argument('--nproc', required=True, type=int, help='Number of processes. If running on a single node, this is the number of alt mtl realizations')


parser.add_argument('-rep', '--reproducing', action='store_true', dest='reproducing', default=False, help = 'WARNING: THIS FLAG SHOULD ONLY BE USED FOR DEBUGGING. Pass this flag to confirm to the alt mtl code that you are trying to reproduce real MTLs. This option should (must?) be used in conjunction with --shuffleSubpriorities.', required = False)
parser.add_argument('-mock', '--mock', dest = 'mock', default=False, action='store_true', help = 'set flag if running pipeline on mocks.')

args = parser.parse_args()
print(args)


def procFunc(nproc):
    altmtldir = os.path.join(args.altMTLBaseDir,'Univ{:03d}'.format(nproc))

    #reading tile tracker to check if all entries are complete (this could be made more lightweight certainly)
    tt = np.array(Table.read(os.path.join(altmtldir,'mainsurvey-{}obscon-TileTracker.ecsv'.format(args.obscon)))) 

    #Tile tracker is only updated when the existing tile tracker has all entries complete and the skip update flag is not set
    if not np.all(tt['DONEFLAG']):
        amt.updateTileTracker(altmtldir, args.endDate)

    amt.loop_alt_ledger(obscon = args.obscon, survey = args.survey, mtldir = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/', zcatdir = '/global/cfs/cdirs/desi/spectro/redux/daily/', altmtlbasedir = args.altMTLBaseDir, ndirs = None, numobs_from_ledger = True, secondary = False, getosubp = False, quickRestart = False, multiproc = True, nproc = nproc, singleDate = False, redoFA = False, mock = args.mock, targets = None, debug = False, verbose = False, reproducing = args.reproducing)

    return

#indices of univ to update
inds = range(args.nproc)
    
p = Pool()
result = p.map(procFunc,inds)