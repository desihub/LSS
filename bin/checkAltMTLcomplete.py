import os
import numpy as np
import argparse
from astropy.table import Table

parser = argparse.ArgumentParser(prog = 'checkAltMTLcomplete', description = 'Checks that all realizations in an alt mtl directory are complete (all tiletracker entries done)')
parser.add_argument('-a', '--altMTLBaseDir', dest='altMTLBaseDir', required=True, type=str, help = 'the path to the location where alt MTLs are stored, up to, but not including survey and obscon information.')
parser.add_argument('-o', '--obscon', dest='obscon', default='DARK', help = 'observation conditions, either BRIGHT or DARK.', required = False, type = str)
parser.add_argument('--nreal', required=False, type=int, help='Number of processes. If running on a single node, this is the number of alt mtl realizations', default = 128)

args = parser.parse_args()

base_dir = os.path.join(args.altMTLBaseDir,'Univ{:03d}')

num_complete = 0

for i in range(128):
    
    tt = np.array(Table.read(os.path.join(base_dir.format(i),'mainsurvey-{}obscon-TileTracker.ecsv'.format(args.obscon))))

    ledgers = os.listdir('/pscratch/sd/d/desica/Y3Run2{}/Univ000/main/{}/'.format(args.obscon,args.obscon.lower()))

    if not np.all(tt['DONEFLAG']):
        print('Univ{:03d}: {}/{} entries are done'.format(i,np.sum(tt['DONEFLAG'])))
    else:
        num_complete += 1

if num_complete == 128:
    print('All realizations are complete')