#standard python
import sys
import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import fitsio
import glob
import argparse
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt
#from desihub
#from desitarget import targetmask
#from regressis, must be installed
#from regressis import DR9Footprint
#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.common_tools as common

from LSS.globals import main
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   

try:
	if os.environ['NERSC_HOST'] == 'cori':
		scratch = 'CSCRATCH'
	elif os.environ['NERSC_HOST'] == 'perlmutter':
		scratch = 'PSCRATCH'
	else:
		print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
		#sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 
except:
    print('no NERSC_HOST')
    scratch ='HOME'

parser = argparse.ArgumentParser()
parser.add_argument("--col_name", help="name of the column to normalize",default='WEIGHT_RF')
parser.add_argument("--tracer", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=os.environ[scratch])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--verspec",help="version for redshifts",default='iron')


args = parser.parse_args()
print(args)

basedir = args.basedir
version = args.version
specrel = args.verspec

maindir = basedir +'/'+args.survey+'/LSS/'
ldirspec = maindir+specrel+'/'

dirout = ldirspec+'LSScats/'+version+'/'
dirin = dirout

fname = dirin+args.tracer+'_full_HPmapcut.dat.fits'
indata = Table(fitsio.read(fname))

selgoodz = common.goodz_infull(args.tracer[:3],indata,zcol='Z_not4clus')

selngc = common.splitGC(indata)
seln = indata['PHOTSYS'] == 'N'

regl = ['N','NGCnotN','SGC']

if args.tracer[:3] == 'QSO':
    regl = ['N','Snotdes','des']
    zrl = [(0.8,1.3),(1.3,2.1),(2.1,3.5)]
    seldes = common.select_regressis_DES(indata)
if args.tracer[:3] == 'LRG':
    zrl = [(0.4,0.6),(0.6,0.8),(0.8,1.1)]
if args.tracer[:3] == 'ELG':
    zrl = [(0.8,1.1),(1.1,1.6)]
if args.tracer[:3] == 'BGS':
    zrl = [(0.1,0.4)]

for reg in regl:
    if reg == 'N':
        selreg = seln
    if reg == 'NGCnotN':
        selreg = selngc & ~seln
    if reg == 'SGC':
        selreg = ~selngc
    if reg == 'Snotdes':
        selreg = ~seln & ~seldes
    if reg == 'des':
        selreg = seldes
    for zr in zrl:
        selzr = indata['Z_not4clus'] > zr[0]
        selzr &= indata['Z_not4clus'] < zr[1]
        meanw = np.mean(indata[selgoodz&selreg&selzr][args.col_name])
        print('mean '+args.col_name+' was '+str(meanw)+' for '+str(zr)+ ' and '+reg)
        indata[args.col_name][selgoodz&selreg&selzr] /= meanw
        meanw = np.mean(indata[selgoodz&selreg&selzr][args.col_name])
        print('mean '+args.col_name+' now '+str(meanw)+' for '+str(zr)+ ' and '+reg+' (should be 1)')


common.write_LSS(indata,fname)

