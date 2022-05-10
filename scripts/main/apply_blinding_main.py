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
from desitarget import targetmask
#from regressis, must be installed
from regressis import DR9Footprint
#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
import LSS.common_tools as common
import LSS.imaging.select_samples as ss
from LSS.globals import main
import LSS.blinding_tools as blind
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='main')
parser.add_argument("--verspec",help="version for redshifts",default='everest')
parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')


args = parser.parse_args()
print(args)

type = args.type
basedir = args.basedir
version = args.version
specrel = args.verspec

notqso = ''
if args.notqso == 'y':
    notqso = 'notqso'

print('blinding catalogs for tracer type '+type+notqso)
    

if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    prog = 'BRIGHT'

else:
    prog = 'DARK'

progl = prog.lower()

#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/'+args.survey+'/LSS/'

ldirspec = maindir+specrel+'/'

dirin = ldirspec+'LSScats/'+version+'/'

dirout = ldirspec+'LSScats/'+version+'/blinded/'
if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)    

data = Table(fitsio.read(dirin+type+notqso+'zdone_full.dat.fits'))
outf = dirout + type+notqso+'zdone_full.dat.fits'
w0 = -0.95
wa = 0.3
blind.apply_zshift_DE(data,outf,w0=w0,wa=wa,zcol='Z_not4clus')