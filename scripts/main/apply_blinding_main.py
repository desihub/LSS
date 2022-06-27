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
#import LSS.main.cattools as ct
#import LSS.common_tools as common
#import LSS.imaging.select_samples as ss
#from LSS.globals import main
import LSS.blinding_tools as blind
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir_in", help="base directory for input, default is location for official catalogs",default='/global/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--basedir_out", help="base directory for output, default is C(P)SCRATCH",default=os.environ[scratch])
parser.add_argument("--version", help="catalog version",default='EDAbeta')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='DA02')
parser.add_argument("--verspec",help="version for redshifts",default='guadalupe')
parser.add_argument("--notqso",help="if y, do not include any qso targets",default='n')
parser.add_argument("--baoblind",help="if y, do the bao blinding shift",default='n')
parser.add_argument("--rsdblind",help="if y, do the bao blinding shift",default='n')


args = parser.parse_args()
print(args)

type = args.type
#basedir = args.basedir
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
maindir = args.basedir_in +'/'+args.survey+'/LSS/'

ldirspec = maindir+specrel+'/'

dirin = ldirspec+'LSScats/'+version+'/'

dirout = args.basedir_out+'LSScats/'+version+'/blinded/'

if not os.path.exists(args.basedir_out+'LSScats/'):
    os.mkdir(args.basedir_out+'LSScats/')
    print('made '+args.basedir_out+'LSScats/')    

if not os.path.exists(args.basedir_out+'LSScats/'+version):
    os.mkdir(args.basedir_out+'LSScats/'+version)
    print('made '+args.basedir_out+'LSScats/'+version)    


if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)    

w0 = -0.95
wa = 0.3

regl = ['_S','_N']
if args.baoblind == 'y':
	data = Table(fitsio.read(dirin+type+notqso+'_full.dat.fits'))
	outf = dirout + type+notqso+'_full.dat.fits'
	blind.apply_zshift_DE(data,outf,w0=w0,wa=wa,zcol='Z_not4clus')

if args.rsdblind == 'y':	
	for reg in regl:
		fnd = dirout+type+notqso+reg+'_clustering.dat.fits'
		fndr = dirout+type+notqso+reg+'_clustering.MGrsd.dat.fits'
		data = Table(fitsio.read(fnd))
		data_real = Table(fitsio.read(fndr))
		
		out_file = fnd
		blind.apply_zshift_RSD(data,data_real,out_file,fgrowth_fid=0.8,fgrowth_blind=0.9)
		
