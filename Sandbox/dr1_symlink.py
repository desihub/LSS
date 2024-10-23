import os
import glob
import sys
import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--survey",default='Y1')
parser.add_argument("--pubrelease",default='dr1')
parser.add_argument("--specrel",default='iron')
parser.add_argument("--outroot",default=None)
args = parser.parse_args()

if args.outroot is None:
    args.outroot = os.getenv('SCRATCH')

indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/'
pubdir = args.outroot+'/desi/survey/catalogs/'+args.pubrelease+'/'

print('root directory is '+pubdir)
versions = ['v1.2','v1.5','v1.5pip']
for ver in versions:
    verdir = pubdir+'/LSS/'+args.specrel+'/LSScats/'+ver
    if not os.path.exists(verdir):
        os.makedirs(verdir)

amtldir = pubdir+'/LSS/altmtl'
if not os.path.exists(amtldir):
    os.makedirs(amtldir)

fnames = ['LSS/collisions-BRIGHT.fits','LSS/collisions-DARK.fits','tiles-BRIGHT.fits','tiles-DARK.fits']
for fname in fnames:
    os.system('ln -s '+indir+fname+' ' +pubdir+fname)

os.system('ln -s '+indir+'/LSS/altmtl/JL_Y1Run2DARK '+pubdir+'/LSS/altmtl/DARK')
os.system('ln -s '+indir+'/LSS/altmtl/JL_Y1Run3BRIGHT '+pubdir+'/LSS/altmtl/BRIGHT')

progl = ['DARK','BRIGHT']
ranl = np.arange(0,18)

for rann in ranl:
    randir = pubdir+'/LSS/random'+str(rann)
    if not os.path.exists(randir):
        os.makedirs(randir)
    
    for prog in progl:
        fname = 'LSS/random'+str(rann)+'/pota-'+prog+'.fits'
        os.system('ln -s '+indir+fname+' ' +pubdir+fname)
        fname = 'LSS/'+args.specrel+'/rancomb_'+str(rann)+prog.lower()+'wdupspec_zdone.fits'
        os.system('ln -s '+indir+fname+' ' +pubdir+fname)

fname = 'LSS/'+args.specrel+'/unique_badfibers.txt'
os.system('ln -s '+indir+fname+' ' +pubdir+fname)

tracers = ['BGS_ANY','BGS_BRIGHT','LRG','ELG_LOPnotqso','QSO']
for tr in tracers:
    fname = 'LSS/'+args.specrel+'/datcomb_'+tr+'_tarspecwdup_zdone.fits'
    os.system('ln -s '+indir+fname+' ' +pubdir+fname)

for prog in progl:
    fname = 'LSS/'+args.specrel+'/datcomb_'+prog.lower()+'_spec_zdone.fits'
    os.system('ln -s '+indir+fname+' ' +pubdir+fname)
    fname = 'LSS/'+args.specrel+'/datcomb_'+prog.lower()+'_zmtl_zdone.fits'
    os.system('ln -s '+indir+fname+' ' +pubdir+fname)

fname = 'LSS/'+args.specrel+'/emlin_catalog.fits'
os.system('ln -s '+indir+fname+' ' +pubdir+fname)


for ver in versions:
	fnames = glob.glob(indir+'/LSS/'+args.specrel+'/LSScats/'+ver+'/*.fits')
	for fname in fnames:
		os.system('ln -s '+fname +' ' +pubdir+fname.replace(indir,''))
	fnames = glob.glob(indir+'/LSS/'+args.specrel+'/LSScats/'+ver+'/*nz.txt')
	for fname in fnames:
		os.system('ln -s '+fname +' ' +pubdir+fname.replace(indir,''))
	dirname = '/LSS/'+args.specrel+'/LSScats/'+ver+'/hpmaps'
	os.system('ln -s '+indir+dirname+' ' +pubdir+dirname)
	