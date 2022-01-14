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
from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desitarget import targetmask
from desimodel.footprint import is_point_in_desi
from multiprocessing import Pool



#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   
import LSS.mkCat_singletile.fa4lsscat as fa

parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--minr", help="minimum number for random files",default=0)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 18 are available (use parallel script for all)",default=18) 
parser.add_argument("--nproc",help="number of processors to use",default=32)

#parser.add_argument("--rann",help='the number for the input random file (0-17)',default=0)



args = parser.parse_args()
print(args)

type = args.type
basedir = args.basedir
rm = int(args.minr)
rx = int(args.maxr)


if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    pr = 'BRIGHT'
    pdir = 'bright'
else:
    pr = 'DARK'
    pdir = 'dark'

pd = pdir

mt = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
wd = mt['SURVEY'] == 'main'
wd &= mt['ZDONE'] == 'true'
wd &= mt['FAPRGRM'] == pdir
mtld = mt[wd]
print('found '+str(len(mtld))+' '+pdir+' time main survey tiles with zdone true')

tiles4comb = Table()
tiles4comb['TILEID'] = mtld['TILEID']
tiles4comb['ZDATE'] = mtld['LASTNIGHT']

#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/main/LSS/'




if not os.path.exists(maindir+'/logs'):
    os.mkdir(maindir+'/logs')
    print('made '+maindir+'/logs')

if not os.path.exists(maindir+'/LSScats'):
    os.mkdir(maindir+'/LSScats')
    print('made '+maindir+'/LSScats')




randir = maindir+'random'
#logf.write('using random files '+str(rm)+ ' through '+str(rx)+' (this is python, so max is not inclusive)\n')
for i in range(rm,rx):
    if not os.path.exists(maindir+'random'+str(i)):
        os.mkdir(maindir+'random'+str(i))
        print('made '+str(i)+' random directory')


#construct a table with the needed tile information
if len(mtld) > 0:
    tilel = []
    ral = []
    decl = []
    mtlt = []
    fal = []
    obsl = []
    pl = []
    fver = []
    fahal = []
    
    #for tile,pro in zip(mtld['TILEID'],mtld['PROGRAM']):
    for tile in mtld['TILEID']:
        ts = str(tile).zfill(6)
        try:
            fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
            tilel.append(tile)
            ral.append(fht['TILERA'])
            decl.append(fht['TILEDEC'])
            mtlt.append(fht['MTLTIME'])
            fal.append(fht['FA_RUN'])
            obsl.append(fht['OBSCON'])
            fav = fht['FA_VER']
            if np.isin(fav,['2.2.0.dev2811','2.3.0','2.3.0.dev2838']):#2.3.0 confirmed to work for these
                fver.append('2.3.0')
            else:
                fver.append(fav)    
            #try:
            #    faha = fht['FA_HA']
            #except:
            #    faha = 0
            #    print(tile,'no FA_HA in this tile header')        
            #pl.append(pro)
            pl.append(pr)
        except:
            print('failed to find and/or get info for tile '+ts)    
    ta = Table()
    ta['TILEID'] = tilel
    ta['RA'] = ral
    ta['DEC'] = decl
    ta['MTLTIME'] = mtlt
    ta['FA_RUN'] = fal
    ta['OBSCON'] = obsl
    ta['PROGRAM'] = pl
    #ta['FA_HA'] = fahal
    #ta['FA_VER'] = fver
    print(np.unique(fver,return_counts=True))
    #wfv = (np.array(fver) == faver)
    #mtld =  mtld[wfv]
    #ta = ta[wfv]
else:
    print('no done tiles in the MTL')

print(len(ta))

dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'  

for rann in range(rm,rx):
	fnr = dirrt+'/randoms-1-'+str(rann)+'.fits'
	rt = fitsio.read(fnr,columns=['RA','DEC','TARGETID','MASKBITS','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z'])
	print('read input random file '+fnr)

	def doran(ii):

		nti = int(len(ta)/rx)+1
		tim = nti*ii
		tix = nti*(ii+1)
		if tix < len(ta):
			tiles = ta[tim:tix]
		else:
			tiles = ta[tim:]
		ct.randomtiles_main_fromran(tiles,rt,rann )

	N = args.nproc
	p = Pool(N)
	inds = []
	for i in range(0,N):
		inds.append(i)
	p.map(doran,inds)

            
     

#if __name__ == '__main__':
	#N = int(sys.argv[2])
