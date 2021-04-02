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

sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.SV2.cattools as ct
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   
import LSS.mkCat_singletile.fa4lsscat as fa


mdir = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/sv2/bright/' #location of ledgers
tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.53.0/targets/sv2/resolve/bright/' #location of targets
mtld = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/mtl-done-tiles.ecsv') #log of tiles completed for mtl

sv2dir = '/global/cfs/cdirs/desi/survey/catalogs/SV2/LSS/'

from desitarget.sv2 import sv2_targetmask
type = 'BGS_ANY'
tarbit = int(np.log2(sv2_targetmask.desi_mask[type]))


if not os.path.exists(sv2dir+'/logs'):
    os.mkdir(sv2dir+'/logs')
    print('made '+sv2dir+'/logs')

if not os.path.exists(sv2dir+'/LSScats'):
    os.mkdir(sv2dir+'/LSScats')
    print('made '+sv2dir+'/LSScats')

dirout = sv2dir+'LSScats/'+version+'/'
if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)


randir = sv2dir+'random'
rm = 0
rx = 18
#logf.write('using random files '+str(rm)+ ' through '+str(rx)+' (this is python, so max is not inclusive)\n')
for i in range(rm,rx):
    if not os.path.exists(sv2dir+'random'+str(i)):
        os.mkdir(sv2dir+'random'+str(i))
        print('made '+str(i)+' random directory')


#construct a table with the needed tile information
tilel = []
ral = []
decl = []
mtlt = []
fal = []
obsl = []
pl = []
for tile,pro in zip(mtld['TILEID'],mtld['PROGRAM']):
    fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/0'+str(tile)[:2]+'/fiberassign-0'+str(tile)+'.fits.gz')
    tilel.append(tile)
    ral.append(fht['TILERA'])
    decl.append(fht['TILEDEC'])
    mtlt.append(fht['MTLTIME'])
    fal.append(fht['FA_RUN'])
    obsl.append(fht['OBSCON'])
    pl.append(pro)
ta = Table()
ta['TILEID'] = tilel
ta['RA'] = ral
ta['DEC'] = decl
ta['MTLTIME'] = mtlt
ta['FA_RUN'] = fal
ta['OBSCON'] = obsl
ta['PROGRAM'] = pl

mktileran = False
runfa = False
mkfulld = True

if mktileran:
    ct.randomtiles_allSV2(ta)
    
if runfa:
    for ii in range(0,len(mtld)):
        tile = mtld['TILEID'][ii]
        fbah = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/0'+str(tile)[:2]+'/fiberassign-0'+str(tile)+'.fits.gz')
        dt = fbah['FA_RUN']
        ttemp = Table(ta[ii])
        ttemp['OBSCONDITIONS'] = 516
        ttemp['IN_DESI'] = 1
        ttemp.write('tiletemp.fits',format='fits', overwrite=True)
        for i in range(rm,rx):
            testfbaf = randir+str(i)+'/fba-0'+str(tile)+'.fits'
            if os.path.isfile(testfbaf):
                print('fba file already made')
            else:   
                
                fa.getfatiles(randir+str(i)+'/tilenofa-'+str(tile)+'.fits','tiletemp.fits',dirout=randir+str(i)+'/',dt = dt)

if mkfulld:
    for tile,zdate in zip(mtld['TILEID'],mtld['ZDATE'])
        ffd = dirout+type+str(tile)+'_full.dat.fits'
        tspec = ct.combspecdata(tile,zdate)
        pdict,goodloc = ct.goodlocdict(tspec)
        fbaf = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/0'+str(tile)[:2]+'/fiberassign-0'+str(tile)+'.fits.gz'
        wt = ta['TILEID'] == tile
        tars = read_targets_in_tiles(mdir,ta[wt],mtl=True)
        tfa = ct.gettarinfo_type(fbaf,tars,goodloc,tarbit,pdict,tp=tp)
        tout = join(tfa,tspec,keys=['TARGETID','LOCATION'],join_type='left') #targetid should be enough, but all three are in both and should be the same
        print(tout.dtype.names)
        wz = tout['ZWARN']*0 == 0
        wzg = tout['ZWARN'] == 0
        print('there are '+str(len(tout[wz]))+' rows with spec obs redshifts and '+str(len(tout[wzg]))+' with zwarn=0')
        
        tout.write(ffd,format='fits', overwrite=True) 
        print('wrote matched targets/redshifts to '+ffd)
        logf.write('made full data files\n')

    