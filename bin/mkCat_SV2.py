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

#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.SV2.cattools as ct
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   
import LSS.mkCat_singletile.fa4lsscat as fa

version = 'test'


mdir = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/sv2/bright/' #location of ledgers
tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.53.0/targets/sv2/resolve/bright/' #location of targets
mtld = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/mtl-done-tiles.ecsv') #log of tiles completed for mtl
imbits = [1,5,6,7,8,9,11,12,13]

sv2dir = '/global/cfs/cdirs/desi/survey/catalogs/SV2/LSS/'

from desitarget.sv2 import sv2_targetmask
type = 'BGS_ANY'
tarbit = int(np.log2(sv2_targetmask.desi_mask[type]))

if type == 'BGS_ANY':
    pr = 'BRIGHT'
else:
    pr = 'DARK'
wp = mtld['PROGRAM'] == pr
mtld = mtld[wp]


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
mkdtiles = False
combd = False
combr = False
mkfulldat = False
mkfullran = False
mkclusdat = False
mkclusran = True

if mktileran:
    ct.randomtiles_allSV2(ta,imin=rm,imax=rx)
    
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

if mkdtiles:
    for tile,zdate in zip(mtld['TILEID'],mtld['ZDATE']):
        ffd = dirout+type+str(tile)+'_full.dat.fits'
        tspec = ct.combspecdata(tile,zdate)
        pdict,goodloc = ct.goodlocdict(tspec)
        fbaf = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/0'+str(tile)[:2]+'/fiberassign-0'+str(tile)+'.fits.gz'
        wt = ta['TILEID'] == tile
        tars = read_targets_in_tiles(mdir,ta[wt],mtl=True)
        tars = inflate_ledger(tars,tdir)
        tars = tars[[b for b in list(tars.dtype.names) if b != 'Z']]
        tars = tars[[b for b in list(tars.dtype.names) if b != 'ZWARN']]
        tars = tars[[b for b in list(tars.dtype.names) if b != 'PRIORITY']]
        tars = join(tars,tspec,keys=['TARGETID'],join_type='left')
        tout = ct.gettarinfo_type(fbaf,tars,goodloc,tarbit,pdict)
        #tout = join(tfa,tspec,keys=['TARGETID','LOCATION'],join_type='left') #targetid should be enough, but all three are in both and should be the same
        print(tout.dtype.names)
        wz = tout['ZWARN']*0 == 0
        wzg = tout['ZWARN'] == 0
        print('there are '+str(len(tout[wz]))+' rows with spec obs redshifts and '+str(len(tout[wzg]))+' with zwarn=0')
        
        tout.write(ffd,format='fits', overwrite=True) 
        print('wrote matched targets/redshifts to '+ffd)
        #logf.write('made full data files\n')

if combd:
    print(len(mtld['TILEID']))
    ct.combtiles(mtld['TILEID'],dirout,type)    


if combr:
    print(len(mtld['TILEID']))
    for i in range(rm,rx):
        ct.combran(mtld,i,randir)
        
        
if mkfulldat:
	ct.mkfulldat(dirout+type+'Alltiles_full.dat.fits',imbits,tdir)
	#get_tilelocweight()
	#logf.write('ran get_tilelocweight\n')
	#print('ran get_tilelocweight\n')

if mkfullran:
    for ii in range(rm,rx):
        outf = dirout+type+'Alltiles_'+str(ii)+'_full.ran.fits'
        ct.mkfullran(randir,ii,imbits,outf)
    #logf.write('ran mkfullran\n')
    #print('ran mkfullran\n')

#needs to happen before randoms so randoms can get z and weights
if mkclusdat:
    ct.mkclusdat(dirout+type+'Alltiles_')
    #logf.write('ran mkclusdat\n')
    #print('ran mkclusdat\n')

if mkclusran:
    for ii in range(rm,rx):
        ct.mkclusran(dirout+type+'Alltiles_',ii)
    #logf.write('ran mkclusran\n')
    #print('ran mkclusran\n')
    
if mkNbar:
	e2e.mkNbar(target_type,program,P0=P0,omega_matter=omega_matter,truez=truez)
	logf.write('made nbar\n')
	print('made nbar\n')

if fillNZ:
	e2e.fillNZ(target_type,program,P0=P0,truez=truez)	
	logf.write('put NZ and weight_fkp into clustering catalogs\n')    
	print('put NZ and weight_fkp into clustering catalogs\n')
        