'''
gather redshift info across all observations for a given target type; for now from a single tile
'''

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


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--tile", help="observed tile to use") #eventually remove this and just gather everything
args = parser.parse_args()
type = args.type
tile = args.tile


if type == 'LRG':
    tarbit = 0 #targeting bit
if type == 'QSO':
    tarbit = 2
if type == 'ELG':
    tarbit = 1

print('gathering type,tile')
print(type,tile)
tp = 'SV1_DESI_TARGET'
print('targeting bit,  target program type; CHECK THEY ARE CORRECT!')
print(tarbit,tp)

#location of inputs
coaddir = '/global/cscratch1/sd/rongpu/desi/sv1/single_exp_coadd_blanc/'+tile
#subsets = [x[0][len(coaddir):].strip('/') for x in os.walk(coaddir)] #something must work better than this, but for now...


#outputs
svdir = '/project/projectdirs/desi/users/ajross/catalogs/SV/'
version = 'testRP/'
dirout = svdir+'redshift_comps/'+version
outf = dirout +'/'+tile+'_'+type+'zinfo.fits'

if not os.path.exists(svdir+'redshift_comps'):
    os.mkdir(svdir+'redshift_comps')
    print('made '+svdir+'redshift_comps random directory')

if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)
  
expf = '/global/cfs/cdirs/desi/users/raichoor/fiberassign-sv1/sv1-exposures.fits'  
exposures = fitsio.read(expf) #this will be used in depth calculations   

fls = glob.glob(coaddir+'/zbest-*.fits') 

exps = []
for fl in fls:
    try:
        fitsio.read(fl,ext='FIBERMAP')
        ep = int(fl[-10:-5])
        exps.append(ep)
    except:
        print('no fibermap for exposure '+str(ep))    

ss = 0 #use to switch from creating to concatenating
for ep in exps:
	print('going through exposure '+str(ep))
	specs = []
	info = exposures[exposures['EXPID'] == ep]
	bd = info['B_DEPTH'][0]
	rd = info['R_DEPTH'][0]
	zd = info['Z_DEPTH'][0]
	tspec = Table.read(coaddir+'/zbest-000'+str(ep)+'.fits',hdu='ZBEST')  
	tf = Table.read(coaddir+'/zbest-000'+str(ep)+'.fits',hdu='FIBERMAP')
	tspec = join(tspec,tf,keys=['TARGETID'])
	tspec['B_DEPTH'] = bd
	tspec['R_DEPTH'] = rd  
	tspec['Z_DEPTH'] = zd  
	
	wtype = ((tspec[tp] & 2**tarbit) > 0)
	print(str(len(tspec))+' total entries '+str(len(tspec[wtype]))+' that are '+type+' entries with '+str(len(np.unique(tspec[wtype]['TARGETID'])))+' unique target IDs')
	tspec = tspec[wtype]
	tspec['subset'] = str(ep)
	if ss == 0:
		tspect = tspec
		ss = 1
	else:
		tspect = vstack([tspect,tspec])
	print('there are now '+str(len(tspect)) +' entries with '+str(len(np.unique(tspect['TARGETID'])))+' unique target IDs')    

#add in deep

fl = '/project/projectdirs/desi/users/ajross/catalogs/SV/redshift_comps/test/'+str(tile)+'_'+type+'zinfo.fits'
td = Table.read(fl)
td = td[td['subset']=='deep']

kc = []
for col in tspect.columns:
    kc.append(col)

td.keep_columns(kc)

tspect = vstack([tspect,td])


tspect.sort('TARGETID')
tspect.write(outf,format='fits', overwrite=True) 
print('wrote to '+outf)
        

