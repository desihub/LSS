'''
gather single exposure redshift info across all observations for a given target type
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

#sys.path.append('../py')

#from this package
import LSS.zcomp.zinfo as zi


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--release", help="what spectro release to use, e.g. blanc or daily",default='blanc') #eventually remove this and just gather everything
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
args = parser.parse_args()
print(args)

type = args.type
basedir = args.basedir
release = args.release
version = args.version



from desitarget.sv1 import sv1_targetmask

tarbit = -1
if type != 'All':
    tarbit = int(np.log2(sv1_targetmask.desi_mask[type]))

print('gathering all tile data for type '+type +' in '+release)

tp = 'SV1_DESI_TARGET'
print('targeting bit,  target program type; CHECK THEY ARE CORRECT!')
print(tarbit,tp)



#outputs
#basedir for official catalogs'/global/cfs/cdirs/desi/survey/catalogs
svdir = basedir+'/SV1/'

dirout = svdir+'redshift_comps/'+release+'/'+version+'/'+type

if not os.path.exists(svdir):
    os.mkdir(svdir)
    print('made '+svdir+' directory')

if not os.path.exists(svdir+'redshift_comps'):
    os.mkdir(svdir+'redshift_comps')
    print('made '+svdir+'redshift_comps directory')

if not os.path.exists(svdir+'redshift_comps/'+release):
    os.mkdir(svdir+'redshift_comps/'+release)
    print('made '+svdir+'redshift_comps/'+release+' directory')

if not os.path.exists(svdir+'redshift_comps/'+release+'/'+version):
    os.mkdir(svdir+'redshift_comps/'+release+'/'+version)
    print('made '+svdir+'redshift_comps/'+release+'/'+version+' directory')

if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)

if not os.path.exists(svdir+'/redshift_comps/logs'):
    os.mkdir(svdir+'/redshift_comps/logs')
    print('made '+svdir+'/redshift_comps/logs')

#set up log file
logfn = svdir+'/redshift_comps/logs/log'+datetime.now().isoformat()+'.txt'
logf = open(logfn,'w')
print('a log of what was run is going to '+logfn)

logf.write('running gatherSV_zinfo_1exps.py from '+os.getcwd()+'\n\n')
logf.write('arguments were:\n')
logf.write(str(args)+'\n')

  
expf = '/global/cfs/cdirs/desi/survey/observations/SV1/sv1-exposures.fits'  
exposures = fitsio.read(expf) #this will be used in depth calculations  
gt = ['BGS+MWS', 'ELG', 'QSO+ELG', 'QSO+LRG']
#location of inputs
tiledir = '/global/cfs/cdirs/desi/spectro/redux/'+release+'/tiles'

tiles = np.unique(exposures['TILEID'])
print('looking for data in these tiles:')
print(tiles)

tilew = []
for tile in tiles:
    tt = np.unique(exposures['TARGETS'][exposures['TILEID']==tile])[0]
    if np.isin(tt,gt): #that tile used cmx target bits
        tile = str(tile)
        coaddir = '/global/cfs/cdirs/desi/spectro/redux/'+release+'/tiles/'+tile+'/exposures/'
        print('going through tile '+tile)
        
        outf = dirout +'/'+tile+'_'+type+'zinfo'#_1exp.fits'
        try:
            fitsio.FITS(outf+'_1exp.fits')
            print(outf+' exists already')
            tilew.append(tile)
        except:
            a = zi.comb_exps_vert(tarbit,tp,tile,coaddir,exposures,outf,dirout)
            
            logf.write('compiled data for tile '+str(tile)+' written to '+outf+'\n')
            if a:
                tilew.append(tile)
        else:
            print('did not find data in '+release +' for tile '+tile)    

#combine all the tiles

dt = Table.read(dirout +'/'+tilew[0]+'_'+type+'zinfo_1exp.fits')
for i in range(1,len(tilew)):
    dtn = Table.read(dirout +'/'+tilew[i]+'_'+type+'zinfo_1exp.fits')
    dt = vstack([dt,dtn])

dt.sort('TARGETID')
outfall = dirout +'/alltiles_'+type+'zinfo_1exp.fits'
dt.write(outfall,format='fits', overwrite=True) 
print('wrote to '+outfall)
logf.write('combined all tiles, written to '+outfall)