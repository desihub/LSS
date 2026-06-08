'''
gather redshift info across all observations for a given target type
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
parser.add_argument("--release", help="what spectro release to use, e.g. blanc or daily",default='cascades') #eventually remove this and just gather everything
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

logf.write('running gatherSV_zinfo_alltiles_inpar.py from '+os.getcwd()+'\n\n')
logf.write('arguments were:\n')
logf.write(str(args)+'\n')

  
expf = '/global/cfs/cdirs/desi/survey/observations/SV1/sv1-exposures.fits'  
exposures = fitsio.read(expf) #this will be used in depth calculations  
gt = ['BGS+MWS', 'ELG', 'QSO+ELG', 'QSO+LRG','BACKUP','SSV']
#location of inputs
tiledir = '/global/cfs/cdirs/desi/spectro/redux/'+release+'/tiles/'

tiles = np.unique(exposures['TILEID'])
print('looking for data in these tiles:')
print(tiles)

mfn = svdir+'/redshift_comps/logs/missingexposures.txt'
fo = open(svdir+'/redshift_comps/logs/missingexposures.txt','w')
fo.close()

tilew = []
#for tile in tiles:

def get_tilezinfo(tile):
    tt = np.unique(exposures['TARGETS'][exposures['TILEID']==tile])[0]
    if np.isin(tt,gt): #that tile used cmx target bits
        tile = str(tile)
        coaddir = '/global/cfs/cdirs/desi/spectro/redux/'+release+'/tiles/'+tile
        subsets = [x[0][len(coaddir):].strip('/') for x in os.walk(coaddir)] #something must work better than this, but for now...
        print(subsets)
        if len(subsets) > 1:
            #print(subsets)
            print('going through tile '+tile)
            outf = dirout +'/'+tile+'_'+type+'zinfo.fits'
            if os.path.isfile(outf): 
                print(outf+' exists already')
                #tilew.append(tile)
                a = True

            else:
                a = zi.comb_subset_vert(tarbit,tp,subsets,tile,coaddir,exposures,outf,tt,mfn=mfn)
                logf.write('compiled data for tile '+str(tile)+' written to '+outf+'\n')
                #if a:
                    #tilew.append(tile)
                #    return True
                #else:
                #    return False    

            if a:
                print('adding info from hd5 files')
                outfall = dirout +'/'+tile+'_'+type+'zinfo_wh5.fits'
                if os.path.isfile(outfall): 
                    print(outfall+' exists already')                
                else:
                    dt = Table.read(outf)
                    cols = ['z','zwarn','chi2','deltachi2','spectype','subtype']
                    for i in range(1,5):
    
                        dt['z_'+str(i)]=np.zeros(len(dt))
                        dt['zwarn_'+str(i)]=np.zeros(len(dt))
                        dt['chi2_'+str(i)]=np.zeros(len(dt))
                        dt['deltachi2_'+str(i)]=np.zeros(len(dt))
                        dt['spectype_'+str(i)] = 'GALAXY'
                        dt['subtype_'+str(i)] = 'GALAXY'
                    for ii in range(0,len(dt)):
                        ln = dt[ii]
           
                        #if ln['RZR'] != 'N':
                        #   zfitdir = '/global/cfs/cdirs/desi/users/rongpu/redux/cascades/'+ln['RZR']+'/'+str(ln['TILEID'])
                        #else:
                        zfitdir = tiledir+str(ln['TILEID'])+'/'+ln['subset']+'/'    
            
                        fl = zfitdir+'/redrock-'+str(ln['PETAL_LOC'])+'-'+str(ln['TILEID'])+'-'+ln['subset']+'.h5'
            
                        zfits = zi.get_zfits(fl,ln['TARGETID'])
                        for jj in range(1,5):
                            for col in cols:
                                dt[col+'_'+str(jj)][ii] = zfits[jj][col]
                
                    dt.write(outfall,format='fits', overwrite=True) 
                    print('wrote to '+outfall)
                return a


        else:
            print('did not find data in '+release +' for tile '+tile)   
            return False 
      
if __name__ == '__main__':
    from multiprocessing import Pool
    import sys
    #N = int(sys.argv[2])
    N = 32
    p = Pool(N)

    expf = '/global/cfs/cdirs/desi/survey/observations/SV1/sv1-exposures.fits'  
    exps = fitsio.read(expf)

    tiles = np.unique(exps['TILEID'])
    print('going through '+str(len(tiles))+' tiles')

    for j in range(0,len(tiles),N):
        #get_tilezinfo(tiles[j])
        inds = []
        for i in range(j,j+N):
            if i == len(tiles):
                break
            inds.append(tiles[i])
        p.map(get_tilezinfo,inds)


    #combine all the tiles

    dt = Table.read(dirout +'/'+str(tiles[0])+'_'+type+'zinfo_wh5.fits')
    dt['TILEID'] = int(tiles[0])
    for i in range(1,len(tiles)):
        tf = dirout +'/'+str(tiles[i])+'_'+type+'zinfo_wh5.fits'
        if os.path.isfile(tf):    
            dtn = Table.read(dirout +'/'+str(tiles[i])+'_'+type+'zinfo_wh5.fits')
            dtn['TILEID'] = int(tiles[i])
            dt = vstack([dt,dtn])
        else:
            print('did not find file for tile '+str(tiles[i]))  

    dt.sort('TARGETID')
    col2remove = ['NUMEXP','NUMTILE','LAMBDA_REF','OBJTYPE','NUMTARGET','FIBERFLUX_IVAR_G','FIBERFLUX_IVAR_R','FIBERFLUX_IVAR_Z','DESI_TARGET','BGS_TARGET','MWS_TARGET','HPXPIXEL','NUM_TILEID','NUM_FIBER']
    for col in col2remove:
        try:
            dt.remove_columns([col])
        except:
            print('didnt find column to remove '+col)
    outfall = dirout +'/alltiles_'+type+'zinfo_wh5.fits'
    dt.write(outfall,format='fits', overwrite=True) 
    print('wrote to '+outfall)
    logf.write('combined all tiles, written to '+outfall)
    