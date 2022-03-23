'''
gather redshift info across all observations for a given target type
this just works on Rongpu's Cascades reruns
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
parser.add_argument("--release", help="either 3x_depth or 4x_depth",default='3x_depth') #eventually remove this and just gather everything
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

logf.write('running gatherSV_zinfo_alltiles_RZx.py from '+os.getcwd()+'\n\n')
logf.write('arguments were:\n')
logf.write(str(args)+'\n')

  
expf = '/global/cfs/cdirs/desi/survey/observations/SV1/sv1-exposures.fits'  
exposures = fitsio.read(expf) #this will be used in depth calculations  
gt = ['BGS+MWS', 'ELG', 'QSO+ELG', 'QSO+LRG','BACKUP','SSV']
#location of inputs
tiledir = '/global/cfs/cdirs/desi/users/rongpu/redux/cascades/'+release+'/'

tiles = np.unique(exposures['TILEID'])
print('looking for data in these tiles:')
print(tiles)

mfn = svdir+'/redshift_comps/logs/missingexposures.txt'
fo = open(svdir+'/redshift_comps/logs/missingexposures.txt','w')
fo.close()

tilew = []
for tile in tiles:
    tt = np.unique(exposures['TARGETS'][exposures['TILEID']==tile])[0]
    if np.isin(tt,gt): #make sure tiles are type we want
        tile = str(tile)
        #if tile != '80607':
        coaddir = '/global/cfs/cdirs/desi/users/rongpu/redux/cascades/'+release+'/'+tile
        #subsets = [x[0][len(coaddir):].strip('/') for x in os.walk(coaddir)] #something must work better than this, but for now...
        a = glob.glob(coaddir+'/*.fits')
        b = []
        for i in range(0,len(a)):
            b.append(a[i][-13:-5])        
        subsets = np.unique(b)
        print(subsets)

        if len(subsets) > 0:
            #print(subsets)
            print('going through tile '+tile)
            outf = dirout +'/'+tile+'_'+type+'zinfo.fits'
            if os.path.isfile(outf): 
                print(outf+' exists already')
                tilew.append(tile)

            else:           
                a = zi.comb_subset_vert(tarbit,tp,subsets,tile,coaddir,exposures,outf,tt,mfn=mfn,md='RZ')
                logf.write('compiled data for tile '+str(tile)+' written to '+outf+'\n')
                if a:
                    tilew.append(tile)

        else:
            print('did not find data in '+release +' for tile '+tile)    
      

#combine all the tiles

dt = Table.read(dirout +'/'+tilew[0]+'_'+type+'zinfo.fits')
dt['TILEID'] = int(tilew[0])
for i in range(1,len(tilew)):
    dtn = Table.read(dirout +'/'+tilew[i]+'_'+type+'zinfo.fits')
    dtn['TILEID'] = int(tilew[i])
    dt = vstack([dt,dtn], metadata_conflicts='silent')
    print(tilew[i],len(dt))

dt.sort('TARGETID')
col2remove = ['NUMEXP','NUMTILE','LAMBDA_REF','OBJTYPE','NUMTARGET','FIBERFLUX_IVAR_G','FIBERFLUX_IVAR_R','FIBERFLUX_IVAR_Z','DESI_TARGET','BGS_TARGET','MWS_TARGET','HPXPIXEL','NUM_TILEID','NUM_FIBER']
for col in col2remove:
    try:
        dt.remove_columns([col])
    except:
        print('didnt fine column to remove '+col)


outfall = dirout +'/alltiles_'+type+'zinfo.fits'
dt.write(outfall,format='fits', overwrite=True) 
print('wrote to '+outfall)
logf.write('combined all tiles, written to '+outfall)

# types = ['QSO','ELG','LRG','BGS_ANY']#,'MWS']
# 
# tiles = {'LRG':[80605,80609],'ELG':[80606,80608],'QSO':[80605,80607,80609],'BGS_ANY':[80613]}
# dates = {'LRG':[210224,21030],'ELG':[210218,210208],'QSO':[210223,210214,210210],'BGS_ANY':[210202]}
# 
# 
# dirvi = '/global/cfs/cdirs/desi/sv/vi/TruthTables/Blanc/'
# svdir = basedir+'/SV1/'
# dirz = svdir+'redshift_comps/'+release+'/'+version+'/'

#for i in range(0,len(types)):
#    tp =types[i]
# tp = type
# tilet = tiles[tp]
# datet = dates[tp]
# gt = []
# for it in range(0,len(tilet)):
#     date = str(datet[it])
#     tile = str(tilet[it])
#     tt=Table.read(dirvi+tp[:3]+'/'+'desi-vi_'+tp[:3]+'_tile'+tile+'_nightdeep_merged_all_'+date+'.csv',format='pandas.csv')
#     tt.keep_columns(['TARGETID','best_z','best_quality','best_spectype','all_VI_issues','all_VI_comments','merger_comment','N_VI'])
#     try:
#         tz = Table.read(dirz+'/'+tp+'/'+tile+'_'+tp+'zinfo.fits')
#         tj = join(tz,tt,join_type='left',keys='TARGETID')
#         tj['N_VI'].fill_value = 0
#         tj['N_VI'] = tj['N_VI'].filled() #should easily be able to select rows with N_VI > 0 to get desired info
#         tj['TILEID'] = tile
#         tj.write(dirz+'/'+tp+'/'+tile+'_'+tp+'zinfo_wVI.fits',format='fits',overwrite=True)
#         print('wrote file with VI info to '+dirz+'/'+tp+'/'+tile+'_'+tp+'zinfo_wVI.fits')
#         gt.append(tile)
#     except:
#         print('didnt find data for tile '+tile) 
# print(gt)
# #if len(tilet) > 1:
# dt = Table.read(dirz+'/'+tp+'/'+str(gt[0])+'_'+tp+'zinfo_wVI.fits')
# for it in range(1,len(gt)):
#     dtn = Table.read(dirz+'/'+tp+'/'+str(gt[it])+'_'+tp+'zinfo_wVI.fits')
#     dt = vstack([dt,dtn])
# 
# print(np.unique(dt['TILEID']))
# cols = ['z','zwarn','chi2','deltachi2','spectype','subtype']
# for i in range(1,5):
#     
#     dt['z_'+str(i)]=np.zeros(len(dt))
#     dt['zwarn_'+str(i)]=np.zeros(len(dt))
#     dt['chi2_'+str(i)]=np.zeros(len(dt))
#     dt['deltachi2_'+str(i)]=np.zeros(len(dt))
#     dt['spectype_'+str(i)] = 'GALAXY'
#     dt['subtype_'+str(i)] = 'GALAXY'
# for ii in range(0,len(dt)):
#     ln = dt[ii]
#     zfitdir = tiledir+str(ln['TILEID'])
#     zfits = zi.get_zfits(ln['TILEID'],ln['PETAL_LOC'],ln['subset'],ln['TARGETID'],zfitdir)
#     for jj in range(1,5):
#         for col in cols:
#             dt[col+'_'+str(jj)][ii] = zfits[jj][col]
#     if ii%1000 == 0:
#         print(ii)
# 
# #dt.sort('TARGETID')
# outfall = dirz +'/'+tp+'/allVItiles_'+tp+'zinfo_wVI.fits'
# dt.write(outfall,format='fits', overwrite=True) 
# print('wrote to '+outfall)
#             
