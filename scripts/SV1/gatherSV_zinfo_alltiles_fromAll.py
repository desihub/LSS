'''
get redshift info across all observations splitting them out by target type
then get the VI info
assumes gatherSV_zinfo_alltiles.py has been run for type all
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
#parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--release", help="e.g., cascades or blanc",default='cascades') #eventually remove this and just gather everything
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
args = parser.parse_args()
print(args)

#type = args.type
basedir = args.basedir
release = args.release
version = args.version



from desitarget.sv1 import sv1_targetmask

#outputs
#basedir for official catalogs'/global/cfs/cdirs/desi/survey/catalogs
svdir = basedir+'/SV1/'


#set up log file
logfn = svdir+'/redshift_comps/logs/log'+datetime.now().isoformat()+'.txt'
logf = open(logfn,'w')
print('a log of what was run is going to '+logfn)

logf.write('running gatherSV_zinfo_alltiles_fromAll.py from '+os.getcwd()+'\n\n')
logf.write('arguments were:\n')
logf.write(str(args)+'\n')


#tarbit = -1
#if type != 'All':
#    tarbit = int(np.log2(sv1_targetmask.desi_mask[type]))

#print('gathering all tile data for type '+type +' in '+release)

#tp = 'SV1_DESI_TARGET'
#print('targeting bit,  target program type; CHECK THEY ARE CORRECT!')
#print(tarbit,tp)

allrzf = svdir+'redshift_comps/'+release+'/'+version+'/All/alltiles_Allzinfo_wrz.fits'
if os.path.isfile(allrzf): 
    print('combined file '+allrzf + 'exists')
    tz = Table.read(allrzf)
else:    
    tz = Table.read(svdir+'redshift_comps/'+release+'/'+version+'/All/alltiles_Allzinfo.fits')
    sl = np.full(len(tz),'N',dtype='U20')
    tz['RZR'] = sl
    rzdirs = ['3x_depth','4x_depth','single_exposures']

    for rzdir in rzdirs:
        rzf = svdir+'redshift_comps/'+rzdir+'/'+version+'/All/alltiles_Allzinfo.fits'

        fz = Table.read(rzf)
        fz['RZR'] = rzdir
        #fzs = np.array(fz['subset'],dtype='U20')
        #fz['subset'] = np.core.defchararray.add(fzs,rzdir)

        tz = vstack([tz,fz])
        if rzdir == '3x_depth':
            print(np.unique(fz['RZR']))
            print(np.unique(tz['RZR']))
        print(len(tz))    

    tz.write(svdir+'redshift_comps/'+release+'/'+version+'/All/alltiles_Allzinfo_wrz.fits',overwrite=True,format='fits')

types = ['ELG','LRG','BGS_ANY','MWS_ANY','QSO']

vitiles = {'LRG':[80605,80609],'ELG':[80606,80608,80610],'QSO':[80605,80607,80609],'BGS_ANY':[80613]}
vidates = {'LRG':[210224,21030],'ELG':[210218,210208,210308],'QSO':[210223,210214,210210],'BGS_ANY':[210202]}
dirvi = '/global/cfs/cdirs/desi/sv/vi/TruthTables/Blanc/'

expf = '/global/cfs/cdirs/desi/survey/observations/SV1/sv1-exposures.fits'  
exposures = fitsio.read(expf) #this will be used in depth calculations  

tiledir = '/global/cfs/cdirs/desi/spectro/redux/'+release+'/tiles/'

for tp in types:

    dirout = svdir+'redshift_comps/'+release+'/'+version+'/'+tp

    if not os.path.exists(dirout):
        os.mkdir(dirout)
        print('made '+dirout)
        
    tarbit = int(np.log2(sv1_targetmask.desi_mask[tp]))
    
    wt = (tz['SV1_DESI_TARGET'] & sv1_targetmask.desi_mask[tp]) > 0
    dt = tz[wt]

    outfall = dirout +'/alltiles_'+tp+'zinfo.fits'
    dt.write(outfall,overwrite=True,format='fits')
    print(len(dt))
    #fitsio.write(outfall,dt,clobber=True)
    print('wrote '+outfall)
    
    if tp != 'MWS_ANY':
        gt = []
        tilet = vitiles[tp]
        datet = vidates[tp]
        ft = Table.read(dirout +'/alltiles_'+tp+'zinfo.fits')
        dirz = svdir+'redshift_comps/'+release+'/'+version+'/'
        for it in range(0,len(tilet)):
            date = str(datet[it])
            tile = str(tilet[it])
            print(tile)

            wz = ft['TILEID'] == int(tile)
            
            fz = ft[wz]
#             for rzdir in rzdirs:
#                 rzf = svdir+'redshift_comps/'+rzdir+'/'+version+'/'+tp+'/'+tile+'_'+tp+'zinfo.fits'
#                 if os.path.isfile(rzf):
#                     print('found '+rzf)
#                     fz = Table.read(rzf)
#                     fz['subset'] = np.core.defchararray.add(f['subset'],rzdir)
#                     print(np.unique(fz['subset']))
#                     tz = vstack([tz,fz])
#                     print(np.unique(tz['subset']))
#                 
#             print(len(tz))
            tt=Table.read(dirvi+tp[:3]+'/'+'desi-vi_'+tp[:3]+'_tile'+tile+'_nightdeep_merged_all_'+date+'.csv',format='pandas.csv')
            tt.keep_columns(['TARGETID','best_z','best_quality','best_spectype','all_VI_issues','all_VI_comments','merger_comment','N_VI'])
            print(len(tt),len(fz))
            tj = join(fz,tt,join_type='left',keys='TARGETID')
            print(len(tj))
            tj['N_VI'].fill_value = 0
            tj['N_VI'] = tj['N_VI'].filled() #should easily be able to select rows with N_VI > 0 to get desired info
            #tj['TILEID'] = tile
            tj.write(dirz+'/'+tp+'/'+tile+'_'+tp+'zinfo_wVI.fits',format='fits',overwrite=True)
            print('wrote file with VI info to '+dirz+'/'+tp+'/'+tile+'_'+tp+'zinfo_wVI.fits')
            gt.append(tile)


        print(gt)
        dt = Table.read(dirz+'/'+tp+'/'+str(gt[0])+'_'+tp+'zinfo_wVI.fits')
        for it in range(1,len(gt)):
            dtn = Table.read(dirz+'/'+tp+'/'+str(gt[it])+'_'+tp+'zinfo_wVI.fits')
            dt = vstack([dt,dtn])

        print(np.unique(dt['TILEID']))
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
           
            if ln['RZR'] != 'N':
                zfitdir = '/global/cfs/cdirs/desi/users/rongpu/redux/cascades/'+ln['RZR']+'/'+str(ln['TILEID'])
            else:
                zfitdir = tiledir+str(ln['TILEID'])+'/'+ln['subset']+'/'    
            
            fl = zfitdir+'/redrock-'+str(ln['PETAL_LOC'])+'-'+str(ln['TILEID'])+'-'+ln['subset']+'.h5'
            
            zfits = zi.get_zfits(fl,ln['TARGETID'])
            for jj in range(1,5):
                for col in cols:
                    dt[col+'_'+str(jj)][ii] = zfits[jj][col]
            if ii%1000 == 0:
                print(ii)

        #dt.sort('TARGETID')
        outfall = dirz +'/'+tp+'/allVItiles_'+tp+'zinfo_wVI.fits'
        dt.write(outfall,format='fits', overwrite=True) 
        print('wrote to '+outfall)

