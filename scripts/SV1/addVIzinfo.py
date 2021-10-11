import os
import sys
import numpy as np
import argparse
import fitsio
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt
#sys.path.append('../py')

#from this package
import LSS.zcomp.zinfo as zinfo


#assumes gathSV_zinfo.py has been run for all tracer types for arguments below

parser = argparse.ArgumentParser()
parser.add_argument("--release", help="what spectro release to use, e.g. blanc or daily",default='blanc') #eventually remove this and just gather everything
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version to load from",default='test')
args = parser.parse_args()
print(args)

basedir = args.basedir
release = args.release
version = args.version


types = ['QSO','ELG','LRG','BGS_ANY']#,'MWS']

tiles = {'LRG':[80605,80609],'ELG':[80606,80608],'QSO':[80605,80607,80609],'BGS_ANY':[80613]}
dates = {'LRG':[210224,21030],'ELG':[210218,210208],'QSO':[210223,210214,210210],'BGS_ANY':[210202]}


dirvi = '/global/cfs/cdirs/desi/sv/vi/TruthTables/Blanc/'
svdir = basedir+'/SV1/'
dirz = svdir+'redshift_comps/'+release+'/'+version+'/'

for i in range(0,len(types)):
    tp =types[i]
    tilet = tiles[tp]
    datet = dates[tp]
    for it in range(0,len(tilet)):
        date = str(datet[it])
        tile = str(tilet[it])
        tt=Table.read(dirvi+tp[:3]+'/'+'desi-vi_'+tp[:3]+'_tile'+tile+'_nightdeep_merged_all_'+date+'.csv',format='pandas.csv')
        tt.keep_columns(['TARGETID','best_z','best_quality','best_spectype','all_VI_issues','all_VI_comments','merger_comment','N_VI'])
        tz = Table.read(dirz+'/'+tp+'/'+tile+'_'+tp+'zinfo.fits')
        tj = join(tz,tt,join_type='left',keys='TARGETID')
        tj['N_VI'].fill_value = 0
        tj['N_VI'] = tj['N_VI'].filled() #should easily be able to select rows with N_VI > 0 to get desired info
        tj.write(dirz+'/'+tp+'/'+tile+'_'+tp+'zinfo_wVI.fits',format='fits',overwrite=True)
        print('wrote file with VI info to '+dirz+'/'+tp+'/'+tile+'_'+tp+'zinfo_wVI.fits')
    if len(tilet) > 1:
        dt = Table.read(dirz+'/'+tp+'/'+str(tilet[0])+'_'+tp+'zinfo_wVI.fits')
        for it in range(1,len(tilet)):
            dtn = Table.read(dirz+'/'+tp+'/'+str(tilet[it])+'_'+tp+'zinfo_wVI.fits')
            dt = vstack([dt,dtn])

    
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
            zfits = zinfo.get_zfits(ln['TILEID'],ln['PETAL_LOC'],ln['subset'],ln['TARGETID'],release)
            for jj in range(1,5):
                for col in cols:
                    dt[col+'_'+str(jj)][ii] = zfits[jj][col]
            if ii%1000 == 0:
                print(ii)

        #dt.sort('TARGETID')
        outfall = dirz +'/'+tp+'/allVItiles_'+tp+'zinfo_wVI.fits'
        dt.write(outfall,format='fits', overwrite=True) 
        print('wrote to '+outfall)
                