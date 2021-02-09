import os
import argparse
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt


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


types = ['ELG','LRG','BGS_ANY']#'QSO',
tiles = ['80608','80609','80613']#'na',]
dates = ['210208','21030','210202']

dirvi = '/global/cfs/cdirs/desi/sv/vi/TruthTables/Blanc/'
svdir = basedir+'/SV1/'
dirz = svdir+'redshift_comps/'+release+'/'+version+'/'

for i in range(0,len(types)):
    tp =types[i]
    tile = tiles[i]
    date = dates[i]
    tt=Table.read(dirvi+tp[:3]+'/'+'desi-vi_'+tp[:3]+'_tile'+tile+'_nightdeep_merged_all_'+date+'.csv',format='pandas.csv')
    tt.keep_columns(['TARGETID','best_z','best_quality','best_spectype','all_VI_issues','all_VI_comments','merger_comment','N_VI'])
    tz = Table.read(dirz+'/'+tp+'/'+tile+'_'+tp+'zinfo.fits')
    tj = join(tz,tt,join_type='left',keys='TARGETID')
    tj['N_VI'].fill_value = 0
    tj['N_VI'] = tj['N_VI'].filled() #should easily be able to select rows with N_VI > 0 to get desired info
    tj.write(dirz+'/'+tp+'/'+tile+'_'+tp+'zinfo_wVI.fits',format='fits',overwrite=True)
    print('wrote file with VI info to '+dirz+'/'+tp+'/'+tile+'_'+tp+'zinfo_wVI.fits')
