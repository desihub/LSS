#adapts Mike Wilson's notebook
import glob
import numpy as np
import astropy.io.fits as fits
import argparse
import fitsio
import os

from   astropy.table import Table, join, unique, vstack
from   desiutil.log import get_logger
import ephem
import astropy.units as u

from   desisurvey.config import Configuration
from   astropy.time import Time
from   astropy.table import Table

import LSS.common_tools as common

config = Configuration()

mayall = ephem.Observer()
mayall.lat = config.location.latitude().to(u.rad).value
mayall.lon = config.location.longitude().to(u.rad).value
mayall.elevation = config.location.elevation().to(u.m).value




parser = argparse.ArgumentParser()
parser.add_argument("--outdir", help="base directory for output, default is SCRATCH",default=os.environ['SCRATCH'])
parser.add_argument("--survey", help="main or sv3",default='main')
parser.add_argument("--release", help="e.g., Y1 or DA2",default='main')
parser.add_argument("--prog", help="dark or bright",default='dark')
parser.add_argument("--verspec",help="version for redshift pipeline",default='daily')
parser.add_argument("--test",help="if yes, test a small fraction of the exposures",default='n')
parser.add_argument("--zcatver",help="version for redshift catalogs (starts with verspec jura)",default='')
parser.add_argument("--par",help="whether to process in parallel",default='n')

args = parser.parse_args()

import logging
# create logger
logname = 'get_speccon'
logger = logging.getLogger(logname)
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

logger.info('script is starting')


sw = args.survey
if args.survey == 'sv3':
    sw = 'SV3'

outf = args.outdir+'/specobscon_'+args.prog+'.fits'

datadir   = '/global/cfs/cdirs/desi/spectro/redux/'+args.verspec+'/'
exposures = fitsio.read(datadir + '/exposures-'+args.verspec+'.fits')
if args.test == 'y':
    exposures = exposures[:10]
exposures = Table(exposures)    
nexp = len(exposures)
#if args.test == 'y':
#    nexp = 10
exposures['MOON_ILLUM'] = np.zeros(nexp)

moon = ephem.Moon(mayall)
for ii in range(0,nexp):

    t = Time(exposures[ii]['MJD'], format='mjd')
    moon.compute(t.datetime)

    moon_illum = moon.moon_phase
    exposures[ii]['MOON_ILLUM'] = moon_illum
    
common.printlog('added moon illumination, median is:'+str(np.median(exposures['MOON_ILLUM'])),logger)



addcols = ['HUMIDITY','ZD','ETCTRANS', 'ETCTHRUB', 'ETCSKY', 'ACQFWHM','SLEWANGL','MOONSEP','PMIRTEMP', 'TAIRTEMP','PARALLAC','ROTOFFST','TURBRMS','WINDSPD','WINDDIR']

for col in addcols:
    exposures[col] = np.ones(nexp)*-99


for ii in range(0,nexp):
    es = str(exposures[ii]['EXPID']).zfill(8)
    efn = '/dvs_ro/cfs/cdirs/desi/spectro/data/'+str(exposures[ii]['NIGHT'])+'/'+es+'/desi-'+es+'.fits.fz'
    hh = fitsio.read_header(efn,ext=1)
    if ii//100 == ii/100:
        common.printlog('at exposure '+str(ii)+ ' out of '+str(nexp),logger)
    for col in addcols:
        try:
            exposures[ii][col] = hh[col]
        except:
            pass

for col in addcols:
    selnull = exposures[col] == -99
    common.printlog('fraction null:',logger)
    common.printlog(col+' '+str(len(exposures[selnull])/len(exposures)),logger)               

ocol = ['MOON_ILLUM','EXPID', 'SEEING_ETC', 'AIRMASS', 'EBV', 'TRANSPARENCY_GFA', 'SEEING_GFA', 'SKY_MAG_AB_GFA', 'SKY_MAG_G_SPEC', 'SKY_MAG_R_SPEC', 'SKY_MAG_Z_SPEC', 'EFFTIME_SPEC','EXPTIME','EFFTIME_ETC','EFFTIME_DARK_GFA', 'EFFTIME_BRIGHT_GFA']
tcol = addcols + ocol
exposures = exposures[tcol]

if args.verspec == 'daily':
    dcat = fitsio.read(args.basedir +'/'+sw+'/LSS/'+args.verspec+'/datcomb_'+args.prog+'_spec_zdone.fits')
else:
    dcat = fitsio.read(datadir+'/zcatalog/'+args.zcatver+'/ztile-'+args.survey+'-'+args.prog+'-'+'cumulative.fits')
tids = np.unique(dcat['TILEID'])

mt = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
wd = mt['SURVEY'] == args.survey
wd &= mt['FAPRGRM'] == args.prog
wd &= np.isin(mt['TILEID'],tids)
mtd = mt[wd]

tiles4comb = Table()
tiles4comb['TILEID'] = mtd['TILEID'].astype(int)
tiles4comb['ZDATE'] = mtd['LASTNIGHT']

common.printlog('numbers of tiles, should match:',logger)
common.printlog(str(len(tids))+','+str(len(tiles4comb)),logger)

coadd_fpaths = []

for ii in range(0,len(tiles4comb)):
    '''
    Retrieve coadd paths for all tiles  
    '''
    
    fpath = '{}/tiles/cumulative/{:d}/{:d}'.format(datadir, tiles4comb['TILEID'][ii], tiles4comb['ZDATE'][ii])
    
    
    # Here we grab the path for each coadd under cumulative/tileid/zdate
    fpaths = sorted(glob.glob(fpath + '/' + 'coadd-?-{:d}-thru{}.fits'.format(tiles4comb['TILEID'][ii], tiles4comb['ZDATE'][ii])))
    
    
    coadd_fpaths += [x.replace('global','dvs_ro') for x in fpaths]

print(coadd_fpaths[:12])

def process_coadd(coadd_fpath):
    '''
    Retrieve the input expids for each location on a given (tileid, thru_night).
    
    Note:  
        assuming input expids may differ by location due to quality cuts.  We 
        run round this later by processing simulatenously all locs with the same
        input expids.
    '''
    tileid     = coadd_fpath.split('/')[-3]
    thru_night = coadd_fpath.split('/')[-2]
 
    coadd      = Table.read(coadd_fpath, hdu='EXP_FIBERMAP')
    # coadd
    
    # expids, cnts = np.unique(coadd['EXPID'], return_counts=True) 
    
    # print(len(expids))
    
    condition_cat = coadd['TARGETID', 'LOCATION']
    condition_cat = unique(condition_cat)
    condition_cat.sort('LOCATION')

    condition_cat['TILEID']     = tileid
    condition_cat['THRU_NIGHT'] = thru_night
    condition_cat['IN_EXPIDS']  = 'x' * 50
    
    locs, cnts     = np.unique(condition_cat['LOCATION'].data, return_counts=True)

    assert cnts.max() == 1
    assert np.all(locs == condition_cat['LOCATION'].data)
    
    for i, loc in enumerate(locs):
        coadd_loc  = coadd[(coadd['LOCATION'] == loc) & (coadd['FIBERSTATUS'] == 0)]
      
        loc_expids = '-'.join(np.unique(coadd_loc['EXPID'].data).astype(str).tolist())
        
        condition_cat['IN_EXPIDS'][i] = loc_expids
          
        # print(i, loc_expids)
    common.printlog('done with '+coadd_fpath,logger)        
    return condition_cat

to_process = coadd_fpaths
if args.par == 'y':
    condition_cat = []
    from concurrent.futures import ProcessPoolExecutor            
    with ProcessPoolExecutor() as executor:
        for con_cat in executor.map(process_coadd, to_process):
           condition_cat.append(con_cat)
else:
    condition_cat = [process_coadd(x) for x in to_process]
condition_cat = vstack(condition_cat,metadata_conflicts='silent')

unique_in_expids = np.unique(condition_cat['IN_EXPIDS'].data).tolist()
                
unique_in_expids.remove('')

update_cols   = list(exposures.dtype.names)
update_cols.remove('EXPID')
sum_cols = ['EFFTIME_SPEC','EXPTIME','EFFTIME_ETC','EFFTIME_DARK_GFA', 'EFFTIME_BRIGHT_GFA']
for col in sum_cols:
    update_cols.remove(col)


for col in update_cols:
    condition_cat[col] = -99.
 
for col in sum_cols:
    condition_cat[col] = -99.

    
for in_expids in unique_in_expids:    
    expids        = np.array(in_expids.split('-')).astype(int)

    # Get the exposure conditions for this set of expids.   
    in_exposures  = exposures[np.isin(exposures['EXPID'].data, expids)]
    
    # print(expids)
    # print(in_exposures)
       
    mean_function = lambda x: np.average(x, weights=in_exposures['EFFTIME_SPEC'])

    sum_function = lambda x: np.sum(x)
    
    # Weighted mean of the condition table for this exp. set (weights are efftime_spec)
    if np.max(in_exposures['EFFTIME_SPEC']) > 0: 
        in_exposures               = in_exposures.groups.aggregate(mean_function)
        in_exposures_sum           = in_exposures.groups.aggregate(sum_function)
        #  To be extra sure, we could include matches to TILEID and thru night.     
        to_update                  = condition_cat['IN_EXPIDS'] == in_expids
                
        for col in update_cols:
            condition_cat[col].data[to_update] = in_exposures[col].data[0]
        for col in sum_cols:
            condition_cat[col].data[to_update] = in_exposures_sum[col].data[0]
    else:
         common.printlog('EFFTIME_SPEC all 0: {}'.format(in_expids),logger)   
    common.printlog('Processed: {}'.format(in_expids),logger)


condition_cat.write(outf, format='fits', overwrite=True)
