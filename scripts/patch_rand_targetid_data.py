import LSS.common_tools as common
from astropy.table import Table, join, vstack
from astropy.coordinates import SkyCoord
import astropy.units as u
import argparse
import numpy as np
from numpy import ma
import glob
import fitsio


parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="base directory for catalogs",default='/dvs_ro/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--tracers", help="all runs all for given survey",default='all')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--blinded",help="blinded or unblinded catalogs",default='unblinded')
parser.add_argument("--nran",help="number of randoms to process",default=1,type=int)
parser.add_argument("--par",help="process randoms in parallel?",default='n')
args = parser.parse_args()



indir = args.basedir+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'+args.blinded+'/'
regl = ['NGC','SGC']
tracers = ['LRG','ELG_LOPnotqso','QSO','BGS_BRIGHT-21.5']
if args.tracers != 'all':
    tracers = [args.tracers]
for tracer in tracers:
    # load catalogs
    datal = []
    for reg in regl:
        data_cat_fn = indir +tracer+'_'+reg+'_clustering.dat.fits'
        data = Table.read(data_cat_fn,memmap=True)
        data.rename_column('TARGETID', 'TARGETID_DATA')
        data.keep_columns(['Z','WEIGHT_COMP','TARGETID_DATA','WEIGHT_SYS'])
        datal.append(data)
    data = vstack(datal)
    
    def _fix_tid_ran(rann):
        for reg in regl:
            ran_cat_fn = indir +tracer+'_'+reg+'_'+str(rann)+'_clustering.ran.fits'
            print(f"loading {ran_cat_fn}")
            ran = Table.read(ran_cat_fn,memmap=True)
            # remove current random TARGETID_DATA
            ran.remove_column('TARGETID_DATA')
            # join catalogs
            join_cat = join(ran, data, keys = ['Z','WEIGHT_COMP','WEIGHT_SYS'], join_type ='left')

            len_mask = len(join_cat[join_cat['TARGETID_DATA'] == ma.masked])
            if len_mask == 0 and len(ran) == len(join_cat):
                common.write_LSS(join_cat,ran_cat_fn.replace('dvs_ro','global'))
                
                print(tracer,reg,str(rann))
                print('SUCCESS')
            else:
                print(tracer,reg,str(rann))
                print('FAILURE!')
                print(f"number of null TARGETID_DATA in joined catalog = {len_mask}")
                print(f"length of original random catalog = {len(ran)}")
                print(f"length of new joined catalog  ... = {len(join_cat)}")

    inds = np.arange(0,args.nran)
    if args.par == 'n':
        for rn in inds:
            _fix_tid_ran(rn)
    if args.par == 'y':
        from multiprocessing import Pool
        with Pool() as pool:
            res = pool.map(_fix_tid_ran, inds)
    
    #for reg in regl:
    #    for rann in range(0,args.nran):
    #        ran_cat_fn = indir +tracer+'_'+reg+'_'+str(rann)+'_clustering.ran.fits'
    #        ran = Table.read(ran_cat_fn,memmap=True)
    #        # remove current random TARGETID_DATA
    #        ran.remove_column('TARGETID_DATA')
    #        # join catalogs
    #        join_cat = join(ran, data, keys = ['Z','WEIGHT_COMP','WEIGHT_SYS'], join_type ='left')

    #        len_mask = len(join_cat[join_cat['TARGETID_DATA'] == ma.masked])
    #        if len_mask == 0 and len(ran) == len(join_cat):
    #            common.write_LSS(join_cat,ran_cat_fn.replace('dvs_ro','global'))
    #            
    #            print(tracer,reg,str(rann))
    #            print('SUCCESS')
    #        else:
    #            print(tracer,reg,str(rann))
    #            print('FAILURE!')
    #            print(f"number of null TARGETID_DATA in joined catalog = {len_mask}")
    #            print(f"length of original random catalog = {len(ran)}")
    #            print(f"length of new joined catalog  ... = {len(join_cat)}")




    

    
    