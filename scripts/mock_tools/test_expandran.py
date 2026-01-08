import numpy as np
from astropy.table import Table, join, vstack
import os,sys
from astropy.io import fits
import fitsio
import time

import LSS.common_tools as common
import h5py
import hdf5plugin

def expand_ran(rann,tracer='LRG',reg='NGC',in_dir='/dvs_ro/cfs/cdirs/desi/mocks/cai/LSS/DA2/mocks/holi_v1/altmtl201/loa-v1/mock201/LSScats/',orig_ran_dir='/dvs_ro/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2/',prog='dark',rancols=['TARGETID','RA','DEC'],datacols=['TARGETID','Z']):
    t0 = time.time()
    in_ran = common.read_hdf5_blosc(orig_ran_dir+prog+'_'+str(rann)+'_full_noveto.ran.h5',columns=rancols)
    in_table = common.read_hdf5_blosc(in_dir+tracer+'_'+reg+'_'+str(rann)+'_clustering.ran.h5',columns=['TARGETID','TARGETID_DATA','WEIGHT','NX'])
    tran = time.time()
    print(str(rann)+' read original randoms;'+str(tran-t0))
    
    olen = len(in_table)
    tids, in_ind, orig_ind = np.intersect1d(in_table['TARGETID'], in_ran['TARGETID'], return_indices=True)
    in_table = in_table[in_ind]
    print(np.array_equal(tids,in_table['TARGETID']))
    in_ran = in_ran[orig_ind]
    for col in rancols:
        if col != 'TARGETID':
            in_table[col] = in_ran[col]
    #in_table = join(in_table,in_ran,keys=['TARGETID'])
    t1 = time.time()
    print(str(rann)+' joined to original randoms;'+str(t1-t0))
    del in_ran
    regl = ['NGC','SGC']
    datal = []
    for reg in regl:
        datal.append(common.read_hdf5_blosc(in_dir+tracer+'_'+reg+'_clustering.dat.h5',columns=datacols))
    in_data = vstack(datal)
    t2 = time.time()
    print(str(rann)+' stacked data;'+str(t2-t0))
    del datal    
    in_data.rename_column('TARGETID', 'TARGETID_DATA')
    #in_table = join(in_table,in_data,keys=['TARGETID_DATA'])
    sorted_idx = np.argsort(in_data['TARGETID_DATA'])
    idx_in_sorted = np.searchsorted(in_data['TARGETID_DATA'], in_table['TARGETID_DATA'], sorter=sorted_idx)
    indices = sorted_idx[idx_in_sorted]
    for col in datacols:
        if col != 'TARGETID':
            in_table[col] = in_data[col][indices]

    t3 = time.time()
    print(str(rann)+' done;'+str(t3-t0))
    print(olen,len(in_table))
    return in_table

if __name__ == '__main__':
    #test = expand_ran(0)
    #if par:
    from multiprocessing import Pool
    #    import sys
        #N = int(sys.argv[2])
        #N = 32
    #    N = rx-rm#+1
        #p = Pool(N)
    #    inds = []
    #    for i in range(rm,rx):
    #        inds.append(i)
        #with sharedmem.MapReduce() as pool:
        #pool = sharedmem.MapReduce(np=6)
        #with pool:
    inds = np.arange(0,18)
    with Pool() as pool:
            #def reduce(ii, r):
            #    logger.info('chunk done '+str(ii))
            #    return r
           pool.map(expand_ran,inds)#,reduce=reduce)

        #p.map(doran,inds)
    #else:
    #    for i in range(rm,rx):
    #       doran(i)
