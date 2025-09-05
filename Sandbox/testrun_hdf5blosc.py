#!/usr/bin/env python

"""
Test run converting randoms
"""

import os, sys, glob
import numpy as np
from astropy.table import Table

import fitsio
import numpy as np
import h5py
import hdf5plugin

import LSS.common_tools as common


# create logger
import logging
logname = 'LSSran'
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

common.printlog('starting script',logger)
tracer = 'LRG' 
ddir = '/pscratch/sd/a/ajross/DA2/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl5/loa-v1/mock5/LSScats/'
inds = np.arange(18) #18 random files
regions = ['NGC','SGC']
for reg in regions:
	flin = ddir + tracer + '_'+reg
	def _parfun(rannum):
		indat = fitsio.read(flin+'_'+str(rannum)+'_clustering.ran.fits')
		out_fn = flin+'_'+str(rannum)+'_clustering.ran.h5'
		common.write_LSShdf5_scratchcp(indat, out_fn,logger=logger)


	from multiprocessing import Pool
	with Pool() as pool:
		res = pool.map(_parfun, inds)
common.printlog('finished script',logger)