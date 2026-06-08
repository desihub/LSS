'''
will write out fiberassignment files for each tile with the FASSIGN, FTARGETS, FAVAIL HDUS
these are what are required to determine the geometry of what fiberassign thinks could have been observed and also match to actual observations (though FASSIGN is not really necessary)
targetf is file with all targets to be run through
tile files and output directories are generated based on the information in the run directories
'''                                

e2ein = '/global/homes/m/mjwilson/desi/survey-validation/svdc-spring2020e-onepercent/' #Mike's most recent
e2eout = '/project/projectdirs/desi/users/ajross/catalogs/e2eoneper/' #I will keep over-writing most recent here until there are outputs to really be tested


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

#import from catalog code
from e2etools import mke2etiles

#import from fiberassign
from fiberassign.targets import (TargetsAvailable)
from fiberassign.utils import option_list, GlobalTimers
from fiberassign.hardware import load_hardware
from fiberassign.tiles import load_tiles, Tiles
from fiberassign.targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                                 TARGET_TYPE_SUPPSKY,
                                 TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                                 Targets, TargetsAvailable, TargetTree,
                                 LocationsAvailable, load_target_file)
from fiberassign.assign import (Assignment, write_assignment_fits,
                                write_assignment_ascii, merge_results,
                                read_assignment_fits_tile)                                 

#desimodel
import desimodel
import desimodel.io as dmio




srun = int(sys.argv[1]) #starting run
nrun = int(sys.argv[2]) #number of runs
program = str(sys.argv[3])

#test
run = srun
fafls = glob.glob(e2ein+'run/quicksurvey/'+program+'/'+str(run)+'/fiberassign/fiberassign*')
hd = fitsio.read_header(fafls[0])
dt = hd['FA_RUN']


#targetf = '/project/projectdirs/desi/users/ajross/catalogs/minisv2/random/random_mtl.fits' #2e8 randoms from dr8 randoms, with columns needs for fa, 
targetf = '/project/projectdirs/desi/users/ajross/catalogs/e2eoneper/'+program+'/randoms/randoms_mtl_cuttod.fits' #above file, cut to ~e2e area with significant padding
#random are only for DARK right now, need to work something out for BRIGHT and GRAY

#use fiberassign tools to read in randoms to be assigned
tgs = Targets()
load_target_file(tgs,targetf)
print('loaded target file '+targetf)
tree = TargetTree(tgs, 0.01)




for run in range(srun,srun+nrun):
	#make the tile file for this run
	mke2etiles(run,program=program)
	tilef = e2eout+'e2etiles_run'+str(run)+'.fits'


	randir = e2eout+program+'/randoms/'+str(run)
	if os.path.isdir(randir):
		pass
	else:	
		os.mkdir(randir)

	fafls = glob.glob(e2ein+'run/quicksurvey/'+program+'/'+str(run)+'/fiberassign/fiberassign*')
	hd = fitsio.read_header(fafls[0])
	dt = hd['FA_RUN']

	hw = load_hardware(rundate=dt)

	tiles = load_tiles(tiles_file=tilef)
	tgsavail = TargetsAvailable(hw, tgs, tiles, tree)
	favail = LocationsAvailable(tgsavail)
	#del tree
	asgn = Assignment(tgs, tgsavail, favail)
	asgn.assign_unused(TARGET_TYPE_SCIENCE)
	write_assignment_fits(tiles, asgn, out_dir=randir, all_targets=True)
	print('wrote assignment files to '+randir)	