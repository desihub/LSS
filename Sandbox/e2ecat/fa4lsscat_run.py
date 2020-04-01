import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import fitsio
import desimodel
import glob

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
import desimodel.io as dmio

import sys
import os
test
run = int(sys.argv[1])

e2ein = '/global/homes/m/mjwilson/desi/survey-validation/svdc-spring2020b-onepercent/'
e2eout = '/project/projectdirs/desi/users/ajross/catalogs/e2eoneper/'
randir = e2eout+'randoms/'+str(run)
if os.path.isdir(randir):
	pass
else:	
	os.mkdir(randir)

fafls = glob.glob(e2ein+'run/quicksurvey/'+str(run)+'/fiberassign/*')
hd = fitsio.read_header(fafls[0])
dt = hd['FA_RUN']
#runtime = datetime.strptime(dt,'%Y-%m-%dT%H:%M:%S')

targetf = '/project/projectdirs/desi/users/ajross/catalogs/minisv2/random/random_mtl.fits' #2e8 randoms from dr8 randoms, with columns needs for fa

#make the tile file for this run
from cattools import mke2etiles
mke2etiles(run)
tilef = e2eout+'e2etiles_run'+str(run)+'.fits'

# First get the starting focalplane from desimodel
fp, exclude, state, tmstr = dmio.load_focalplane(runtime)


'''
will write out fiberassignment files for each tile with the FASSIGN, FTARGETS, FAVAIL HDUS
these are what are required to determine the geometry of what fiberassign thinks could have been observed and also match to actual observations (though FASSIGN is not really necessary)
targetf is file with all targets to be run through
tilef lists the tiles to "assign"
dirout is the directory where this all gets written out !make sure this is unique for every different target!
'''                                
tgs = Targets()
load_target_file(tgs,targetf)
print('loaded target file '+targetf)
tree = TargetTree(tgs, 0.01)
hw = load_hardware(rundate=dt)
tiles = load_tiles(tiles_file=tilef)
tgsavail = TargetsAvailable(hw, tgs, tiles, tree)
favail = LocationsAvailable(tgsavail)
del tree
asgn = Assignment(tgs, tgsavail, favail)
asgn.assign_unused(TARGET_TYPE_SCIENCE)
write_assignment_fits(tiles, asgn, out_dir=randir, all_targets=True)
print('wrote assignment files to '+randir)	