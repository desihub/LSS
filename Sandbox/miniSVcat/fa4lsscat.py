import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import fitsio
import desimodel

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

minisvdir = '/project/projectdirs/desi/users/ajross/catalogs/minisv2/'
randir = minisvdir+'random/'

dt = '2020-03-10T00:00:00'
#runtime = datetime.utcnow()

# First get the starting focalplane from desimodel
#fp, exclude, state, tmstr = dmio.load_focalplane(runtime)

def getfatiles(targetf,tilef,dirout=randir):
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
	#hw = load_hardware(focalplane=(fp, exclude, state))
	hw = load_hardware(rundate=dt)
	tiles = load_tiles(tiles_file=tilef)
	tgsavail = TargetsAvailable(hw, tgs, tiles, tree)
	favail = LocationsAvailable(tgsavail)
	del tree
	asgn = Assignment(tgs, tgsavail, favail)
	asgn.assign_unused(TARGET_TYPE_SCIENCE)
	write_assignment_fits(tiles, asgn, out_dir=dirout, all_targets=True)
	print('wrote assignment files to '+dirout)	