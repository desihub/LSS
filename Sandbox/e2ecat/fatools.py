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

def mkfa(targetf,tiles,rd,fad,srun,nrun):
	'''
	make the fiberassign files needed
	'''
	#targetf = e2eout +program+'/randoms_mtl_cuttod.fits' #above file, cut to ~e2e area with significant padding
	#random are only for DARK right now, need to work something out for BRIGHT and GRAY

	#use fiberassign tools to read in randoms to be assigned
	tgs = Targets()
	load_target_file(tgs,targetf)
	print('loaded target file '+targetf)
	tree = TargetTree(tgs, 0.01)

	for run in range(srun,srun+nrun):
		#make the tile file for this run
		#e2e.mke2etiles(run,program=program)
		tilef = tiles+str(run)+'.fits'
		#+str(run)+'.fits'


		randir = rd +str(run)
		#+str(run)
		if os.path.isdir(randir):
			ofls = glob.glob(randir)
			for fl in ofls:
				os.remove(fl) #remove the old files if we are rerunning
		else:	
			os.mkdir(randir)
		
		fafls = glob.glob(fad+str(run)+'/fiberassign/fiberassign*')
		#+str(run)+'/fiberassign/fiberassign*')
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
