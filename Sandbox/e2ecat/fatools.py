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

def mkfa(targetf,tileff,rd,fad,srun,nrun,DESIMODEL):
	'''
	make the fiberassign files needed
	targetf is the target file (e.g., an mtl file)
	tilef is the root string for the tile files produced for each epoch
	rd is the output directory
	fad is the directory for the data fiberassign files
	srun is the initial epoch
	nrun is the number of epochs
	DESIMODEL is the directory for where to find the focal plane model for running these
	'''
	os.environ['DESIMODEL'] = DESIMODEL
	#targetf = e2eout +program+'/randoms_mtl_cuttod.fits' #above file, cut to ~e2e area with significant padding
	

	#use fiberassign tools to read in randoms to be assigned
	tgs = Targets()
	load_target_file(tgs,targetf)
	print('loaded target file '+targetf)
	tree = TargetTree(tgs, 0.01)

	for run in range(srun,srun+nrun):
		#make the tile file for this run
		#e2e.mke2etiles(run,program=program)
		tilef = tileff+str(run)+'.fits'
		#+str(run)+'.fits'


		randir = rd +str(run)
		#+str(run)
		if os.path.isdir(randir):
			ofls = glob.glob(randir+'/*')
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
