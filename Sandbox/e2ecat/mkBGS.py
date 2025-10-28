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

print('Need to account for "legacy" LRG redshifts')

#import from catalog code
import e2etools as e2e


E2EDIR = '/global/homes/m/mjwilson/desi/survey-validation/svdc-spring2020f-onepercent/'
#E2EDIR = os.environ['E2EDIR']
print('end to end directory is')
print(E2EDIR)

#directories for inputs and output
e2ein  = E2EDIR                     # Most recent: f. 
e2eout = E2EDIR + 'run/catalogs/'


targroot  = '/project/projectdirs/desi/target/catalogs/dr8/0.31.1/targets/main/resolve/targets-dr8'
ranf      = '/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randomsall/randoms-inside-dr8-0.31.0-all.fits' #DR8 imaging randoms file

#now the routines in e2etools that need these will have them as globals
e2e.setglobals(e2ein,e2eout,targroot,ranf)

elgandlrgbits = [1,5,6,7,8,9,11,12,13] #the combination of mask bits proposed for LRGs and ELGs, for simplicity using them again for quasars
#run through steps to make LRG catalogs

type = 60 #target bit for BGS
program = 'bright'
#epochs
srun = 0

# number of epochs
nrun = 5


#list of independent tasks to perform
mkrandoms = False #make randoms specific for type/observing program
farandoms = False #run randoms through fiberassign; doesn't need to be done if already done for LRGs
combran = False #concatenate random files and match randoms from FAVAIL back to full info using targetID; doesn't need to be done if already done for LRGs
matchran = False
combtar = False #concatenate target files; doesn't need to be done if already done for LRGs 
matchtar = False #match targets to mtl info and to zcat info; doesn't need to be done if already done for LRGs
plotntile = False
mkfullran = True #make "full" catalog for randoms
mkfulldat = True #make "full" catalog for data
mkclusdat = True #make "clustering" catalog for data
mkclusran = True #make clustering catalog for randoms



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

if mkrandoms:
	e2e.mkran_type(type,program)

if farandoms:
	targetf = e2eout +program+'/randoms_mtl_cuttod.fits' #above file, cut to ~e2e area with significant padding
	#random are only for DARK right now, need to work something out for BRIGHT and GRAY

	#use fiberassign tools to read in randoms to be assigned
	tgs = Targets()
	load_target_file(tgs,targetf)
	print('loaded target file '+targetf)
	tree = TargetTree(tgs, 0.01)

	for run in range(srun,srun+nrun):
		#make the tile file for this run
		e2e.mke2etiles(run,program=program)
		tilef = e2eout+program+'/e2etiles_run'+str(run)+'.fits'


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
		
if combran: #concatenate random files, match to mtl
	e2e.combran(srun,nrun,program)

if matchran:
	e2e.matchran(program)
 
if combtar: #concatenate target files, match to mtl and zcat info
	e2e.combtargets(srun,nrun,program)

if matchtar:	
	rmax=nrun-1
	e2e.matchtar(program,rmax)
	e2e.matchzcattar(program,rmax)
	
if plotntile:
	e2e.plotrntile(program)
	e2e.plottntile(program)
	e2e.plotzprobvsntile(program,60)
	
if mkfullran:
    e2e.mkfullran('BGS','bright',elgandlrgbits)

if mkfulldat:
    e2e.mkfulldat('BGS','bright',elgandlrgbits)

#needs to happen before randoms so randoms can get z and weights
if mkclusdat:
    e2e.mkclusdat('BGS','bright')

if mkclusran:
    e2e.mkclusran('BGS','bright')

   	
		
