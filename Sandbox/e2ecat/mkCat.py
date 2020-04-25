'''
one executable to make any of the four kinds of tracer, just setting type
target_type can be LRG, QSO, BGS, or ELG
'''

#target_type = 'BGS'
target_type = 'ELG'

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
import e2etools as e2e
import fatools as fa

DESIMODEL = '/global/homes/m/mjwilson/desi/survey-validation/svdc-spring2020f-onepercent/desimodel/'
os.environ['DESIMODEL'] = DESIMODEL
E2EDIR = '/global/homes/m/mjwilson/desi/survey-validation/svdc-spring2020f-onepercent/' # Most recent: f.
#E2EDIR = os.environ['E2EDIR']
print('end to end directory is')
print(E2EDIR)

#directories for inputs and output
e2ein  = E2EDIR                      
e2eout = E2EDIR + 'run/catalogs/'

if os.path.isdir(E2EDIR + 'run/catalogs/logfiles'):
	pass
else:
	os.mkdir(E2EDIR + 'run/catalogs/logfiles')	

	

logf = open(e2eout+'/logfiles/mkCat'+str(datetime.today())+'.log','w')
logf.write('mkCat.py run for type '+str(type)+'\n')
logf.write('inputs root '+e2ein+'\n')
logf.write('output root '+e2ein+'\n')

targroot  = '/project/projectdirs/desi/target/catalogs/dr8/0.31.1/targets/main/resolve/targets-dr8'
ranf      = '/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randomsall/randoms-inside-dr8-0.31.0-all.fits' #DR8 imaging randoms file

#now the routines in e2etools that need these will have them as globals
e2e.setglobals(e2ein,e2eout,targroot,ranf)

elgandlrgbits = [1,5,6,7,8,9,11,12,13] #the combination of mask bits proposed for LRGs and ELGs, for simplicity using them again for quasars
#run through steps to make LRG catalogs

if target_type == 'BGS':
	print('Need to account for "legacy" LRG redshifts')
	type = 60 #target bit for BGS
	program = 'bright'
	#epochs
	srun = 0
	# number of epochs
	nrun = 5
	imbits =  elgandlrgbits #mask bits for imaging

if target_type == 'LRG':
	type = 0 #target bit for LRG
	program = 'dark'
	#epochs
	srun = 0
	# number of epochs
	nrun = 7
	imbits =  elgandlrgbits #mask bits for imaging

if target_type == 'QSO':
	type = 2 #target bit for QSO
	program = 'dark'
	#epochs
	srun = 0
	# number of epochs
	nrun = 7
	imbits =  elgandlrgbits #mask bits for imaging

if target_type == 'ELG':
	type = 1 #target bit for ELG
	program = 'gray'
	#epochs
	srun = 0
	# number of epochs
	nrun = 7
	imbits =  elgandlrgbits #mask bits for imaging


logf.write('used following settings:\n')
logf.write("\
type = "+str(type)+"\n\
program = "+program+"\n\
srun = "+str(srun)+"\n\
nrun = "+str(nrun)+"\n\
imbits =  "+str(imbits)+"\n\
\n\
")


#list of independent tasks to perform
mkrandoms = False #make randoms specific for type/observing program
farandoms = False #run randoms through fiberassign; doesn't need to be done if already done for LRGs
combran = False #concatenate random files and match randoms from FAVAIL back to full info using targetID; doesn't need to be done if already done for LRGs
matchran = False
combtar = False #concatenate target files; doesn't need to be done if already done for LRGs 
matchtar = False #match targets to mtl info and to zcat info; doesn't need to be done if already done for LRGs
plotntile = False
plotzeff = False
plottilehist = False
mkfullran = True #make "full" catalog for randoms
mkfulldat = True #make "full" catalog for data
mkclusdat = True #make "clustering" catalog for data
mkclusran = True #make clustering catalog for randoms
plotfoot = True
plottilecomp = True




if mkrandoms:
	e2e.mkran_type(type,program)
	logf.write('ran mkrandoms\n')

if farandoms:
	for run in range(srun,srun+nrun):
		#make the tile file for this run
		e2e.mke2etiles(run,program=program)

	targetf = e2eout +program+'/randoms_mtl_cuttod.fits' #above file, cut to ~e2e area with significant padding
	tiles = e2eout+program+'/e2etiles_run'
	rd = e2eout+program+'/randoms/'
	if program == 'gray':
		fad = e2ein+'run/quicksurvey/dark/'
	else:
		fad = e2ein+'run/quicksurvey/'+program+'/'
	fa.mkfa(targetf,tiles,rd,fad,srun,nrun,DESIMODEL)
	logf.write('ran farandoms\n')
		
if combran: #concatenate random files, match to mtl
	e2e.combran(srun,nrun,program)
	logf.write('ran combran\n')

if matchran:
	e2e.matchran(program)
	logf.write('ran matchran\n')
 
if combtar: #concatenate target files, match to mtl and zcat info
	e2e.combtargets(srun,nrun,program)
	logf.write('ran combtargets\n')

if matchtar:	
	rmax=nrun-1
	e2e.matchtar(program,rmax)
	e2e.matchzcattar(program,rmax)
	logf.write('ran matchtar\n')
	
if plotntile:
	e2e.plotrntile(program)
	e2e.plottntile(program)
	logf.write('ran plotntile\n')

if plotzeff:
	e2e.plotzprobvsntile(program,type)
	logf.write('ran plotzeff\n')
	
if plottilehist:
	e2e.comphistNT(program)	
	
if mkfullran:
    e2e.mkfullran(target_type,program,imbits)
    logf.write('ran mkfullran\n')

if mkfulldat:
    e2e.mkfulldat(target_type,program,imbits)
    logf.write('ran mkfulldat\n')

#needs to happen before randoms so randoms can get z and weights
if mkclusdat:
    e2e.mkclusdat(target_type,program)
    logf.write('ran mkclusdat\n')

if mkclusran:
    e2e.mkclusran(target_type,program)
    logf.write('ran mkclusran\n')

if plotfoot:
	e2e.plotcompdr_full(target_type,program)  
	
if plottilecomp:
	e2e.plotcompvsntile(target_type,program)	 	
		
