'''
one executable to make any of the four kinds of tracer, just setting type
target_type can be LRG, QSO, BGS, or ELG
'''



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

sys.path.append("../")

try:
	target_type = str(sys.argv[1])
except:	
	target_type = 'LRG'

truez=False

if truez:
	print('using true redshifts, whether assigned or not')


#import from catalog code
import e2etools as e2e
import fatools as fa

ver='g'
DESIMODEL = '/global/homes/m/mjwilson/desi/survey-validation/svdc-spring2020'+ver+'-onepercent/desimodel/'
os.environ['DESIMODEL'] = DESIMODEL
E2EDIR = '/global/homes/m/mjwilson/desi/survey-validation/svdc-spring2020'+ver+'-onepercent/' # Most recent: g.
#E2EDIR = os.environ['E2EDIR']
print('end to end directory is')
print(E2EDIR)

#directories for inputs and output
e2ein  = E2EDIR

test = True
if test:
	e2eout = '/project/projectdirs/desi/users/ajross/catalogs/test/'
else:                      
	e2eout = E2EDIR + 'run/catalogs/'
	
print('output directory is '+e2eout)	

#make these directories if they do not exist already
if os.path.isdir(e2eout + 'logfiles'):
	pass
else:
	os.mkdir(e2eout + 'logfiles')	

if os.path.isdir(e2eout + 'bright'):
	pass
else:
	os.mkdir(e2eout + 'bright')	

if os.path.isdir(e2eout + 'dark'):
	pass
else:
	os.mkdir(e2eout + 'dark')	

if os.path.isdir(e2eout + 'gray'):
	pass
else:
	os.mkdir(e2eout + 'gray')	

if os.path.isdir(e2eout + 'bright/randoms'):
	pass
else:
	os.mkdir(e2eout + 'bright/randoms')	

if os.path.isdir(e2eout + 'dark/randoms'):
	pass
else:
	os.mkdir(e2eout + 'dark/randoms')	

if os.path.isdir(e2eout + 'gray/randoms'):
	pass
else:
	os.mkdir(e2eout + 'gray/randoms')	
	

logf = open(e2eout+'/logfiles/mkCat'+str(datetime.today())+'.log','w')
logf.write('using e2e onepercent version '+ver)
logf.write('mkCat.py run for type '+str(type)+'\n')
logf.write('inputs root '+e2ein+'\n')
logf.write('output root '+e2ein+'\n')

targroot  = '/project/projectdirs/desi/target/catalogs/dr8/0.31.1/targets/main/resolve/targets-dr8'
ranf      = '/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randomsall/randoms-inside-dr8-0.31.0-all.fits' #DR8 imaging randoms file
ranfmtl = '/project/projectdirs/desi/users/ajross/catalogs/minisv2/random/random_mtl.fits' #random file with 2e8 rows and columns needed for fiberassign

#now the routines in e2etools that need these will have them as globals
e2e.setglobals(e2ein,e2eout,targroot,ranf,ranfmtl)

elgandlrgbits = [1,5,6,7,8,9,11,12,13] #the combination of mask bits proposed for LRGs and ELGs, for simplicity using them again for quasars
#run through steps to make LRG catalogs

omega_matter = 0.31 # for Nbar and then FKP weights

if ver == 'g':
	nrundark = 15
	nrunbright = 9

if ver == 'f':
	nrundark = 7
	nrunbright = 5


if target_type in ['BGS','BGS_BRIGHT','BGS_FAINT','BGS_BRIGHT_HIP','BGS_FAINT_HIP']:
	print('Need to account for "legacy" LRG redshifts')
	type = 60 #target bit for BGS
	program = 'bright'
	#epochs
	srun = 0
	# number of epochs
	nrun = nrunbright
	imbits =  elgandlrgbits #mask bits for imaging
	P0 = 2500.

if target_type == 'LRG':
	type = 0 #target bit for LRG
	program = 'dark'
	#epochs
	srun = 0
	# number of epochs
	nrun = nrundark
	imbits =  elgandlrgbits #mask bits for imaging
	P0 = 10000

if target_type == 'QSO':
	type = 2 #target bit for QSO
	program = 'dark'
	#epochs
	srun = 0
	# number of epochs
	nrun = nrundark
	imbits =  elgandlrgbits #mask bits for imaging
	P0 = 7000

if target_type == 'ELG':
	type = 1 #target bit for ELG
	program = 'gray'
	#epochs
	srun = 0
	# number of epochs
	nrun = nrundark
	imbits =  elgandlrgbits #mask bits for imaging
	P0 = 4000


logf.write('used following settings:\n')
logf.write("\
type = "+str(type)+"\n\
program = "+program+"\n\
srun = "+str(srun)+"\n\
nrun = "+str(nrun)+"\n\
imbits =  "+str(imbits)+"\n\
truez =  "+str(truez)+"\n\
P0 =  "+str(P0)+"\n\
omega_matter =  "+str(omega_matter)+"\n\
\n\
")


#list of independent tasks to perform
countfavail = False
cutran = False #cut big random file to only occupy one percent footprint
mkrandoms = False #make randoms specific for type/observing program
farandoms = False #run randoms through fiberassign; doesn't need to be done for QSO if already done for LRGs
combran = False #concatenate random files and match randoms from FAVAIL back to full info using targetID; doesn't need to be done if already done for LRGs
matchran = False
combtar = False #concatenate target files; doesn't need to be done if already done for LRGs 
matchtar = False #match targets to mtl info and to zcat info; doesn't need to be done if already done for LRGs
plotntile = False
plotzeff = False
plottilehist = False
mkfulldat = False #make "full" catalog for data
mkprob = True #add fraction with good z at tileloc to full data
mkfullran = True #make "full" catalog for randoms
mkclusdat = True #make "clustering" catalog for data
mkclusran = True #make clustering catalog for randoms
mkNbar = True
fillNZ = True
plotfoot = False
plottilecomp = False
plotfatiledr = False


if cutran:
	e2e.cutran(ver)
	logf.write('ran cutran\n')

if mkrandoms:
	e2e.mkran_type(type,program)
	logf.write('ran mkrandoms\n')

if farandoms:
	for run in range(srun,srun+nrun):
		#make the tile file for this run
		e2e.mke2etiles(run,program=program,ver=ver)

	targetf = e2eout +program+'/randoms_mtl_cuttod.fits' #above file, cut to ~e2e area with significant padding
	tiles = e2eout+program+'/e2etiles_run'
	rd = e2eout+program+'/randoms/'
	if program == 'gray':
		fad = e2ein+'run/quicksurvey/dark/'
	else:
		fad = e2ein+'run/quicksurvey/'+program+'/'
	fa.mkfa(targetf,tiles,rd,fad,srun,nrun,DESIMODEL)
	logf.write('ran farandoms\n')

if countfavail:
	e2e.count_tarfavail(srun,nrun,program)

		
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
	print(rmax)
	if truez == False:
		e2e.matchtar(program,rmax)
		#rmax = nrun-1
		e2e.matchzcattar(program,rmax)
	else:
		e2e.matchzcattar_nofa(program)
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
	
if mkfulldat:
	print(truez)
	e2e.mkfulldat(target_type,program,imbits,truez=truez)
	logf.write('ran mkfulldat\n')
	print('ran mkfulldat\n')

if mkprob:
	e2e.get_tilelocweight(target_type,program)
	logf.write('ran get_tilelocweight\n')
	print('ran get_tilelocweight\n')

if mkfullran:
    e2e.mkfullran(target_type,program,imbits,truez=truez)
    logf.write('ran mkfullran\n')
    print('ran mkfullran\n')

#needs to happen before randoms so randoms can get z and weights
if mkclusdat:
    print(truez)
    e2e.mkclusdat(target_type,program,truez=truez)
    logf.write('ran mkclusdat\n')
    print('ran mkclusdat\n')

if mkclusran:
    e2e.mkclusran(target_type,program,truez=truez)
    logf.write('ran mkclusran\n')
    print('ran mkclusran\n')
    
if mkNbar:
	e2e.mkNbar(target_type,program,P0=P0,omega_matter=omega_matter,truez=truez)
	logf.write('made nbar\n')
	print('made nbar\n')

if fillNZ:
	e2e.fillNZ(target_type,program,P0=P0,truez=truez)	
	logf.write('put NZ and weight_fkp into clustering catalogs\n')    
	print('put NZ and weight_fkp into clustering catalogs\n')

if plotfoot:
	e2e.plotcompdr_full(target_type,program)  
	
if plottilecomp:
	e2e.plotcompvsntile(target_type,program)	

if plotfatiledr:	
	e2e.compfavail_dr(53478,epoch=14,program=program) 	
	e2e.compfavail_dr(53457,epoch=14,program=program) 	
		
