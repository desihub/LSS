'''
python functions to do various useful date processing/manipulation
'''
import os
import numpy as np
import fitsio
import glob
import astropy.io.fits as fits
from astropy.table import Table,vstack,unique,join#,setdiff
from matplotlib import pyplot as plt
import desimodel.footprint
import desimodel.focalplane #


#targroot  = '/project/projectdirs/desi/target/catalogs/dr8/0.31.1/targets/main/resolve/targets-dr8'
#ranf      = '/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randomsall/randoms-inside-dr8-0.31.0-all.fits' #DR8 imaging randoms file

# E2E
#e2ein     = os.environ['E2EDIR']
#e2ein = '/global/homes/m/mjwilson/desi/survey-validation/svdc-spring2020f-onepercent/'
#e2eout    = e2ein + 'run/catalogs/'



def setglobals(e2einv,e2eoutv,targrootv,ranfv):
	'''
	set paths/files names for important inputs without relying on environment variables
	goal is to always have explicit paths and file names set in the main
	makes it easy to trace back to exact location of everything while ensuring consistency
	has to be a better way, but this is working right now
	'''
	global e2ein
	e2ein = e2einv #where survey files were generated, e.g., '/global/homes/m/mjwilson/desi/survey-validation/svdc-spring2020f-onepercent/'
	global e2eout
	e2eout = e2eoutv #where catalogs are output, e.g., e2ein + 'run/catalogs/'
	global targroot
	targroot = targrootv #where target files are e.g., '/project/projectdirs/desi/target/catalogs/dr8/0.31.1/targets/main/resolve/targets-dr8'
	global ranf
	ranf = ranfv #where randoms are, e.g., '/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randomsall/randoms-inside-dr8-0.31.0-all.fits'
	#global bits
	#bits = bitsv #the imaging mask bits to use

	print('end to end directory is')
	print(e2ein)


def mkran4fa(N=2e8,fout='random_mtl.fits'):
        '''
        cut imaging random file to first N entries and add columns necessary for fiberassignment routines
        this made the initial big one for dark tiles
        '''
        dirout=e2eout+'random/'
        rall = fitsio.read(ranf,rows=np.arange(N))
        rmtl = Table()
        for name in rall.dtype.names:
                rmtl[name] = rall[name]
        rmtl['TARGETID'] = np.arange(len(rall))
        rmtl['DESI_TARGET'] = np.ones(len(rall),dtype=int)*2
        rmtl['SV1_DESI_TARGET'] = np.ones(len(rall),dtype=int)*2
        rmtl['NUMOBS_INIT'] = np.zeros(len(rall),dtype=int)
        rmtl['NUMOBS_MORE'] = np.ones(len(rall),dtype=int)
        rmtl['PRIORITY'] = np.ones(len(rall),dtype=int)*3400
        rmtl['OBSCONDITIONS'] = np.ones(len(rall),dtype=int)
        rmtl['SUBPRIORITY'] = np.random.random(len(rall))
        rmtl.write(dirout+fout,format='fits', overwrite=True)

def mkran_type(dt,program):
	if program == 'bright':
		oc = 2
	if program == 'gray':
		oc = 1
	f = Table.read(e2eout+'dark/randoms_mtl_cuttod.fits') #this has things set for running dark time tiles
	f['DESI_TARGET'] = np.ones(len(f),dtype=int)*2**dt
	f['OBSCONDITIONS'] = np.ones(len(f),dtype=int)*2**oc
	if program == 'bright':
		f['BGS_TARGET'] = np.ones(len(f),dtype=int)*2
	f.write(e2eout+program+'/randoms_mtl_cuttod.fits',format='fits', overwrite=True)

def combran(srun=0,nrun=7,program='dark'):
	dir0 = e2eout+program+'/randoms/'+str(srun)+'/'
	outf=program+'/randoms/randoms_oneper.fits'
	fafls0 = glob.glob(dir0+'fba-*.fits')
	fah = fitsio.read_header(fafls0[0])
	tile = fah['TILEID']
	#
	#exps = fitsio.read(e2ein+'run/survey/complete_exposures_surveysim_fix.fits')
	if program == 'gray':
		exps = fitsio.read(e2ein+'run/quicksurvey/dark/epochs-dark.fits')
		programf = 'dark'
	else:
		exps = fitsio.read(e2ein+'run/quicksurvey/'+program+'/epochs-'+program+'.fits')
		programf = program
	w = exps['TILEID'] == tile
	if len(exps[w]) > 1:
			return 'NEED to deal with multiple exposures of same tile'

	#if exps[w]['EPOCH'][0] == srun:#[0]:
	#       pass
	#else:
	#       return 'first tile was not observed in assigned epoch, fix code'
	expid = exps[w]['EXPID'][0]     
	ep = exps[w]['EPOCH'][0]
	#if expid < 100:
	#       zer = '000000'
	fmap = fitsio.read(e2ein+'run/quicksurvey/'+programf+'/'+str(ep)+'/fiberassign/fibermap-'+str(expid).zfill(8)+'.fits')
	#fmap['FIBERSTATUS'] = 0
	#print('set fiberstatus all to 0; fix this once propagated to zcat')
	wloc = fmap['FIBERSTATUS'] == 0
	gloc = fmap[wloc]['LOCATION']
	fa = Table.read(fafls0[0],hdu='FAVAIL')
	wg = np.isin(fa['LOCATION'],gloc)
	fg = fa[wg]
	fgu = unique(fg,keys='TARGETID')
	print(str(len(fgu))+' unique randoms')
	aa = np.chararray(len(fgu),unicode=True,itemsize=100)
	aa[:] = str(tile)
	fgu['TILE'] = aa
	
	for i in range(1,len(fafls0)):
		fah = fitsio.read_header(fafls0[i])
		tile = fah['TILEID']
		w = exps['TILEID'] == fah['TILEID']
		#if exps[w]['EPOCH'][0] == srun:


		if len(exps[w]) > 1:
				return 'NEED to deal with multiple exposures of same tile'
		if len(exps[w]) > 0:
			expid = exps[w]['EXPID'][0]     
			ep = exps[w]['EPOCH'][0]
			fmap = fitsio.read(e2ein+'run/quicksurvey/'+programf+'/'+str(ep)+'/fiberassign/fibermap-'+str(expid).zfill(8)+'.fits')
			#fmap['FIBERSTATUS'] = 0
			#print('set fiberstatus all to 0; fix this once propagated to zcat')

			wloc = fmap['FIBERSTATUS'] == 0
			gloc = fmap[wloc]['LOCATION']
			fa = Table.read(fafls0[i],hdu='FAVAIL')
			wg = np.isin(fa['LOCATION'],gloc)
			fg = fa[wg]
			fgun = unique(fg,keys='TARGETID')
			aa = np.chararray(len(fgun),unicode=True,itemsize=100)
			aa[:] = str(tile)
			fgun['TILE'] = aa

			#fgun['TILE'] = str(tile)
			#print(len(fg),len(gloc))
			fv = vstack([fgu,fgun])
			#print(len(fv))
			fgo = fgu
			fgu = unique(fv,keys='TARGETID')
			#fguc = setdiff((fgun,fgu))
			dids = np.isin(fgun['TARGETID'],fgo['TARGETID']) #get the rows with target IDs that were duplicates in the new file
			didsc = np.isin(fgu['TARGETID'],fgun['TARGETID'][dids]) #get the row in the concatenated table that had dup IDs
			aa = np.chararray(len(fgu['TILE']),unicode=True,itemsize=20)
			aa[:] = '-'+str(tile)
			#rint(aa)
			ms = np.core.defchararray.add(fgu['TILE'][didsc],aa[didsc])
			print(ms)
			fgu['TILE'][didsc] = ms #add the tile info
			print(str(len(fgu))+' unique randoms')
		#else:
		#       print(str(tile)+' not observed in assigned epoch')      
	print(np.unique(fgu['TILE']))
	#return('ended test')
	print('run '+str(srun) +' done')
	for run in range(srun+1,srun+nrun):
		dirr =  e2eout+program+'/randoms/'+str(run)+'/'
		faflsr = glob.glob(dirr+'fba-*.fits')
		for i in range(0,len(faflsr)):
			fah = fitsio.read_header(faflsr[i])
			tile = fah['TILEID']
			w = exps['TILEID'] == fah['TILEID']
			#if exps[w]['EPOCH'][0] == run:

	
			if len(exps[w]) > 1:
				return 'NEED to deal with multiple exposures of same tile'
			if len(exps[w]) > 0:
				expid = exps[w]['EXPID'][0]     
				ep = exps[w]['EPOCH'][0]
				fmap = fitsio.read(e2ein+'run/quicksurvey/'+programf+'/'+str(ep)+'/fiberassign/fibermap-'+str(expid).zfill(8)+'.fits')
				#fmap['FIBERSTATUS'] = 0
				#print('set fiberstatus all to 0; fix this once propagated to zcat')

				wloc = fmap['FIBERSTATUS'] == 0
				gloc = fmap[wloc]['LOCATION']
				fa = Table.read(faflsr[i],hdu='FAVAIL')
				wg = np.isin(fa['LOCATION'],gloc)
				fg = fa[wg]
				#print(len(fg),len(gloc))
				fgun = unique(fg,keys='TARGETID')
				aa = np.chararray(len(fgun),unicode=True,itemsize=100)
				aa[:] = str(tile)
				fgun['TILE'] = aa

				fv = vstack([fgu,fgun])
				#print(len(fv))
				fgo = fgu
				fgu = unique(fv,keys='TARGETID')
				dids = np.isin(fgun['TARGETID'],fgo['TARGETID']) #get the rows with target IDs that were duplicates in the new file
				didsc = np.isin(fgu['TARGETID'],fgun['TARGETID'][dids]) #get the row in the concatenated table that had dup IDs
				aa = np.chararray(len(fgu['TILE']),unicode=True,itemsize=20)
				aa[:] = '-'+str(tile)
				#rint(aa)
				ms = np.core.defchararray.add(fgu['TILE'][didsc],aa[didsc])
				print(ms)
				fgu['TILE'][didsc] = ms #add the tile info

				print(str(len(fgu))+' unique randoms')
				#else:
				#       print(str(tile)+' not observed in assigned epoch')      

		print(np.unique(fgu['TILE']))
		print('run '+str(run) +' done')
	fgu['PROGRAM'] = program
	if program == 'gray':
		program = 'dark'
		for run in range(srun,srun+nrun):
			dirr =  e2eout+program+'/randoms/'+str(run)+'/'
			faflsr = glob.glob(dirr+'fba-*.fits')
			for i in range(0,len(faflsr)):
				fah = fitsio.read_header(faflsr[i])
				tile = fah['TILEID']
				w = exps['TILEID'] == fah['TILEID']
				#if exps[w]['EPOCH'][0] == run:

	
				if len(exps[w]) > 1:
					return 'NEED to deal with multiple exposures of same tile'
				if len(exps[w]) > 0:
					expid = exps[w]['EXPID'][0]     
					ep = exps[w]['EPOCH'][0]
					fmap = fitsio.read(e2ein+'run/quicksurvey/'+program+'/'+str(ep)+'/fiberassign/fibermap-'+str(expid).zfill(8)+'.fits')
					#fmap['FIBERSTATUS'] = 0
					#print('set fiberstatus all to 0; fix this once propagated to zcat')

					wloc = fmap['FIBERSTATUS'] == 0
					gloc = fmap[wloc]['LOCATION']
					fa = Table.read(faflsr[i],hdu='FAVAIL')
					wg = np.isin(fa['LOCATION'],gloc)
					fg = fa[wg]
					#print(len(fg),len(gloc))
					fgun = unique(fg,keys='TARGETID')
					aa = np.chararray(len(fgun),unicode=True,itemsize=100)
					aa[:] = str(tile)
					fgun['TILE'] = aa
					fgun['PROGRAM'] = program

					fv = vstack([fgu,fgun])
					#print(len(fv))
					fgo = fgu
					fgu = unique(fv,keys='TARGETID')
					dids = np.isin(fgun['TARGETID'],fgo['TARGETID']) #get the rows with target IDs that were duplicates in the new file
					didsc = np.isin(fgu['TARGETID'],fgun['TARGETID'][dids]) #get the row in the concatenated table that had dup IDs
					aa = np.chararray(len(fgu['TILE']),unicode=True,itemsize=20)
					aa[:] = '-'+str(tile)
					#rint(aa)
					ms = np.core.defchararray.add(fgu['TILE'][didsc],aa[didsc])
					print(ms)
					fgu['TILE'][didsc] = ms #add the tile info

					print(str(len(fgu))+' unique randoms')
					#else:
					#       print(str(tile)+' not observed in assigned epoch')      

			print(np.unique(fgu['TILE']))
			print('run '+str(run) +' done')
	up = np.unique(fgu['PROGRAM'])
	print(up)
	for p in up:
		w = fgu['PROGRAM'] == p
		print(str(len(fgu[w]))+ ' '+p)	
	fgu.write(e2eout+outf,format='fits', overwrite=True)    

def combtargets(srun=0,nrun=7,program='dark'):
	'''
	Catalog of all TARGETIDs from e.g. parent MTL that could have been assigned
	to at least one GOOD fiber in any tile, where GOOD is defined by FIBERSTATUS
	in the FIBERMAP.
	'''

	# Glob fiberassign files in epoch 'srun'and read the first in list. 
	programf = program
	if program == 'gray':
		programf = 'dark'
	dir0   = e2ein+'run/quicksurvey/'+programf+'/'+str(srun)+'/fiberassign/'
	outf   = program+'/targets_oneper.fits'
	fafls0 = glob.glob(dir0+'fiberassign-*.fits')
	tile = fitsio.read_header(fafls0[0])['TILEID']
	

	#exps = fitsio.read(e2ein+'run/survey/complete_exposures_surveysim_fix.fits')

	# Exposure file complete with epoch that a given tile was completed in. 
	exps = fitsio.read(e2ein+'run/quicksurvey/'+programf+'/epochs-'+programf+'.fits')

	if program == 'dark':
		we = exps['PROGRAM'] == b'DARK'
		exps = exps[we]


	if program == 'gray':
		we = exps['PROGRAM'] == b'GRAY'
		exps = exps[we]

	print(len(exps))
	# Exposure info. for this tile. 
	w = exps['TILEID'] == tile
	i = 1
	while len(exps[w]) == 0 and i < len(fafls0):
		print(i,len(fafls0))
		tile = fitsio.read_header(fafls0[i])['TILEID']
		w = exps['TILEID'] == tile
		i += 1
	
	while i == len(fafls0):
		srun += 1
		dir0   = e2ein+'run/quicksurvey/'+programf+'/'+str(srun)+'/fiberassign/'
		fafls0 = glob.glob(dir0+'fiberassign-*.fits')
		tile = fitsio.read_header(fafls0[0])['TILEID']
		w = exps['TILEID'] == tile
		i = 1
		while len(exps[w]) == 0 and i < len(fafls0):
			tile = fitsio.read_header(fafls0[i])['TILEID']
			w = exps['TILEID'] == tile
			i += 1

	if len(exps[w]) > 1:
			return 'NEED to deal with multiple exposures of same tile'

	# TO DO:  What happens if an assigned tile was not completed?
	#         This is caught for the other tiles below. 

	#if exps[w]['EPOCH'][0] == srun:#[0]:
	#       pass
	#else:
	#       return 'first tile was not observed in assigned epoch, fix code'

	# EXPID & completion EPOCH for this first assignment. 
	expid = exps[w]['EXPID'][0]     
	ep    = exps[w]['EPOCH'][0]

	#if expid < 100:
	#       zer = '000000'

	# Find the fibermap for this first assigned tile using completion epoch.
	fmap = fitsio.read(e2ein+'/run/quicksurvey/'+programf+'/'+str(ep)+'/fiberassign/fibermap-'+str(expid).zfill(8)+'.fits')
	#fmap['FIBERSTATUS'] = 0
	#print('set fiberstatus all to 0; fix this once propagated to zcat')

	# Good fibers. 
	wloc = fmap['FIBERSTATUS'] == 0
	gloc = fmap[wloc]['LOCATION']

	# Now get the FAVIL info. for the locations with good fibers. 
	fa = Table.read(fafls0[0],hdu='FAVAIL')
	wg = np.isin(fa['LOCATION'],gloc)
	fg = fa[wg]
	fgu = unique(fg,keys='TARGETID')
	print(str(len(fgu))+' unique targets')
	aa = np.chararray(len(fgu),unicode=True,itemsize=100)
	aa[:] = str(tile)
	fgu['TILE'] = aa

	# Stack on the other fiberassign files. 
	for i in range(1,len(fafls0)):
		fah = fitsio.read_header(fafls0[i])
		tile = fah['TILEID']
		w = exps['TILEID'] == fah['TILEID']

		# The assigned tile was completed. 
		if len(exps[w]) > 0:
			#if exps[w]['EPOCH'][0] == srun:

			if len(exps[w]) > 1:
					return 'NEED to deal with multiple exposures of same tile'
			expid = exps[w]['EXPID'][0]     
			ep = exps[w]['EPOCH'][0]
			fmap = fitsio.read(e2ein+'run/quicksurvey/'+programf+'/'+str(ep)+'/fiberassign/fibermap-'+str(expid).zfill(8)+'.fits')
			#fmap['FIBERSTATUS'] = 0
			#print('set fiberstatus all to 0; fix this once propagated to zcat')

			wloc = fmap['FIBERSTATUS'] == 0
			gloc = fmap[wloc]['LOCATION']
			fa = Table.read(fafls0[i],hdu='FAVAIL')
			wg = np.isin(fa['LOCATION'],gloc)
			fg = fa[wg]
			fgun = unique(fg,keys='TARGETID')
			aa = np.chararray(len(fgun),unicode=True,itemsize=100)
			aa[:] = str(tile)
			fgun['TILE'] = aa

			#fgun['TILE'] = str(tile)
			#print(len(fg),len(gloc))
			fv = vstack([fgu,fgun])
			#print(len(fv))
			fgo = fgu
			fgu = unique(fv,keys='TARGETID')
			#fguc = setdiff((fgun,fgu))
			dids = np.isin(fgun['TARGETID'],fgo['TARGETID']) #get the rows with target IDs that were duplicates in the new file
			didsc = np.isin(fgu['TARGETID'],fgun['TARGETID'][dids]) #get the row in the concatenated table that had dup IDs
			aa = np.chararray(len(fgu['TILE']),unicode=True,itemsize=20)
			aa[:] = '-'+str(tile)
			#rint(aa)
			ms = np.core.defchararray.add(fgu['TILE'][didsc],aa[didsc])
			print(ms)
			fgu['TILE'][didsc] = ms #add the tile info
			print(str(len(fgu))+' unique targets')
			#else:
			#       print(str(tile)+' not observed in assigned epoch')      
	print(np.unique(fgu['TILE']))
	#return('ended test')
	print('run '+str(srun) +' done')

	for run in range(srun+1,srun+nrun):
		dirr =  e2ein+'run/quicksurvey/'+programf+'/'+str(run)+'/fiberassign/'
		faflsr = glob.glob(dirr+'fiberassign-*.fits')
		for i in range(0,len(faflsr)):
			fah = fitsio.read_header(faflsr[i])
			tile = fah['TILEID']
			w = exps['TILEID'] == fah['TILEID']
			if len(exps[w]) > 0:
				#if exps[w]['EPOCH'][0] == run:


				if len(exps[w]) > 1:
						return 'NEED to deal with multiple exposures of same tile'
				expid = exps[w]['EXPID'][0]     
				ep = exps[w]['EPOCH'][0]
				fmap = fitsio.read(e2ein+'run/quicksurvey/'+programf+'/'+str(ep)+'/fiberassign/fibermap-'+str(expid).zfill(8)+'.fits')
				#fmap['FIBERSTATUS'] = 0
				#print('set fiberstatus all to 0; fix this once propagated to zcat')

				wloc = fmap['FIBERSTATUS'] == 0
				gloc = fmap[wloc]['LOCATION']
				fa = Table.read(faflsr[i],hdu='FAVAIL')
				wg = np.isin(fa['LOCATION'],gloc)
				fg = fa[wg]
				#print(len(fg),len(gloc))
				fgun = unique(fg,keys='TARGETID')
				aa = np.chararray(len(fgun),unicode=True,itemsize=100)
				aa[:] = str(tile)
				fgun['TILE'] = aa

				fv = vstack([fgu,fgun])
				#print(len(fv))
				fgo = fgu
				fgu = unique(fv,keys='TARGETID')
				dids = np.isin(fgun['TARGETID'],fgo['TARGETID']) #get the rows with target IDs that were duplicates in the new file
				didsc = np.isin(fgu['TARGETID'],fgun['TARGETID'][dids]) #get the row in the concatenated table that had dup IDs
				aa = np.chararray(len(fgu['TILE']),unicode=True,itemsize=20)
				aa[:] = '-'+str(tile)
				#rint(aa)
				ms = np.core.defchararray.add(fgu['TILE'][didsc],aa[didsc])
				print(ms)
				fgu['TILE'][didsc] = ms #add the tile info

				print(str(len(fgu))+' unique targets')
				#else:
				#       print(str(tile)+' not observed in assigned epoch')      

		print(np.unique(fgu['TILE']))
		print('run '+str(run) +' done')
	fgu['PROGRAM'] = program
	if program == 'gray':
		program = 'dark'
		exps = fitsio.read(e2ein+'run/quicksurvey/'+programf+'/epochs-'+programf+'.fits')

		if program == 'dark':
			we = exps['PROGRAM'] == b'DARK'
			exps = exps[we]

		for run in range(srun,srun+nrun):

			dirr =  e2ein+'run/quicksurvey/'+program+'/'+str(run)+'/fiberassign/'
			faflsr = glob.glob(dirr+'fiberassign-*.fits')
			for i in range(0,len(faflsr)):
				fah = fitsio.read_header(faflsr[i])
				tile = fah['TILEID']
				w = exps['TILEID'] == fah['TILEID']
				if len(exps[w]) > 0:
					#if exps[w]['EPOCH'][0] == run:


					if len(exps[w]) > 1:
							return 'NEED to deal with multiple exposures of same tile'
					expid = exps[w]['EXPID'][0]     
					ep = exps[w]['EPOCH'][0]
					fmap = fitsio.read(e2ein+'run/quicksurvey/'+program+'/'+str(ep)+'/fiberassign/fibermap-'+str(expid).zfill(8)+'.fits')
					#fmap['FIBERSTATUS'] = 0
					#print('set fiberstatus all to 0; fix this once propagated to zcat')

					wloc = fmap['FIBERSTATUS'] == 0
					gloc = fmap[wloc]['LOCATION']
					fa = Table.read(faflsr[i],hdu='FAVAIL')
					wg = np.isin(fa['LOCATION'],gloc)
					fg = fa[wg]
					#print(len(fg),len(gloc))
					fgun = unique(fg,keys='TARGETID')
					aa = np.chararray(len(fgun),unicode=True,itemsize=100)
					aa[:] = str(tile)
					fgun['TILE'] = aa
					fgun['PROGRAM'] = program

					fv = vstack([fgu,fgun])
					#print(len(fv))
					fgo = fgu
					fgu = unique(fv,keys='TARGETID')
					dids = np.isin(fgun['TARGETID'],fgo['TARGETID']) #get the rows with target IDs that were duplicates in the new file
					didsc = np.isin(fgu['TARGETID'],fgun['TARGETID'][dids]) #get the row in the concatenated table that had dup IDs
					aa = np.chararray(len(fgu['TILE']),unicode=True,itemsize=20)
					aa[:] = '-'+str(tile)
					#rint(aa)
					ms = np.core.defchararray.add(fgu['TILE'][didsc],aa[didsc])
					print(ms)
					fgu['TILE'][didsc] = ms #add the tile info

					print(str(len(fgu))+' unique targets')
					#else:
					#       print(str(tile)+' not observed in assigned epoch')      
	up = np.unique(fgu['PROGRAM'])
	print(up)
	for p in up:
		w = fgu['PROGRAM'] == p
		print(str(len(fgu[w]))+ ' '+p)	

	fgu.write(e2eout+outf,format='fits', overwrite=True)    

def matchzcatmtl(srun,nrun,program='dark'):
        '''
        Read the last mtl and (largest) zcat in the survey.  Join on TARGETID keeping all entries
        in mtl.
        '''
        
        outf = program+'/mtlzcat'+program+'.fits'
        rmax = srun+nrun-1

        mtl  = Table.read(e2ein+'run/quicksurvey/'+program+'/'+str(rmax)+'/mtl-'+program+'.fits')
        zc   = Table.read(e2ein+'run/quicksurvey/'+program+'/'+str(rmax)+'/zcat-'+program+'.fits')

        # BUG:  Table names should be mtl, zcat?
        mtlj = join(mtl,zc,keys=['TARGETID'], table_names=['mtl','zcat'], join_type='left')

        #       for i in range(srun+1,srun+nrun):
        #               mtl = Table.read(e2ein+'run/quicksurvey/'+str(i)+'/mtl.fits')
        #               zc = Table.read(e2ein+'run/quicksurvey/'+str(i)+'/zcat.fits')
        #               mtlji = join(zc,mtl,keys=['TARGETID'],table_names=['zcat','mtl'])
        #               mtlj = vstack([mtlj,mtlji])

        w    = mtlj['ZWARN'] == 0

        print('number of obs, number of good redshifts:')
        print(len(mtlj),len(mtlj[w]))   

        mtlj.write(e2eout+outf,format='fits', overwrite=True)   

def matchzcattar(program='dark',rmax=6):
        '''

        '''
        
        outf = program+'/targets_oneper_jmtl_jzcat.fits'
        programf = program
        if program == 'gray':
        	programf = 'dark'
        mtl  = Table.read(e2eout+program+'/targets_oneper_jmtl.fits')
        zc   = Table.read(e2ein+'run/quicksurvey/'+programf+'/'+str(rmax)+'/zcat-'+programf+'.fits')

        mtlj = join(mtl,zc,keys=['TARGETID'],table_names=['mtl', 'zcat'],join_type='left')

        w    = mtlj['ZWARN'] == 0

        print('number of targets, number of good redshifts:')
        print(len(mtlj),len(mtlj[w]))   

        mtlj.write(e2eout+outf, format='fits', overwrite=True)
        
def matchzcattar_nofa(program='dark'):
        '''

        '''

        outf = program+'/targets_oneper_jmtl_jztrue.fits'
        programf = program
        if program == 'gray':
        	programf = 'dark'
        mtl  = Table.read(e2eout+program+'/targets_oneper_jmtl.fits')
        
        zc   = Table.read(e2ein+'run/targets/truth-'+programf+'.fits') 
        zc['Z'] = zc['TRUEZ']
        
        mtlj = join(mtl,zc,keys=['TARGETID'],table_names=['mtl', 'zcat'],join_type='left')
        
        w    = (mtlj['Z'] > 0) & (mtlj['Z'] < 5)
        print('number of targets, number of good redshifts:')
        print(len(mtlj),len(mtlj[w]))   

        mtlj.write(e2eout+outf, format='fits', overwrite=True)



def plotcompvsntile(type,program='dark'):
	'''

	'''        
	dz  = fitsio.read(e2eout+ program+'/'+type+'_oneper_full.dat.fits')
	rz  = fitsio.read(e2eout+ program+'/'+type+'_oneper_full.ran.fits')
	normt = len(rz)/len(dz)

	ntl = np.unique(rz['NTILE'])
	zfl = []
	zftl = []
	ttl = []

	for nt in ntl:
		w    = dz['NTILE'] == nt
		ntar = len(dz[w])

		# (dz['NUMOBS_MORE_mtl'] == 0)
		wz   = w & (dz['ZWARN'] == 0)
		nz   = len(dz[wz])
		wr = rz['NTILE'] == nt
		nran = len(rz[wr])
		print(nt,nz,ntar,nran)
		if nran > 100:
			zfl.append(nz/nran*normt)
			zftl.append(nz/ntar)
			ttl.append(nt)
	plt.plot(ttl,zfl,'k-',label='compared to randoms')
	plt.plot(ttl,zftl,'r-',label='compared to all targets')
	plt.xlabel('NTILES')
	plt.ylabel('N good z/N random for '+type)
	plt.legend()

	plt.show()      




def plotzprobvsntile(program='dark',type=0):
        '''
        Defaults to LRG.
        '''        
        dz  = fitsio.read(e2eout+program+'/targets_oneper_jmtl_jzcat.fits')

        # Shouldn't need this when setting target bit:  & (dz['NUMOBS_MORE_mtl'] > -1)
        wt  = (dz['DESI_TARGET'] & 2**type > 0)
        dz  = dz[wt]
        ntl = np.unique(dz['NTILE'])
        zfl = []

        for nt in ntl:
                w    = dz['NTILE'] == nt
                ntar = len(dz[w])

                # (dz['NUMOBS_MORE_mtl'] == 0)
                wz   = w & (dz['ZWARN'] == 0)
                nz   = len(dz[wz])
                print(nt,nz,ntar)
                zfl.append(nz/ntar)
        plt.plot(ntl,zfl,'k-')
        plt.xlabel('NTILES')

        if type == 0:
                plt.ylabel('fraction of LRG targets with good z')

        if type == 1:
                plt.ylabel('fraction of ELG targets with good z')

        if type == 60:
                plt.ylabel('fraction of BGS targets with good z')


        plt.show()      

def mkzprobvsntiledic(type,program='dark'):
	dz = fitsio.read(e2eout+ program+'/'+type+'_oneper_full.dat.fits')
	dr = fitsio.read(e2eout+ program+'/'+type+'_oneper_full.ran.fits')
	ntl = np.unique(dr['NTILE'])
	zfl = []
	for nt in ntl:
		w = dz['NTILE'] == nt
		ntar = len(dz[w])
		if ntar > 0:
			wz = w &  (dz['ZWARN'] == 0)
			nz = len(dz[wz])
			print(nt,nz,ntar)
			zfl.append((nt,nz/ntar))
		else:
			zfl.append((nt,1.))
			
	return dict(zfl)

def plotcompdr_full(type,program='dark'):
        r = fitsio.read(e2eout+program+'/'+type+'_oneper_full.ran.fits')
        rarw = r['RA']
        wr = rarw > 180
        rarw[wr] -= 360
        plt.plot(rarw,r['DEC'],'k,',label='randoms')

        d = fitsio.read(e2eout+program+'/'+type+'_oneper_full.dat.fits')
        w = d['ZWARN'] == 0
        dw = d[w]
        radw = dw['RA']
        wr = radw > 180
        radw[wr] -= 360
        #wt = (dw['DESI_TARGET'] & 2**0 > 0) & (dw['DESI_TARGET'] & 2**1 == 0) #select LRG targets that are not ELGs
        plt.plot(radw,dw['DEC'],'r,',label='data')
        plt.xlabel('RA')
        plt.ylabel('DEC')
        plt.legend()
        plt.title('e2e one per cent survey '+type)
        plt.show()

def plotznz_nt(program='dark'):

        d = fitsio.read(e2eout+program+'/tarzcat'+program+'.fits')
        w = (d['DESI_TARGET'] & 2**0 > 0) & (d['DESI_TARGET'] & 2**1 == 0) #select LRG targets that are not ELGs
        dw = d[w]
        radw = dw['RA']
        wr = radw > 180
        radw[wr] -= 360
        ntl = np.unique(dw['NTILE'])
        for nt in ntl:
                wt = dw['NTILE'] == nt
                wz = wt & (dw['ZWARN']==0)
                plt.plot(radw[wt],dw['DEC'][wt],'k,',label='all targets')
                plt.plot(radw[wz],dw['DEC'][wz],'r,',label='all good z')
                plt.title('regions with '+str(nt)+' overlapping tiles')
                plt.xlabel('RA')
                plt.ylabel('DEC')
                plt.legend()
                plt.show()


def comphistNT(program='dark'):
        r = fitsio.read(e2eout+program+'/randoms_oneper_jmtl.fits')
        plt.hist(r['NTILE'],normed=True,histtype='step',color='k',label='randoms',bins=8,range=(0.5,8.5))
        t = fitsio.read(e2eout+program+'/targets_oneper_jmtl_jzcat.fits')
        w = t['NUMOBS_MORE_mtl'] > -1
        t = t[w]
        plt.hist(t['NTILE'],normed=True,histtype='step',color='r',label='targets',bins=8,range=(0.5,8.5))
        plt.legend()
        plt.xlabel('Number of Tiles')
        plt.ylabel('Relative number')
        plt.show()

def plotrntile(program='dark'):
        r = fitsio.read(e2eout+program+'/randoms_oneper_jmtl.fits')
        rarw = r['RA']
        wr = rarw > 180
        rarw[wr] -= 360
        plt.scatter(rarw,r['DEC'],c=r['NTILE'],marker='.',s=.1)

        plt.xlabel('RA')
        plt.ylabel('DEC')
        plt.colorbar()
        plt.title('e2e one per cent survey NTILE for randoms')
        plt.show()

def plottntile(program='dark'):
        r = fitsio.read(e2eout+program+'/targets_oneper_jmtl.fits')
        rarw = r['RA']
        wr = rarw > 180
        rarw[wr] -= 360
        plt.scatter(rarw,r['DEC'],c=r['NTILE'],marker='.',s=.1)

        plt.xlabel('RA')
        plt.ylabel('DEC')
        plt.colorbar()
        plt.title('e2e one per cent survey NTILE for targets')
        plt.show()


def plotzcat_tilecen(pr='dark'):
        if pr == 'dark':
                pt = b'DARK'
        ts = fitsio.read(e2ein+'run/survey/tiles/des.fits')     
        wp = ts['PROGRAM'] == pt
        ts = ts[wp]
        df=e2eout+'mtlzcat'+pr+'.fits'
        d = fitsio.read(df)
        w = d['ZWARN'] == 0
        dw = d[w]
        raw = dw['RA']
        wr = raw > 300
        raw[wr] -= 360
        plt.plot(raw,dw['DEC'],'k,',label='dark time good z',zorder=0)
        
        exps = fitsio.read(e2ein+'run/survey/complete_exposures_surveysim_fix.fits')
        s = 0
        for tile in exps['TILEID']:
                wt = ts['TILEID'] == tile
                if s == 0:
                        plt.plot(ts[wt]['RA'],ts[wt]['DEC'],'ro',label='completed dark time tile centers')
                else:
                        plt.plot(ts[wt]['RA'],ts[wt]['DEC'],'ro')
        plt.show()              

def testfavail(tile,epoch=6,program='dark'):
        mtl = Table.read(e2ein+'run/quicksurvey/'+program+'/'+str(epoch)+'/mtl-'+program+'.fits')       
        tilef = Table.read(e2ein+'run/quicksurvey/'+program+'/'+str(epoch)+'/fiberassign/fiberassign-'+str(tile).zfill(6)+'.fits',hdu='FAVAIL')
        tj = join(tilef,mtl,keys=['TARGETID'],join_type='left')
        print(np.unique(tj['NUMOBS_MORE']))
                
def matchran(program='dark'):
        faran = Table.read(e2eout+program+'/randoms/randoms_oneper.fits')
        faran['NTILE'] = np.char.count(faran['TILE'],'-')
        faran['NTILE'] += 1
        print(max(faran['NTILE']))
        mtlran = Table.read(e2eout+program+'/randoms_mtl_cuttod.fits')
        jran = join(faran,mtlran,keys=['TARGETID'],table_names=['fa', 'mtl'])
        print(len(jran),len(faran),len(mtlran))
        jran.write(e2eout+program+'/randoms_oneper_jmtl.fits',format='fits', overwrite=True)

def matchtar(program='dark',rmax=6):
        '''
        Read targets in the (one-percent) geometry - available to a GOOD fiber - and join to mtl.
        Tracking NTILE:  how many tiles a target COULD have been observed in.  
        '''
        faran           = Table.read(e2eout+program+'/targets_oneper.fits')
        faran['NTILE']  = np.char.count(faran['TILE'],'-')

        # Counting was - based. 
        faran['NTILE'] += 1

        print(max(faran['NTILE']))
        programf = program
        if program == 'gray':
        	programf = 'dark'

        mtlran          = Table.read(e2ein+'run/quicksurvey/'+programf+'/'+str(rmax)+'/mtl-'+programf+'.fits')
        jran            = join(faran,mtlran,keys=['TARGETID'],table_names=['fa', 'mtl'])

        print(len(jran),len(faran),len(mtlran))

        jran.write(e2eout+program+'/targets_oneper_jmtl.fits',format='fits', overwrite=True)

def mkfulldat(type,program,bits,truez=False):
	'''
	take targets, cut them to particular target type and mask for particular target type 
	program is dark,gray, or bright
	type is 'LRG', 'QSO', 'ELG', or 'BGS'
	bits are imaging mask bits to apply
	'''    
	if type == 'LRG':
		#bits = elgandlrgbits 
		tb = 0
	if type == 'ELG':
		#bits = elgandlrgbits 
		tb = 1
	if type == 'QSO':
		tb = 2    
	if type == 'BGS':
		tb = 60    

	print(truez)
	if truez:
		print('running truez full dat')
		tarf = Table.read(e2eout+ program+'/targets_oneper_jmtl_jztrue.fits')
		outf = e2eout+ program+'/'+type+'_oneperztrue_full.dat.fits'
	else:
		tarf = Table.read(e2eout+ program+'/targets_oneper_jmtl_jzcat.fits')
		outf = e2eout+ program+'/'+type+'_oneper_full.dat.fits'
	tarf = cutphotmask(tarf,bits) 
	wt = tarf['DESI_TARGET'] & 2**tb > 0
	tt = tarf[wt]
	
	tt.write(outf,format='fits', overwrite=True)
    
def mkclusdat(type,program,truez=False):
	'''
	take full catalog, cut to ra,dec,z add any weight
	program is dark,gray, or bright
	type is 'LRG', 'QSO', 'ELG', or 'BGS'

	'''    
	if truez:
		ff = Table.read(e2eout+ program+'/'+type+'_oneperztrue_full.dat.fits')
		outf = e2eout+ program+'/'+type+'_oneperztrue_clus.dat.fits'

	else:
		ff = Table.read(e2eout+ program+'/'+type+'_oneper_full.dat.fits')
		outf = e2eout+ program+'/'+type+'_oneper_clus.dat.fits'
		wz = ff['ZWARN'] == 0
		ff = ff[wz]

	ff.keep_columns(['RA','DEC','Z'])
	ff['WEIGHT'] = np.ones(len(ff))
	ff.write(outf,format='fits', overwrite=True)
    
       
    
def mkfullran(type,program,bits):
    '''
    take randoms, mask for particular target type 
    program is dark,gray, or bright
    type is 'LRG', 'QSO', 'ELG', or 'BGS'
    bits are imaging mask bits to apply

    '''    
        
    tarf = Table.read(e2eout+program+'/randoms_oneper_jmtl.fits')
    tarf = cutphotmask(tarf,bits) 
    outf = e2eout+ program+'/'+type+'_oneper_full.ran.fits'
    tarf.write(outf,format='fits', overwrite=True)

def mkclusran(type,program,truez=False):
    from random import random
    '''
    take full catalog, cut to ra,dec,z add any weight
    assign redshifts by randomly sampling data clustering
    '''    

#     if type == 'LRG':
#         #bits = elgandlrgbits 
#         tb = 0
#     if type == 'QSO':
#     	tb = 2    
#     if type == 'BGS':
#     	tb = 60    

    ff = Table.read(e2eout+ program+'/'+type+'_oneper_full.ran.fits')
    
    if truez:
    	ffd = Table.read(e2eout+ program+'/'+type+'_oneperztrue_clus.dat.fits')
    	outf = e2eout+ program+'/'+type+'_oneperztrue_clus.ran.fits'
    else:
    	ffd = Table.read(e2eout+ program+'/'+type+'_oneper_clus.dat.fits')
    	outf = e2eout+ program+'/'+type+'_oneper_clus.ran.fits'
    	zeffdic = mkzprobvsntiledic(type,program=program)
    #ff['WEIGHT'] = zeffdic[ff['NTILE']]
    ff['WEIGHT']= np.ones(len(ff))
    ff['Z'] = np.zeros(len(ff))
    nd = len(ffd)
    for i in range(0,len(ff)):
    	ind = int(random()*nd)
    	ff['Z'][i] = ffd['Z'][ind]
    	if truez == False:
    		ff['WEIGHT'][i] = zeffdic[ff['NTILE'][i]]
    		ff['WEIGHT'][i] *= ffd['WEIGHT'][ind]

    ff.keep_columns(['RA','DEC','Z','WEIGHT'])
    print(np.min(ff['WEIGHT']),np.max(ff['WEIGHT']))
    ff.write(outf,format='fits', overwrite=True)
       

def cutphotmask(aa,bits):
	keep = (aa['NOBS_G']>0) & (aa['NOBS_R']>0) & (aa['NOBS_Z']>0)
	for biti in bits:
		keep &= ((aa['MASKBITS'] & 2**biti)==0)
	aa = aa[keep]
	print(str(len(aa)) +' after imaging veto' )
	return aa
        


def mke2etiles(run,program='dark'):
        dirout=e2eout+program
        if program == 'gray':
        	program = 'dark'
        fout = dirout+'/e2etiles_run'+str(run)+'.fits'
        fafiles = glob.glob(e2ein+'run/quicksurvey/'+program+'/'+str(run)+'/fiberassign/*')
        atl = fitsio.read(e2ein+'run/survey/tiles/des.fits')
        tls = []
        for i in range(0,len(fafiles)):
                tl = fafiles[i].split('-')[-1].strip('.fits')
                tl = int(tl)
                print(tl)
                tls.append(tl)
        w = np.isin(atl['TILEID'],tls)
        rtl = atl[w]
        print(run,len(rtl))
        print('writing to '+fout)
        fitsio.write(fout,rtl,clobber=True)
        
                        

def gathertargets(type):
        fns      = glob.glob(targroot+'*.fits')
        keys = ['RA', 'DEC', 'BRICKNAME','MORPHTYPE','DCHISQ','FLUX_G', 'FLUX_R', 'FLUX_Z','MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z','NOBS_G', 'NOBS_R', 'NOBS_Z','PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z', 'GALDEPTH_G', 'GALDEPTH_R',\
        'GALDEPTH_Z','FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z',\
        'MASKBITS', 'EBV', 'PHOTSYS','TARGETID','DESI_TARGET']
        #put information together, takes a couple of minutes
        ncat     = len(fns)
        mydict   = {}
        for key in keys:
                mydict[key] = []
        if type == 'ELG':
                bit = 1 #target bit for ELGs
        if type == 'LRG':
                bit = 0
        if type == 'QSO':
                bit = 2
        for i in range(0,ncat):
                data = fitsio.read(fns[i],columns=keys)
                data = data[(data['DESI_TARGET'] & 2**bit)>0]
                for key in keys:
                        mydict[key] += data[key].tolist()
                print(i)        
        outf = tardir+type+'allDR8targinfo.fits'
        collist = []
        for key in keys:
                fmt = fits.open(fns[0])[1].columns[key].format
                collist.append(fits.Column(name=key,format=fmt,array=mydict[key]))
                print(key)
        hdu  = fits.BinTableHDU.from_columns(fits.ColDefs(collist))
        hdu.writeto(outf,overwrite=True)
        print('wrote to '+outf)
        
if __name__ == '__main__':
        #combran()      
        #matchran()
        #matchzcatmtl(0,7)
        #plotcompdr()
        #plotrntile()
        #testfavail(47693)
        #testfavail(47714)
        #combtargets()
        #matchtar(rmax=0)
        #plottntile()
        #plotrntile()
        #matchzcattar()
        #plotzprobvsntile(type=0)
        #plotzprobvsntile(type=1)
        #comphistNT()
        #plotznz_nt()

        print('\n\n<Done.\n\n')
