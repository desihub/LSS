'''
python functions to do various useful date processing/manipulation
'''
import numpy as np
import fitsio
import glob
import astropy.io.fits as fits
from astropy.table import Table,vstack,unique,join#,setdiff
from matplotlib import pyplot as plt
import desimodel.footprint
import desimodel.focalplane #

ranf = '/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randomsall/randoms-inside-dr8-0.31.0-all.fits' #DR8 imaging randoms file
minisvdir = '/project/projectdirs/desi/users/ajross/catalogs/minisv2/'
tardir = minisvdir+'targets/'
dircat = minisvdir+'LSScats/'
targroot = '/project/projectdirs/desi/target/catalogs/dr8/0.31.1/targets/main/resolve/targets-dr8'
e2ein = '/global/homes/m/mjwilson/desi/survey-validation/svdc-spring2020f-onepercent/'
e2eout = '/project/projectdirs/desi/users/ajross/catalogs/e2eoneper/'


def mkran4fa(N=2e8,fout='random_mtl.fits',dirout=minisvdir+'random/'):
	'''
	cut imaging random file to first N entries and add columns necessary for fiberassignment routines
	'''
	rall = fitsio.read('/project/projectdirs/desi/target/catalogs/dr8/0.31.0/randomsall/randoms-inside-dr8-0.31.0-all.fits',rows=np.arange(N))
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

def combran(srun=0,nrun=7,program='dark'):
	dir0 = '/project/projectdirs/desi/users/ajross/catalogs/e2eoneper/'+program+'/randoms/'+str(srun)+'/'
	outf=program+'/randoms/randoms_darktime.fits'
	fafls0 = glob.glob(dir0+'fba-*.fits')
	fah = fitsio.read_header(fafls0[0])
	tile = fah['TILEID']
	#
	#exps = fitsio.read(e2ein+'run/survey/complete_exposures_surveysim_fix.fits')
	exps = fitsio.read(e2ein+'run/quicksurvey/'+program+'/epochs-'+program+'.fits')
	w = exps['TILEID'] == tile
	if len(exps[w]) > 1:
		return 'NEED to deal with multiple exposures of same tile'
	
	#if exps[w]['EPOCH'][0] == srun:#[0]:
	#	pass
	#else:
	#	return 'first tile was not observed in assigned epoch, fix code'
	expid = exps[w]['EXPID'][0]	
	ep = exps[w]['EPOCH'][0]
	#if expid < 100:
	#	zer = '000000'
	fmap = fitsio.read(e2ein+'run/quicksurvey/'+program+'/'+str(ep)+'/fiberassign/fibermap-'+str(expid).zfill(8)+'.fits')
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
		expid = exps[w]['EXPID'][0]	
		ep = exps[w]['EPOCH'][0]
		fmap = fitsio.read(e2ein+'run/quicksurvey/'+program+'/'+str(ep)+'/fiberassign/fibermap-'+str(expid).zfill(8)+'.fits')
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
		#	print(str(tile)+' not observed in assigned epoch')	
	print(np.unique(fgu['TILE']))
	#return('ended test')
	print('run '+str(srun) +' done')
	for run in range(srun+1,srun+nrun):
		dirr = 	'/project/projectdirs/desi/users/ajross/catalogs/e2eoneper/'+program+'/randoms/'+str(run)+'/'
		faflsr = glob.glob(dirr+'fba-*.fits')
		for i in range(0,len(faflsr)):
			fah = fitsio.read_header(faflsr[i])
			tile = fah['TILEID']
			w = exps['TILEID'] == fah['TILEID']
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
			#	print(str(tile)+' not observed in assigned epoch')	

		print(np.unique(fgu['TILE']))
		print('run '+str(run) +' done')
	fgu.write(e2eout+outf,format='fits', overwrite=True)	

def matchzcatmtl(srun,nrun,program='dark'):
	outf=program+'/mtlzcat'+program+'.fits'
	rmax = srun+nrun-1
	mtl = Table.read(e2ein+'run/quicksurvey/'+program+'/'+str(rmax)+'/mtl-'+program+'.fits')
	zc = Table.read(e2ein+'run/quicksurvey/'+program+'/'+str(rmax)+'/zcat-'+program+'.fits')
	mtlj = join(mtl,zc,keys=['TARGETID'],table_names=['zcat','mtl'],join_type='left')

	# 	for i in range(srun+1,srun+nrun):
	# 		mtl = Table.read(e2ein+'run/quicksurvey/'+str(i)+'/mtl.fits')
	# 		zc = Table.read(e2ein+'run/quicksurvey/'+str(i)+'/zcat.fits')
	# 		mtlji = join(zc,mtl,keys=['TARGETID'],table_names=['zcat','mtl'])
	# 		mtlj = vstack([mtlj,mtlji])
	w = mtlj['ZWARN'] == 0
	print('number of obs, number of good redshifts:')
	print(len(mtlj),len(mtlj[w]))	
	mtlj.write(e2eout+outf,format='fits', overwrite=True)	

def plotcompdr(program='dark'):
	r = fitsio.read(e2eout+program+'/randoms/randoms_darktime_jmtl.fits')
	rarw = r['RA']
	wr = rarw > 180
	rarw[wr] -= 360
	plt.plot(rarw,r['DEC'],'k,',label='randoms')

	d = fitsio.read(e2eout+program+'/mtlzcat'+program+'.fits')
	w = d['ZWARN'] == 0
	dw = d[w]
	radw = dw['RA']
	wr = radw > 180
	radw[wr] -= 360
	wt = (dw['DESI_TARGET'] & 2**0 > 0) & (dw['DESI_TARGET'] & 2**1 == 0) #select LRG targets that are not ELGs
	plt.plot(radw[wt],dw[wt]['DEC'],'r,',label='data')
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.legend()
	plt.title('e2e one per cent survey LRGs')
	plt.show()

def plotrntile(program='dark'):
	r = fitsio.read(e2eout+program+'/randoms/randoms_darktime_jmtl.fits')
	rarw = r['RA']
	wr = rarw > 180
	rarw[wr] -= 360
	plt.scatter(rarw,r['DEC'],c=r['NTILE'],marker='.')

	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.colorbar()
	plt.title('e2e one per cent survey NTILE for randoms')
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
	
		
def matchran(program='dark'):
	faran =	Table.read(e2eout+program+'/randoms/randoms_darktime.fits')
	faran['NTILE'] = np.char.count(faran['TILE'],'-')
	faran['NTILE'] += 1
	print(max(faran['NTILE']))
	mtlran = Table.read(e2eout+program+'/randoms/randoms_mtl_cuttod.fits')
	jran = join(faran,mtlran,keys=['TARGETID'])
	print(len(jran),len(faran),len(mtlran))
	jran.write(e2eout+program+'/randoms/randoms_darktime_jmtl.fits',format='fits', overwrite=True)
	

def randomtiles(tilef = minisvdir+'msvtiles.fits'):
	tiles = fitsio.read(tilef)
	rt = fitsio.read(minisvdir+'random/random_mtl.fits')
	print('loaded random file')
	indsa = desimodel.footprint.find_points_in_tiles(tiles,rt['RA'], rt['DEC'])
	print('got indexes')
	for i in range(0,len(indsa)):
		tile = tiles['TILEID']
		fname = minisvdir+'random/tilenofa-'+str(tile)+'.fits'
		inds = indsa[i]
		fitsio.write(fname,rt[inds],clobber=True)
		print('wrote tile '+str(tile))

def randomtilesi(tilef = minisvdir+'msvtiles.fits'):
	tiles = fitsio.read(tilef)
	trad = desimodel.focalplane.get_tile_radius_deg()*1.1 #make 10% greater just in case
	print(trad)
	rt = fitsio.read(minisvdir+'random/random_mtl.fits')
	print('loaded random file')	
	
	for i in range(0,len(tiles)):
		tile = tiles['TILEID'][i]
		fname = minisvdir+'random/tilenofa-'+str(tile)+'.fits'
		tdec = tiles['DEC'][i]
		decmin = tdec - trad
		decmax = tdec + trad
		wdec = (rt['DEC'] > decmin) & (rt['DEC'] < decmax)
		print(len(rt[wdec]))
		inds = desimodel.footprint.find_points_radec(tiles['RA'][i], tdec,rt[wdec]['RA'], rt[wdec]['DEC'])
		print('got indexes')
		fitsio.write(fname,rt[wdec][inds],clobber=True)
		print('wrote tile '+str(tile))

def targtilesi(type,tilef = minisvdir+'msvtiles.fits'):
	tiles = fitsio.read(tilef)
	trad = desimodel.focalplane.get_tile_radius_deg()*1.1 #make 10% greater just in case
	print(trad)
	rt = fitsio.read(tardir+type+'allDR8targinfo.fits')
	print('loaded random file')	
	
	for i in range(0,len(tiles)):
		tile = tiles['TILEID'][i]
		fname = tardir+type+str(tile)+'.fits'
		tdec = tiles['DEC'][i]
		decmin = tdec - trad
		decmax = tdec + trad
		wdec = (rt['DEC'] > decmin) & (rt['DEC'] < decmax)
		print(len(rt[wdec]))
		inds = desimodel.footprint.find_points_radec(tiles['RA'][i], tdec,rt[wdec]['RA'], rt[wdec]['DEC'])
		print('got indexes')
		fitsio.write(fname,rt[wdec][inds],clobber=True)
		print('wrote tile '+str(tile))

def mke2etiles(run,dirout=e2eout,program='dark'):
	fout = dirout+'e2etiles_run'+str(run)+'.fits'
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
	
def mkminisvtilef(dirout=minisvdir,fout='msvtiles.fits'):
	'''
	manually make tile fits file for sv tiles
	'''
	msvtiles = Table()
	msvtiles['TILEID'] = np.array([70000,70001,70002,70003,70004,70005,70006],dtype=int)
	msvtiles['RA'] = np.array([119.,133.,168.,214.75,116.,158.,214.75])
	msvtiles['DEC'] = np.array([50.,26.5,27.6,53.4,20.7,25.,53.4])
	msvtiles['PASS'] = np.zeros(7,dtype=int)
	msvtiles['IN_DESI'] = np.ones(7,dtype=int)
	msvtiles['OBSCONDITIONS'] = np.ones(7,dtype=int)*65535
	pa = []
	for i in range(0,7):
		pa.append(b'DARK')
	msvtiles['PROGRAM'] = np.array(pa,dtype='|S6')
	msvtiles.write(dirout+fout,format='fits', overwrite=True)
	
def plotdatran(type,tile,night):
	df = fitsio.read(dircat+type +str(tile)+'_'+night+'_clustering.dat.fits')
	rf = fitsio.read(dircat+type +str(tile)+'_'+night+'_clustering.ran.fits')
	plt.plot(rf['RA'],rf['DEC'],'k,')		     
	if type == 'LRG':
		pc = 'r'
		pt = 'o'
	if type == 'ELG':
		pc = 'b'
		pt = '*'
	plt.scatter(df['RA'],df['DEC'],s=df['WEIGHT']*3,c=pc,marker=pt)
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.title(type + ' '+tile+' '+night)
	plt.savefig('dataran'+type+tile+night+'.png')
	plt.show()
		

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
	plotrntile()
	
	