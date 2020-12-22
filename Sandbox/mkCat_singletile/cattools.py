'''
python functions to do various useful date processing/manipulation
'''
import numpy as np
import fitsio
import glob
import astropy.io.fits as fits
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt
import desimodel.footprint
import desimodel.focalplane
from random import random


def combspecdata(tile,night,coaddir ):
    #put data from different spectrographs together, one table for fibermap, other for z
    specs = []
    #find out which spectrograph have data
    for si in range(0,10):
        try:
            fitsio.read(coaddir+str(tile)+'/'+night+'/zbest-'+str(si)+'-'+str(tile)+'-'+night+'.fits')
            specs.append(si)
        except:
            print('no spectrograph '+str(si)+ ' on night '+night)
    print('spectrographs with data:')
    print(specs)			
    tspec = Table.read(coaddir+str(tile)+'/'+night+'/zbest-'+str(specs[0])+'-'+str(tile)+'-'+night+'.fits',hdu='ZBEST')
    tf = Table.read(coaddir+str(tile)+'/'+night+'/zbest-'+str(specs[0])+'-'+str(tile)+'-'+night+'.fits',hdu='FIBERMAP')
    for i in range(1,len(specs)):
        tn = Table.read(coaddir+str(tile)+'/'+night+'/zbest-'+str(specs[i])+'-'+str(tile)+'-'+night+'.fits',hdu='ZBEST')
        tnf = Table.read(coaddir+str(tile)+'/'+night+'/zbest-'+str(specs[i])+'-'+str(tile)+'-'+night+'.fits',hdu='FIBERMAP')
        tspec = vstack([tspec,tn])
        tf = vstack([tf,tnf])
    tf = unique(tf,keys=['TARGETID'])
    tf.keep_columns(['TARGETID','LOCATION','FIBERSTATUS','PRIORITY'])
    tspec = join(tspec,tf,keys=['TARGETID'])
    print(len(tspec),len(tf))
    #tspec['LOCATION'] = tf['LOCATION']
    #tspec['FIBERSTATUS'] = tf['FIBERSTATUS']
    #tspec['PRIORITY'] = tf['PRIORITY']
    return tspec

def goodlocdict(tf):
    '''
    Make a dictionary to map between location and priority
    '''
    wloc = tf['FIBERSTATUS'] == 0
    print(str(len(tf[wloc])) + ' locations with FIBERSTATUS 0')
    goodloc = tf[wloc]['LOCATION']
    pdict = dict(zip(tf['LOCATION'], tf['PRIORITY'])) #to be used later for randoms
    return pdict,goodloc

def gettarinfo_type(fadir,tile,goodloc,mtlf,tarbit,tp='CMX_TARGET'):
    #get target info
    tfa = Table.read(fadir+'fba-0'+str(tile)+'.fits',hdu='FAVAIL')
    tft = unique(tfa,keys=['TARGETID'])
   
    wgt = (np.isin(tfa['LOCATION'],goodloc)) 
    print(str(len(np.unique(tfa[wgt]['LOCATION']))) + ' good locations')
    print('comparison of number targets, number of targets with good locations')
    print(len(tfa),len(tfa[wgt]))
    tfa = unique(tfa[wgt],keys=['TARGETID'])
    tt = Table.read(mtlf)
    tt.remove_columns(['Z','ZWARN'])
    wtype = ((tt[tp] & 2**tarbit) > 0)
    tt = tt[wtype]
    tfa = join(tfa,tt,keys=['TARGETID'])
    tft = join(tft,tt,keys=['TARGETID'])
    print(str(len(tfa)) +' unique targets with good locations and  at '+str(len(np.unique(tfa['LOCATION'])))+' unique locations and '+str(len(tft))+ ' total unique targets at '+str(len(np.unique(tft['LOCATION']))) +' unique locations ')

    #Mark targets that actually got assigned fibers
    tfall = Table.read(fadir+'fba-0'+str(tile)+'.fits',hdu='FASSIGN')
    
    tfall.keep_columns(['TARGETID','LOCATION'])
    tfa = join(tfa,tfall,keys=['TARGETID'],join_type='left',table_names = ['', '_ASSIGNED'], uniq_col_name='{col_name}{table_name}')
    wgl = np.isin(tfa['LOCATION_ASSIGNED'],goodloc)
    wtype = ((tfa[tp] & 2**tarbit) > 0)
    wtfa = wgl & wtype
    print('number of assigned fibers at good locations '+str(len(tfa[wtfa])))

    wal = tfa['LOCATION_ASSIGNED']*0 == 0
    print('number of assigned fibers '+str(len(tfa[wal])))
    tfa['LOCATION_ASSIGNED'] = np.zeros(len(tfa),dtype=int)
    tfa['LOCATION_ASSIGNED'][wal] = 1
    wal = tfa['LOCATION_ASSIGNED'] == 1
    print('number of assigned fibers '+str(len(tfa[wal])))

    return tfa

def mkfullran(tile,goodloc,pdict,randir):
    ranf = randir+'fba-0'+str(tile)+'.fits'
    f1 = fitsio.read(ranf)
    f2 = fitsio.read(ranf,ext=2)
    f3 = fitsio.read(ranf,ext=3)

    goodranw = np.isin(f3['LOCATION'],goodloc)
    goodranid = np.unique(f3[goodranw]['TARGETID'])

    t2 = Table.read(ranf,hdu=2)
    tj = Table()
    tj['TARGETID'] = f3[goodranw]['TARGETID']
    tj['LOCATION'] = f3[goodranw]['LOCATION']
    tj['FIBER'] = f3[goodranw]['FIBER']
    tj = unique(tj,keys=['TARGETID'])
    t2.remove_columns(['PRIORITY','OBSCONDITIONS','SUBPRIORITY'])
    rant = join(tj,t2,keys=['TARGETID'],join_type='left')    
    #now match back to randoms with all columns

    tall = Table.read(randir+'tilenofa-'+str(tile)+'.fits')
    tall.remove_columns(['NUMOBS_MORE','PRIORITY','OBSCONDITIONS','SUBPRIORITY','NUMOBS_INIT'])

    ranall = join(rant,tall,keys=['TARGETID'],join_type='left')
    print('number of randoms:')
    print(len(ranall))


    ranall['PRIORITY'] = np.vectorize(pdict.__getitem__)(ranall['LOCATION'])
    return ranall

def mkclusdat(ffd,fcd,zfailmd= 'zwarn',weightmd= 'wloc',maskbits=[]):
    dd = fitsio.read(ffd)	
    
    ddm = cutphotmask(dd,maskbits)
    #print(np.unique(dd['ZWARN']))
    maxp = np.max(dd['PRIORITY'])
    if zfailmd == 'zwarn':
        wfail = (dd['ZWARN'] != 999999) & (dd['ZWARN'] > 0)	
        wg = (ddm['ZWARN'] == 0) 
    loc_fail = dd[wfail]['LOCATION']	
    print(' number of redshift failures:')
    print(len(loc_fail))
# 
    
    nl = countloc(ddm)
# 
    
    ddzg = ddm[wg]
# 
    print('clustering catalog will have '+str(len(ddzg))+ ' objects in it')
# 
    ddclus = Table()
    ddclus['RA'] = ddzg['RA']
    ddclus['DEC'] = ddzg['DEC']
    ddclus['Z'] = ddzg['Z']
    
    if weightmd == 'wloc':
        ddclus['WEIGHT'] = assignweights(ddzg,nl)
# 
    print('minimum,maximum weight')
    print(np.min(ddclus['WEIGHT']),np.max(ddclus['WEIGHT']))	
# 
    ddclus.write(fcd,format='fits',overwrite=True)
    print('write clustering data file to '+fcd)
    return maxp,loc_fail

def mkclusran(ffr,fcr,fcd,maxp,loc_fail,maskbits=[]):
    dr = fitsio.read(ffr)
    drm = cutphotmask(dr,maskbits)
# 
    wpr = drm['PRIORITY'] <= maxp
    wzf = np.isin(drm['LOCATION'],loc_fail)
    wzt = wpr & ~wzf
# 
    drmz = drm[wzt]
    print(str(len(drmz))+' after cutting based on failures and priority')
    rclus = Table()
    rclus['RA'] = drmz['RA']
    rclus['DEC'] = drmz['DEC']
    dd = fitsio.read(fcd) 
    #rclus['Z'] = 0
    #rclus['WEIGHT'] = 1
    zl = []
    wl = []
    ndz = 0
    naz = 0
    for ii in range(0,len(rclus)):
        ind = int(random()*len(dd))
        zr = dd[ind]['Z']
        if zr == 0:
            ndz += 1.
        naz += 1    
        wr = dd[ind]['WEIGHT']
        #rclus[ii]['Z'] = zr
        #rclus[ii]['WEIGHT'] = wr
        zl.append(zr)
        wl.append(wr)
    zl = np.array(zl)
    wl = np.array(wl)
    rclus['Z'] = zl
    rclus['WEIGHT'] = wl
    wz = rclus['Z'] == 0
    print(ndz,naz,len(rclus[wz]))
    rclus.write(fcr,format='fits',overwrite=True)
    print('write clustering random file to '+fcr)



def cutphotmask(aa,bits):
    print(str(len(aa)) +' before imaging veto' )
    keep = (aa['NOBS_G']>0) & (aa['NOBS_R']>0) & (aa['NOBS_Z']>0)
    for biti in bits:
        keep &= ((aa['MASKBITS'] & 2**biti)==0)
    aa = aa[keep]
    print(str(len(aa)) +' after imaging veto' )
    return aa

def countloc(aa):
    locs = aa['LOCATION']
    la = np.max(locs)+1
    nl = np.zeros(la)
    for i in range(0,len(aa)):
        nl[locs[i]] += 1
    return nl

def assignweights(aa,nl):
    wts = np.ones(len(aa))
    for i in range(0,len(aa)):
        loc = aa[i]['LOCATION']
        wts[i] = nl[loc]
    return wts	


def mkran4fa(N=None,fout='random_mtl.fits',dirout=''):
	'''
	cut imaging random file to first N (or take all if N is None) entries and add columns necessary for fiberassignment routines
	'''
	if N is None:	
		rall = fitsio.read('/global/cfs/cdirs/desi/target/catalogs/dr9m/0.44.0/randoms/resolve/randoms-1-0.fits')
	else:
		rall = fitsio.read('/global/cfs/cdirs/desi/target/catalogs/dr9m/0.44.0/randoms/resolve/randoms-1-0.fits',rows=np.arange(N))
	    
	print('read '+str(len(rall))+ ' rows from random file')
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
	print('added columns, writing to '+dirout+fout)
	del rall
	rmtl.write(dirout+fout,format='fits', overwrite=True)

def randomtiles(tilef ):
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

def randomtilesi(tilef ,dirout,ii):
	tiles = fitsio.read(tilef)
	trad = desimodel.focalplane.get_tile_radius_deg()*1.1 #make 10% greater just in case
	print(trad)
	rt = fitsio.read('/global/cfs/cdirs/desi/target/catalogs/dr9m/0.44.0/randoms/resolve/randoms-1-'+str(ii)+'.fits')
	#rt = fitsio.read(minisvdir+'random/random_mtl.fits')
	print('loaded random file')	
	
	for i in range(0,len(tiles)):
		tile = tiles['TILEID'][i]
		fname = dirout+str(ii)+'/tilenofa-'+str(tile)+'.fits'
		tdec = tiles['DEC'][i]
		decmin = tdec - trad
		decmax = tdec + trad
		wdec = (rt['DEC'] > decmin) & (rt['DEC'] < decmax)
		print(len(rt[wdec]))
		inds = desimodel.footprint.find_points_radec(tiles['RA'][i], tdec,rt[wdec]['RA'], rt[wdec]['DEC'])
		print('got indexes')
		#fitsio.write(fname,rt[wdec][inds],clobber=True)
		#print('wrote tile '+str(tile))
		#rmtl = Table.read(fname)
		rtw = rt[wdec][inds]
		rmtl = Table(rtw)
		rmtl['TARGETID'] = np.arange(len(rmtl))
		rmtl['DESI_TARGET'] = np.ones(len(rmtl),dtype=int)*2
		rmtl['SV1_DESI_TARGET'] = np.ones(len(rmtl),dtype=int)*2
		rmtl['NUMOBS_INIT'] = np.zeros(len(rmtl),dtype=int)
		rmtl['NUMOBS_MORE'] = np.ones(len(rmtl),dtype=int)
		rmtl['PRIORITY'] = np.ones(len(rmtl),dtype=int)*3400
		rmtl['OBSCONDITIONS'] = np.ones(len(rmtl),dtype=int)*tiles['OBSCONDITIONS'][i]
		rmtl['SUBPRIORITY'] = np.random.random(len(rmtl))
		print('added columns, writing to '+fname)
		rmtl.write(fname,format='fits', overwrite=True)



def ELGtilesi(tilef ):
	tiles = fitsio.read(tilef)
	trad = desimodel.focalplane.get_tile_radius_deg()*1.1 #make 10% greater just in case
	print(trad)
	rt = fitsio.read(minisvdir+'targets/MTL_all_SV0_ELG_tiles_0.37.0.fits')
	print('loaded random file')	
	
	for i in range(3,len(tiles)):
		tile = tiles['TILEID'][i]
		fname = minisvdir+'targets/MTL_TILE_ELG_'+str(tile)+'_0.37.0.fits'
		tdec = tiles['DEC'][i]
		decmin = tdec - trad
		decmax = tdec + trad
		wdec = (rt['DEC'] > decmin) & (rt['DEC'] < decmax)
		print(len(rt[wdec]))
		inds = desimodel.footprint.find_points_radec(tiles['RA'][i], tdec,rt[wdec]['RA'], rt[wdec]['DEC'])
		print('got indexes')
		fitsio.write(fname,rt[wdec][inds],clobber=True)
		print('wrote tile '+str(tile))


def targtilesi(type,tilef ):
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

def mktilef_date(dirout,fout='msvtiles.fits'):
	'''
	make a tile file for a date that Anand made tiles
	TBD
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


	
def mkminisvtilef(dirout,fout='msvtiles.fits'):
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

def mkminisvtilef_SV0(dirout,fout='msv0tiles.fits'):
	'''
	manually make tile fits file for minisv0 tiles
	'''
	msvtiles = Table()
	msvtiles['TILEID'] = np.array([68000,68001,68002,67142,67230],dtype=int)
	msvtiles['RA'] = np.array([214.75,214.76384,202.,204.136476102484,138.997356099811])
	msvtiles['DEC'] = np.array([53.4,53.408,8.25,5.90422737037591,0.574227370375913])
	msvtiles['PASS'] = np.zeros(5,dtype=int)
	msvtiles['IN_DESI'] = np.ones(5,dtype=int)
	msvtiles['OBSCONDITIONS'] = np.ones(5,dtype=int)*65535
	pa = []
	for i in range(0,5):
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