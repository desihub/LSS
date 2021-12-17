'''
python functions to do various useful date processing/manipulation
'''
import numpy as np
import fitsio
import glob
import os
import astropy.io.fits as fits
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt
import desimodel.footprint as foot
import desimodel.focalplane
from random import random
from desitarget.io import read_targets_in_tiles
from desitarget.targetmask import obsmask, obsconditions, zwarn_mask

#from LSS.Cosmo import distance
from LSS.imaging import densvar
from LSS.common_tools import find_znotposs
 


   
def combtile_spec(tiles,outf='',md=''):
    s = 0
    n = 0
    nfail = 0
    if os.path.isfile(outf):
        specd = Table.read(outf)
        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)
    else:
        tmask = np.ones(len(tiles)).astype('bool')    

    for tile,zdate in zip(tiles[tmask]['TILEID'],tiles[tmask]['ZDATE']):
        if md == 'zmtl':
            tspec = combzmtl(tile,zdate)
        else:
            tspec = combspecdata(tile,zdate)
        if tspec:
            tspec['TILEID'] = tile
            if s == 0:
                specd = tspec
                s = 1
            else:
                specd = vstack([specd,tspec],metadata_conflicts='silent')
            specd.sort('TARGETID')
            kp = (specd['TARGETID'] > 0)
            specd = specd[kp]
    
            n += 1
            print(tile,n,len(tiles[tmask]),len(specd)) 
        else:
            print(str(tile)+' failed')
            nfail += 1  
    print('total number of failures was '+str(nfail))
    specd.write(outf,format='fits', overwrite=True)       
 

def combspecdata(tile,zdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/archive/',md='' ):
    #put data from different spectrographs together, one table for fibermap, other for z
    zdate = str(zdate)
    specs = []
    #find out which spectrograph have data
    zfn = 'zbest'
    zhdu = 'ZBEST'
    shdu = 'SCORES'
    if int(zdate) >  20210730:
        zfn = 'redrock'
        zhdu = 'REDSHIFTS'
        #shdu = 'TSNR2' 
        

    for si in range(0,10):
        ff = coaddir+str(tile)+'/'+zdate+'/'+zfn+'-'+str(si)+'-'+str(tile)+'-thru'+zdate+'.fits'
        if os.path.isfile(ff):
            fq = coaddir+str(tile)+'/'+zdate+'/zmtl-'+str(si)+'-'+str(tile)+'-thru'+zdate+'.fits'
            if os.path.isfile(fq):

                specs.append(si)
            else:
                print('did not find '+fq)    
        elif zfn == 'zbest':
            zfnt = 'redrock'
            ff = coaddir+str(tile)+'/'+zdate+'/'+zfnt+'-'+str(si)+'-'+str(tile)+'-thru'+zdate+'.fits'
            if os.path.isfile(ff):
                fq = coaddir+str(tile)+'/'+zdate+'/zmtl-'+str(si)+'-'+str(tile)+'-thru'+zdate+'.fits'
                zfn = zfnt
                zhdu = 'REDSHIFTS'
                if os.path.isfile(fq):

                    specs.append(si)
                else:
                    print('did not find '+fq)    
            else:
                print('did not find '+ff)            
        else:
            print('did not find '+ff)        
    print('spectrographs with data on tile '+str(tile)+':')
    print(specs)            
    if len(specs) == 0:
        return None
    for i in range(0,len(specs)):
        tn = Table.read(coaddir+str(tile)+'/'+zdate+'/'+zfn+'-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu=zhdu)
        tnq = Table.read(coaddir+str(tile)+'/'+zdate+'/zmtl-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits')
        tnf = Table.read(coaddir+str(tile)+'/'+zdate+'/'+zfn+'-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='FIBERMAP')
        tns = Table.read(coaddir+str(tile)+'/'+zdate+'/coadd-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu=shdu)
    
        if i == 0:
           tspec = tn
           tq = tnq
           tf = tnf
           ts = tns
        else:    
            ts = vstack([ts,tns],metadata_conflicts='silent')
            tq = vstack([tq,tnq],metadata_conflicts='silent')
            tspec = vstack([tspec,tn],metadata_conflicts='silent')
            tf = vstack([tf,tnf],metadata_conflicts='silent')
        
    
    tf = unique(tf,keys=['TARGETID'])
    if md == '4combtar': #target files should contain the rest of the info
        tf.keep_columns(['FIBERASSIGN_X','FIBERASSIGN_Y','TARGETID','LOCATION','FIBERSTATUS','PRIORITY','DELTA_X','DELTA_Y','PSF_TO_FIBER_SPECFLUX','EXPTIME','OBJTYPE'])
    tq.keep_columns(['TARGETID','Z_QN','Z_QN_CONF','IS_QSO_QN','ZWARN'])
    tq['ZWARN'].name = 'ZWARN_MTL'
    tspec = join(tspec,tf,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')
    tspec = join(tspec,ts,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')
    tspec = join(tspec,tq,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')

    print(len(tspec),len(tf))
    #tspec['LOCATION'] = tf['LOCATION']
    #tspec['FIBERSTATUS'] = tf['FIBERSTATUS']
    #tspec['PRIORITY'] = tf['PRIORITY']
    return tspec

def combzmtl(tile,zdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/',md='' ):
    #put data from different spectrographs together, one table for fibermap, other for z
    zdate = str(zdate)
    specs = []
    #find out which spectrograph have data
    for si in range(0,10):
        ff = coaddir+str(tile)+'/'+zdate+'/zmtl-'+str(si)+'-'+str(tile)+'-thru'+zdate+'.fits'
        if os.path.isfile(ff):
            specs.append(si)
        else:
            print('did not find '+ff)    
    print('spectrographs with data:')
    print(specs)            
    if len(specs) == 0:
        return None
    for i in range(0,len(specs)):
        tn = Table.read(coaddir+str(tile)+'/'+zdate+'/zmtl-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits')
    
        if i == 0:
           tspec = tn
        else:    
            tspec = vstack([tspec,tn],metadata_conflicts='silent')
        
    
    tspec.keep_columns(['TARGETID','Z_QN','Z_QN_CONF','IS_QSO_QN','ZWARN','ZTILEID'])
    tspec['ZWARN'].name = 'ZWARN_MTL'
    tspec['ZTILEID'].name = 'TILEID'
    return tspec


def combfibmap(tile,zdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/' ):
    #put data from different spectrographs together, one table for fibermap, other for z
    specs = []
    #find out which spectrograph have data
    for si in range(0,10):
        
        #try:
        ff = coaddir+str(tile)+'/'+zdate+'/zbest-'+str(si)+'-'+str(tile)+'-thru'+zdate+'.fits'
        if os.path.isfile(ff):
            #fitsio.read(ff)
            specs.append(si)
        #except:
        #    print('no spectrograph '+str(si)+ ' for tile '+str(tile))
            #print(ff)
    #print('spectrographs with data:')
    #print(specs)            
    if len(specs) == 0:
        return None
    tf = Table.read(coaddir+str(tile)+'/'+zdate+'/zbest-'+str(specs[0])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='FIBERMAP')
    for i in range(1,len(specs)):
        tnf = Table.read(coaddir+str(tile)+'/'+zdate+'/zbest-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='FIBERMAP')
        tf = vstack([tf,tnf],metadata_conflicts='silent')
        
    
    tf = unique(tf,keys=['TARGETID'])
    tf.keep_columns(['FIBERASSIGN_X','FIBERASSIGN_Y','TARGETID','LOCATION','FIBERSTATUS','PRIORITY','DELTA_X','DELTA_Y','PSF_TO_FIBER_SPECFLUX','EXPTIME','OBJTYPE'])
    return tf

def combfibmap_and_scores(tile,zdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/' ):
    #put data from different spectrographs together, one table for fibermap, other for z
    specs = []
    #find out which spectrograph have data
    for si in range(0,10):
        
        #try:
        ff = coaddir+str(tile)+'/'+zdate+'/zbest-'+str(si)+'-'+str(tile)+'-thru'+zdate+'.fits'
        if os.path.isfile(ff):
            #fitsio.read(ff)
            specs.append(si)
        #except:
        #    print('no spectrograph '+str(si)+ ' for tile '+str(tile))
            #print(ff)
    #print('spectrographs with data:')
    #print(specs)            
    if len(specs) == 0:
        return None
    tf = Table.read(coaddir+str(tile)+'/'+zdate+'/zbest-'+str(specs[0])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='FIBERMAP')
    ts = Table.read(coaddir+str(tile)+'/'+zdate+'/coadd-'+str(specs[0])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='SCORES')
    for i in range(1,len(specs)):
        tnf = Table.read(coaddir+str(tile)+'/'+zdate+'/zbest-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='FIBERMAP')
        tf = vstack([tf,tnf],metadata_conflicts='silent')
        try:
            tns = Table.read(coaddir+str(tile)+'/'+zdate+'/coadd-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='SCORES')
            ts = vstack([ts,tns],metadata_conflicts='silent')
        except:
            print('did not find '+coaddir+str(tile)+'/'+zdate+'/coadd-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits')
    
    tf = unique(tf,keys=['TARGETID'])
    tf.keep_columns(['FIBERASSIGN_X','FIBERASSIGN_Y','TARGETID','LOCATION','FIBERSTATUS','PRIORITY','DELTA_X','DELTA_Y','PSF_TO_FIBER_SPECFLUX','EXPTIME','OBJTYPE'])
    tf = join(tf,ts,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')
    return tf


def goodlocdict(tf):
    '''
    Make a dictionary to map between location and priority
    tf should come from combspecdata above
    '''
    wloc = tf['FIBERSTATUS'] == 0
    print(str(len(tf[wloc])) + ' locations with FIBERSTATUS 0')
    goodloc = tf[wloc]['LOCATION']
    pdict = dict(zip(tf['LOCATION'], tf['PRIORITY'])) #to be used later for randoms
    return pdict,goodloc

def cutphotmask(aa,bits):
    print(str(len(aa)) +' before imaging veto' )
    keep = (aa['NOBS_G']>0) & (aa['NOBS_R']>0) & (aa['NOBS_Z']>0)
    for biti in bits:
        keep &= ((aa['MASKBITS'] & 2**biti)==0)
    aa = aa[keep]
    print(str(len(aa)) +' after imaging veto' )
    return aa

def combtiles_wdup(tiles,fout='',tarcol=['RA','DEC','TARGETID','DESI_TARGET','BGS_TARGET','MWS_TARGET','SUBPRIORITY','PRIORITY_INIT','TARGET_STATE','TIMESTAMP','ZWARN','PRIORITY']):
    s = 0
    n = 0
    if os.path.isfile(fout):
        tarsn = Table.read(fout)
        s = 1
        tdone = np.unique(tarsn['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)
    else:
        tmask = np.ones(len(tiles)).astype('bool')    
    for tile in tiles[tmask]['TILEID']:
        ts = str(tile).zfill(6)
        faf = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
        fht = fitsio.read_header(faf)
        mdir = '/global/cfs/cdirs/desi'+fht['MTL'][8:]+'/'
        if mdir == '/global/cfs/cdirs/desi/survey/ops/staging/mtl/main/dark/':
            mdir = '/global/cfs/cdirs/desi/target/catalogs/mtl/1.0.0/mtl/main/dark/'
        if mdir == '/global/cfs/cdirs/desi/survey/ops/staging/mtl/main/bright/':
            mdir = '/global/cfs/cdirs/desi/target/catalogs/mtl/1.0.0/mtl/main/bright/'
        wt = tiles['TILEID'] == tile
        tars = read_targets_in_tiles(mdir,tiles[wt],mtl=True,isodate=fht['MTLTIME'])
        #tars.keep_columns(tarcols)
        tars = tars[[b for b in tarcol]]
        
        tt = Table.read(faf,hdu='POTENTIAL_ASSIGNMENTS')
        tars = join(tars,tt,keys=['TARGETID'])
        tars['TILEID'] = tile
        tars.remove_columns(['ZWARN'])
        if s == 0:
            tarsn = tars
            s = 1
        else:
            tarsn = vstack([tarsn,tars],metadata_conflicts='silent')
        tarsn.sort('TARGETID')
        n += 1
        print(tile,n,len(tiles[tmask]),len(tarsn)) 
    tarsn.write(fout,format='fits', overwrite=True)       

def combtiles_wdup_hp(hpx,tiles,fout='',tarcol=['RA','DEC','TARGETID','DESI_TARGET','BGS_TARGET','MWS_TARGET','SUBPRIORITY','PRIORITY_INIT','TARGET_STATE','TIMESTAMP','ZWARN','PRIORITY']):
    s = 0
    n = 0
    
    tls = foot.pix2tiles(8,[hpx],tiles)
    if os.path.isfile(fout):
        tarsn = Table.read(fout)
        s = 1
        tdone = np.unique(tarsn['TILEID'])
        tmask = ~np.isin(tls['TILEID'],tdone)
    else:
        tmask = np.ones(len(tls)).astype('bool')    
    for tile in tls[tmask]['TILEID']:
        ts = str(tile).zfill(6)
        faf = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
        fht = fitsio.read_header(faf)
        mdir = '/global/cfs/cdirs/desi'+fht['MTL'][8:]+'/'
        if mdir == '/global/cfs/cdirs/desi/survey/ops/staging/mtl/main/dark/':
            mdir = '/global/cfs/cdirs/desi/target/catalogs/mtl/1.0.0/mtl/main/dark/'
        if mdir == '/global/cfs/cdirs/desi/survey/ops/staging/mtl/main/bright/':
            mdir = '/global/cfs/cdirs/desi/target/catalogs/mtl/1.0.0/mtl/main/bright/'
        wt = tls['TILEID'] == tile
        tars = read_targets_in_tiles(mdir,tiles[wt],mtl=True,isodate=fht['MTLTIME'])
        #tars.keep_columns(tarcols)
        tars = tars[[b for b in tarcol]]
        theta, phi = np.radians(90-tars['DEC']), np.radians(tars['RA'])
        tpix = hp.ang2pix(8,theta,phi,nest=True)
        sel = tpix == hpx
        tars = tars[sel]
        tt = Table.read(faf,hdu='POTENTIAL_ASSIGNMENTS')
        tars = join(tars,tt,keys=['TARGETID'])
        tars['TILEID'] = tile
        tars.remove_columns(['ZWARN'])
        if s == 0:
            tarsn = tars
            s = 1
        else:
            tarsn = vstack([tarsn,tars],metadata_conflicts='silent')
        tarsn.sort('TARGETID')
        n += 1
        print(tile,n,len(tls[tmask]),len(tarsn)) 
    tarsn.write(fout,format='fits', overwrite=True)       


def gettarinfo_type(faf,tars,goodloc,pdict,tp='SV3_DESI_TARGET'):
    #get target info
    #in current files on SVN, TARGETS has all of the necessary info on potential assignments
    #no more, so commented out
    #tt = Table.read(faf,hdu='TARGETS')
    #tt.keep_columns(['TARGETID','FA_TARGET','FA_TYPE','PRIORITY','SUBPRIORITY','OBSCONDITIONS'])
    tt = Table.read(faf,hdu='POTENTIAL_ASSIGNMENTS')
    #if len(tt) != len(tfa):
    #    print('!!!mismatch between targets and potential assignments, aborting!!!')
    #    return None
    #tt = join(tt,tfa,keys=['TARGETID'])    

    
    
    wgt = (np.isin(tt['LOCATION'],goodloc)) 
    print(str(len(np.unique(tt[wgt]['LOCATION']))) + ' good locations')
    print('comparison of number targets, number of targets with good locations')
    print(len(tt),len(tt[wgt]))
    
    tt = tt[wgt]
    
    tt = join(tt,tars,keys=['TARGETID'],table_names = ['_AVAIL', ''], uniq_col_name='{col_name}{table_name}')
    
    #Mark targets that actually got assigned fibers
    tfall = Table.read(faf,hdu='FIBERASSIGN')
    
    tfall.keep_columns(['TARGETID','LOCATION','PRIORITY'])
    
    tt = join(tt,tfall,keys=['TARGETID'],join_type='left',table_names = ['', '_ASSIGNED'], uniq_col_name='{col_name}{table_name}')
    wal = tt['LOCATION_ASSIGNED']*0 == 0
    tt['LOCATION'][wal] = tt['LOCATION_ASSIGNED'][wal]
    tt['LOCATION_AVAIL'][wal] = tt['LOCATION_ASSIGNED'][wal]
    #print('differences between assigned locations')
    #print(np.unique(tt['LOCATION_AVAIL'][wal]-tt['LOCATION_ASSIGNED'][wal]))
    #print(tt.columns)


    tt = unique(tt,keys=['TARGETID']) #cut to unique target ids
    #print(tarf)
    #tars = Table.read(tarf)
    #tars.remove_columns(['Z','ZWARN'])#,'PRIORITY','SUBPRIORITY','OBSCONDITIONS'])
    #we want to get these from the zbest file that is specific to the tile and thus when it was observed
    
    
    
    #tfa = unique(tfa[wgt],keys=['TARGETID'])
    #wtype = ((tt[tp] & 2**tarbit) > 0) #don't cut by type here any more
    #tt = tt[wtype]
    
    #tfa = join(tfa,tt,keys=['TARGETID'])
    #tft = join(tft,tt,keys=['TARGETID'])
    #print(str(len(tfa)) +' unique targets with good locations and  at '+str(len(np.unique(tfa['LOCATION'])))+' unique locations and '+str(len(tft))+ ' total unique targets at '+str(len(np.unique(tft['LOCATION']))) +' unique locations ')

    #wgl = np.isin(tfa['LOCATION_ASSIGNED'],goodloc)
    #wtype = ((tfa[tp] & 2**tarbit) > 0)
    #wtfa = wgl & wtype
    #print('number of assigned fibers at good locations '+str(len(tfa[wtfa])))

    wal = tt['LOCATION_ASSIGNED']*0 == 0
    print('number of assigned fibers '+str(len(tt[wal])))
    print('number of unique target id '+str(len(np.unique(tt[wal]['TARGETID']))))
    print('max priority of assigned '+str(np.max(tt[wal]['PRIORITY_ASSIGNED'])))
    #tt[wal]['LOCATION'] = tt[wal]['LOCATION_ASSIGNED']
    #tt[wal]['LOCATION_AVAIL'] = tt[wal]['LOCATION_ASSIGNED']
    #print('are location and location_avail the same for assigned targets?')
    #print(np.array_equal(tt[wal]['LOCATION'], tt[wal]['LOCATION_AVAIL']))
    #print('are location_avail and location_assigned the same for assigned targets?')
    #print(np.array_equal(tt[wal]['LOCATION_ASSIGNED'], tt[wal]['LOCATION_AVAIL']))

    tt['LOCATION_ASSIGNED'] = np.zeros(len(tt),dtype=int)
    tt['LOCATION_ASSIGNED'][wal] = 1
    wal = tt['LOCATION_ASSIGNED'] == 1
    print('number of assigned fibers '+str(len(tt[wal]))+' (check to match agrees with above)')
    wal = tt['LOCATION']*0 == 0
    print('number of locations from z file '+str(len(tt[wal]))+' (check to match agrees with above)')
    #print('are location and location_avail the same for assigned targets?')
    #print(np.array_equal(tt[wal]['LOCATION'], tt[wal]['LOCATION_AVAIL']))
    #tt['PRIORITY_ASSIGNED'] = np.vectorize(pdict.__getitem__)(tt['LOCATION'])

    return tt

    
def get_specdat(indir,pd):
    #indir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+specrel
    zf = indir+'/datcomb_'+pd+'_tarspecwdup_zdone.fits'
    dz = Table.read(zf) 
    selz = dz['ZWARN'] != 999999
    fs = dz[selz]

    #first, need to find locations to veto based data
    nodata = fs["ZWARN_MTL"] & zwarn_mask["NODATA"] != 0
    num_nod = np.sum(nodata)
    print('number with no data '+str(num_nod))
    badqa = fs["ZWARN_MTL"] & zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA") != 0
    num_badqa = np.sum(badqa)
    print('number with bad qa '+str(num_badqa))
    nomtl = nodata | badqa
    wfqa = ~nomtl
    return fs[wfqa]


def count_tiles_better(dr,pd,rann=0,specrel='daily',fibcol='COADD_FIBERSTATUS'):
    '''
    from files with duplicates that have already been sorted by targetid, quickly go 
    through and get the multi-tile information
    dr is either 'dat' or 'ran'
    returns file with TARGETID,NTILE,TILES,TILELOCIDS
    '''
    
    #fs = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+specrel+'/datcomb_'+pd+'_spec_zdone.fits')
    #wf = fs['FIBERSTATUS'] == 0
    #wf = fs[fibcol] == 0
    #nodata = fs["ZWARN_MTL"] & zwarn_mask["NODATA"] != 0
    #num_nod = np.sum(nodata)
    #print('number with no data '+str(num_nod))
    #badqa = fs["ZWARN_MTL"] & zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA") != 0
    #num_badqa = np.sum(badqa)
    #print('number with bad qa '+str(num_badqa))
    #nomtl = nodata & badqa
    #wfqa = ~nomtl
    
    indir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+specrel
    fs = get_specdat(indir,pd)

    stlid = 10000*fs['TILEID'] +fs['LOCATION']
    gtl = np.unique(stlid)
    
    if dr == 'dat':
        fj = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+specrel+'/datcomb_'+pd+'_tarspecwdup_zdone.fits')
        #outf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/datcomb_'+pd+'ntileinfo.fits' 
    if dr == 'ran':
        fj = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+specrel+'/rancomb_'+str(rann)+pd+'wdupspec_zdone.fits')
        
        #outf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random'+str(rann)+'/rancomb_'+pd+'ntileinfo.fits'
    wg = np.isin(fj['TILELOCID'],gtl)  
    fjg = fj[wg]  

    tids = np.unique(fjg['TARGETID'])
    nloc = []#np.zeros(len(np.unique(f['TARGETID'])))
    nt = []
    tl = []
    tli = []
    ti = 0
    i = 0
    while i < len(fjg):
        tls  = []
        tlis = []
        nli = 0
    
        while fjg[i]['TARGETID'] == tids[ti]:
            nli += 1
            tls.append(fjg[i]['TILEID'])
            tlis.append(fjg[i]['TILELOCID'])
            i += 1
            if i == len(fjg):
                break
        nloc.append(nli)
        tlsu = np.unique(tls)
        tlisu = np.unique(tlis)
        nt.append(len(tlsu))
        tl.append("-".join(tlsu.astype(str)))
        tli.append("-".join(tlisu.astype(str)))
      
        if ti%100000 == 0:
            print(ti)
        ti += 1 
    tc = Table()
    tc['TARGETID'] = tids
    tc['NTILE'] = nt
    tc['TILES'] = tl
    tc['TILELOCIDS'] = tli
    
    return tc
    
    
def count_tiles(tiles,catdir,pd,ttp='ALL',imask=False):
    '''
    For list of tileids, simply track the tiles a target shows up as available in
    pd is dark or bright
    just output targetid and tiles, meant to be matched to other processing
    don't worry about what was assigned, purpose is to just count tile overlaps
    '''

    s = 0
    cnt = 0 
    for tile in tiles:
        fl = catdir+ttp+str(tile)+'_full.dat.fits'
        fgun = Table.read(fl)
        if imask:
            wm = fgun['MASKBITS'] == 0
            fgun = fgun[wm]
        fgun['TILELOCID'] = 10000*tile +fgun['LOCATION_AVAIL']
        fgun.keep_columns(['TARGETID','TILELOCID'])
        print(len(fgun),len(np.unique(fgun['TARGETID'])))

        aa = np.chararray(len(fgun),unicode=True,itemsize=100)
        aa[:] = str(tile)

        fgun['TILES'] = aa
        
        ai = np.chararray(len(fgun),unicode=True,itemsize=300)
        tlids = np.copy(fgun['TILELOCID']).astype('<U300')
        fgun['TILELOCIDS'] = tlids

        if s == 0:
            fgu = fgun
            s =1
        else:
            fgo = fgu.copy()
            fgu = vstack([fgu,fgun],metadata_conflicts='silent')
            fgu = unique(fgu,keys='TARGETID')#,keep='last') 
                
            #I think this works when the ordering is the same; things got messed up other places with sorts
            dids = np.isin(fgun['TARGETID'],fgo['TARGETID']) #get the rows with target IDs that were duplicates in the new file
            didsc = np.isin(fgu['TARGETID'],fgun['TARGETID'][dids]) #get the row in the concatenated table that had dup IDs

            aa = np.chararray(len(fgu['TILES']),unicode=True,itemsize=20)
            aa[:] = '-'+str(tile)
            #rint(aa)
            ms = np.core.defchararray.add(fgu['TILES'][didsc],aa[didsc])
            #print(ms)
            fgu['TILES'][didsc] = ms #add the tile info
            aa = np.copy(fgun[dids]['TILELOCIDS'])#np.chararray(len(fgu['TILELOCIDS']),unicode=True,itemsize=100)
            aa[:] = np.core.defchararray.add('-',aa)
            
            #rint(aa)
            ms = np.core.defchararray.add(fgu['TILELOCIDS'][didsc],aa)
            #print(ms)
            fgu['TILELOCIDS'][didsc] = ms #add the tile info


        print(tile,cnt,len(tiles),len(fgu))
        cnt += 1

    fu = fgu
    fl = np.chararray(len(fu),unicode=True,itemsize=100)
    for ii in range(0,len(fu)):
        tl = fu['TILES'][ii]
        tls = tl.split('-')#np.unique()#.astype('int')
        tli = tls[0]
        if len(tls) > 1:
            #tls = tls.astype('int')
            tls.sort()
            tli = tls[0]
            for i in range(1,len(tls)):
                tli += '-'+tls[i]
        #else:
        #    tli = tls
        #print(tli)
        fl[ii] = tli   
    
    fu['TILES'] = fl
    print(np.unique(fu['TILES']))
    
    fu.write(catdir+'Alltiles_'+pd+'_tilelocs.dat.fits',format='fits', overwrite=True)    


def combtiles(tiles,catdir,tp,tmask,tc='SV3_DESI_TARGET',ttp='ALL',imask=False):
    '''
    For list of tileids, combine data generated per tile , taking care of overlaps
    
    '''

    s = 0
    cnt = 0 
    for tile in tiles:
        fl = catdir+ttp+str(tile)+'_full.dat.fits'
        fgun = Table.read(fl)
        if imask:
            wm = fgun['MASKBITS'] == 0
            fgun = fgun[wm]

        if tp != 'dark' and tp != 'bright':
            wt = (fgun[tc] & tmask[tp]) > 0
            fgun = fgun[wt]
        fgun['TILELOCID'] = 10000*tile +fgun['LOCATION_AVAIL']
        fgun['TILELOCID_ASSIGNED'] = np.zeros(len(fgun))
        wm = fgun['LOCATION_ASSIGNED'] == 1
        fgun['TILELOCID_ASSIGNED'][wm] = fgun['TILELOCID'][wm]
        nl,nla = countloc(fgun)
        fgun['ZPOSS'] = np.zeros(len(fgun)).astype(int)
        if tp != 'dark' and tp != 'bright':
            #fgun['LOC_NOTBLOCK'] = np.zeros(len(fgun)).astype(int)
            locsna = []
            for i in range(0,len(nla)):
                if nla[i] == 0 and nl[i] > 0:
                    locsna.append(i)

            print('number of unassigned locations',len(locsna))
            was = ~np.isin(fgun['LOCATION_AVAIL'],locsna)
            #fgun['LOC_NOTBLOCK'][was] = 1
            wg = was
            fgun['ZPOSS'][wg] = 1
            #fgun.sort('ZPOSS')

        #aa = np.chararray(len(fgun),unicode=True,itemsize=100)
        #aa[:] = str(tile)
        fgun['TILE'] = int(tile)
        #fgun['TILES'] = aa
        #tlids = np.copy(fgun['TILELOCID']).astype('<U300')
        #fgun['TILELOCIDS'] = tlids
        
        #print('sum of assigned,# of unique TILELOCID (should match)')
        #print(np.sum(fgun['LOCATION_ASSIGNED'] == 1),len(np.unique(fgun['TILELOCID'])))
        #ai = np.chararray(len(fgun),unicode=True,itemsize=300)
        #
        #

        if s == 0:
            fgu = fgun
            s =1
        else:
            #fgo = fgu.copy()
            fgu = vstack([fgu,fgun],metadata_conflicts='silent')
            #wn = fgu['PRIORITY_ASSIGNED']*0 != 0
            #wn |= fgu['PRIORITY_ASSIGNED'] == 999999
            #print(len(fgu[~wn]),np.max(fgu[~wn]['PRIORITY_ASSIGNED']),'max priority assigned')
            #fgu[wn]['PRIORITY_ASSIGNED'] = 0
            #fgu['sort'] = -1.*fgu['LOCATION_ASSIGNED']*fgu['PRIORITY_ASSIGNED'] #create this column so assigned always show up in order of highest priority
            #wa = fgu['LOCATION_ASSIGNED'] == 1
            #wa &= fgu['PRIORITY_ASSIGNED'] >= 2000 #this was put SV2 to ignore BGS repeats
            #fa = fgu[wa]
            #print(len(fa),len(np.unique(fa['TARGETID'])))
            #fgu.sort('sort')
            #fgu = unique(fgu,keys='TARGETID',keep='last') 
                
            #dids = np.isin(fgun['TARGETID'],fgo['TARGETID']) #get the rows with target IDs that were duplicates in the new file
            #didsc = np.isin(fgu['TARGETID'],fgun['TARGETID'][dids]) #get the row in the concatenated table that had dup IDs
            #print(len(fgu),len(fgo),len(fgun),len(fgu[didsc]),len(fgun[dids]))
            #fgu['TILELOCID'][didsc] = fgun['TILELOCID'][dids] #give the repeats the new tilelocids, since those are the most likely to be available to low priority targets
            #if tp != 'dark' and tp != 'bright':
            #    fgu['LOC_NOTBLOCK'][didsc] = np.maximum(fgu['LOC_NOTBLOCK'][didsc],fgun['LOC_NOTBLOCK'][dids]) 
            #    fgu['ZPOSS'][didsc] = np.maximum(fgu['ZPOSS'][didsc],fgun['ZPOSS'][dids]) 

            #aa = np.chararray(len(fgu['TILES']),unicode=True,itemsize=20)
            #aa[:] = '-'+str(tile)
            #rint(aa)
            #ms = np.core.defchararray.add(fgu['TILES'][didsc],aa[didsc])
            #print(ms)
            #fgu['TILES'][didsc] = ms #add the tile info
            #aa = np.copy(fgun[dids]['TILELOCIDS'])#np.chararray(len(fgu['TILELOCIDS']),unicode=True,itemsize=100)
            #aa[:] = np.core.defchararray.add('-',aa)
            
            #rint(aa)
            #ms = np.core.defchararray.add(fgu['TILELOCIDS'][didsc],aa)
            #print(ms)
            #fgu['TILELOCIDS'][didsc] = ms #add the tile info


        print(tile,cnt,len(tiles))#,np.sum(fgu['LOCATION_ASSIGNED']),len(fgu),len(np.unique(fgu['TILELOCID'])),np.sum(fgu['ZPOSS']))#,np.unique(fgu['TILELOCIDS'])
        cnt += 1

    #fgu['TILES'] = np.copy(fgu['TILE']).astype('<U100')
    #tlids = np.copy(fgu['TILELOCID']).astype('<U300')
    #fgu['TILELOCIDS'] = tlids
    
    tsnrcol = 'TSNR2_'+tp
    if tp == 'ELG_HIP':
        tsnrcol = 'TSNR2_ELG'
    if tp == 'BGS_ANY':
        tsnrcol = 'TSNR2_BGS'    
    wt = (fgu[tsnrcol] == 1e20) | (fgu[tsnrcol]*0 != 0)
    print('number with bad tsnrcol is '+str(len(fgu[wt])))
    fgu[tsnrcol][wt] = 0
    wn = fgu['PRIORITY_ASSIGNED']*0 != 0
    wn |= fgu['PRIORITY_ASSIGNED'] == 999999
    #print(len(fgu[~wn]),np.max(fgu[~wn]['PRIORITY_ASSIGNED']),'max priority assigned')
    fgu[wn]['PRIORITY_ASSIGNED'] = 0
    fgu['sort'] = -1.*fgu['LOCATION_ASSIGNED']*fgu['PRIORITY_ASSIGNED']*fgu[tsnrcol] #create this column so assigned always show up in order of highest priority
   
    
    if tp != 'dark' and tp != 'bright':
        #wa = fgu['LOCATION_ASSIGNED'] == 1
        #print('ZPOSS for LOCATION_ASSIGNED = 1:')
        #print(np.unique(fgu[wa]['ZPOSS']))
        fgu['sort'] = fgu['sort']*fgu['ZPOSS']-fgu['ZPOSS']
        wa = fgu['LOCATION_ASSIGNED'] == 1
        #wp = fgu['ZPOSS']
        loclz,nloclz = np.unique(fgu[wa]['TILELOCID_ASSIGNED'],return_counts=True)
        wp = fgu['ZPOSS'] == 1
        natloc = ~np.isin(fgu[wp]['TILELOCID'],loclz)
        print('number of zposs with tilelocid not showing up in tilelocid_assigned:')
        print(np.sum(natloc))
    fgu.sort('sort')
    #fgu.sort('ZPOSS')
    fu = unique(fgu,keys='TARGETID')#,keep='last')
 
    tidsu = fu['TARGETID']#[wp][natloc]
    tids = fgu['TARGETID']
 
    
    if tp != 'dark' and tp != 'bright':
        
        wa = fu['LOCATION_ASSIGNED'] == 1
        #wp = fgu['ZPOSS']
        loclz,nloclz = np.unique(fu[wa]['TILELOCID_ASSIGNED'],return_counts=True)
        wp = fu['ZPOSS'] == 1
        nalz = ~np.isin(fu['TILELOCID'],loclz)
        natloc = wp & nalz#~np.isin(fu[wp]['TILELOCID'],loclz)
        print('after cutting to unique, number of zposs with tilelocid not showing up in tilelocid_assigned:')
        print(np.sum(natloc))
        tlocs = fgu['TILELOCID']
        ntl = []
        ch = 0
        bl = 0
        print(len(tidsu),len(natloc))
        for ii in range(0,len(tidsu)):
            #if wp[ii] & natloc[ii]:
            if natloc[ii]:
                bl += 1
                tid = tidsu[ii]
                wt = tids == tid
                tls = tlocs[wt]
                s = 0
                for tl in tls:
                    if s == 0:
                        if np.isin(tl,loclz):
                            #wu = fu['TARGETID'] == tid
                            fu[ii]['TILELOCID'] = tl
                            #ntl.append(tl)
                            ch += 1
                            s = 1
            if ii%10000 == 0:
                print(ii,len(tidsu),ch,bl)
        wa = fu['LOCATION_ASSIGNED'] == 1
        #wp = fgu['ZPOSS']
        loclz,nloclz = np.unique(fu[wa]['TILELOCID_ASSIGNED'],return_counts=True)
        wp = fu['ZPOSS'] == 1
        natloc = ~np.isin(fu[wp]['TILELOCID'],loclz)
        print('after cutting to unique and reassignment, number of zposs with tilelocid not showing up in tilelocid_assigned:')
        print(np.sum(natloc))
        
        
    #print(len(np.unique(fgu['TARGETID'])),np.sum(fgu['LOCATION_ASSIGNED']))
    
#     tiles = fgu['TILES']
#     tilesu = fu['TILES']
#     tlids = fgu['TILELOCIDS']
#     tlidsu = fu['TILELOCIDS']
# 
#     for ii in range(0,len(tidsu)): #this takes a long time and something more efficient will be necessary
#         tid = tidsu[ii]#fu[ii]['TARGETID']
#         wt = tids == tid
#         ot = tilesu[ii]
#         otl = tlidsu[ii]
#         tt = tiles[wt]
#         tti = tlids[wt]
#         for tl in tt:
#             if tl != ot:
#                 tilesu[ii] += '-'+str(tl)
#         for ti in tti:
#             if ti != otl:
#                 tlidsu[ii] += '-'+str(ti)
#         if ii%1000 == 0:
#             print(ii)        
#     fu['TILES'] = tilesu
#     fu['TILELOCIDS'] = tlidsu
#     
#     #wa = fu['LOCATION_ASSIGNED'] == 1
#     #wa &= fu['PRIORITY_ASSIGNED'] >= 2000
    print(np.sum(fu['LOCATION_ASSIGNED']))

    #need to resort tile string
#     fl = np.chararray(len(fu),unicode=True,itemsize=100)
#     for ii in range(0,len(fu)):
#         tl = fu['TILES'][ii]
#         tls = tl.split('-')#.astype('int')
#         tli = tls[0]
#         if len(tls) > 1:
#             #tls = tls.astype('int')
#             tls.sort()
#             tli = tls[0]
#             for i in range(1,len(tls)):
#                 tli += '-'+tls[i]
#         #else:
#         #    tli = tls
#         #print(tli)
#         fl[ii] = tli   
#     
#     fu['TILES'] = fl
    #print(np.unique(fu['TILES']))
#     print('number of unique tiles configurations '+str(len(np.unique(fu['TILES']))))
    #fu.write(catdir+tp+'Alltiles_'+pd+'_full.dat.fits',format='fits', overwrite=True)    
    fu.write(catdir+'/datcomb_'+tp+'_Alltiles.fits',format='fits', overwrite=True)

def countloc(aa):
    locs = aa['LOCATION_AVAIL']
    locsa = aa['LOCATION_ASSIGNED']
    la = np.max(locs)+1
    nl = np.zeros(la)
    nla = np.zeros(la)
    for i in range(0,len(aa)):
        nl[locs[i]] += 1
        nla[locs[i]] += locsa[i]
    return nl,nla


def combran_wdup(tiles,rann,randir,tp,lspecdir,specf,keepcols=[]):

    s = 0
    td = 0
    #tiles.sort('ZDATE')
    print(len(tiles))
    delcols = ['DESI_TARGET','BGS_TARGET','MWS_TARGET','SUBPRIORITY','OBSCONDITIONS','PRIORITY_INIT',\
    'NUMOBS_INIT','SCND_TARGET','NUMOBS_MORE','NUMOBS','Z','ZWARN','TARGET_STATE','TIMESTAMP','VERSION','PRIORITY']
    outf = randir+str(rann)+'/rancomb_'+tp+'wdup_Alltiles.fits'

    if os.path.isfile(outf):
        fgu = Table.read(outf)
        #tarsn.keep_columns(['RA','DEC','TARGETID''LOCATION','FIBER','TILEID'])
        s = 1
        tdone = np.unique(fgu['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)
    else:
        tmask = np.ones(len(tiles)).astype('bool')    
    for tile in tiles[tmask]['TILEID']:
        ffa = randir+str(rann)+'/fba-'+str(tile).zfill(6)+'.fits'
        ffna = randir+str(rann)+'/tilenofa-'+str(tile)+'.fits'
        if os.path.isfile(ffa):
            fa = Table.read(ffa,hdu='FAVAIL')
            
            ffna = Table.read(ffna)
            fgun = join(fa,ffna,keys=['TARGETID'])
            #fgun.remove_columns(delcols)
                        
            td += 1
            fgun['TILEID'] = int(tile)
            fgun.keep_columns(['RA','DEC','TARGETID','LOCATION','FIBER','TILEID'])
            if s == 0:
                fgu = fgun
                s = 1
            else:   
                fgu = vstack([fgu,fgun],metadata_conflicts='silent')
            fgu.sort('TARGETID')
            print(tile,td, len(tiles), len(fgun),len(fgu))
        else:
            print('did not find '+ffa)

    if len(tiles[tmask]['TILEID']) > 0:
        fgu.write(outf,format='fits', overwrite=True)
    #specf = Table.read(lspecdir+'datcomb_'+tp+'_spec_zdone.fits')
    specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    specf.keep_columns(keepcols)
    #specf.keep_columns(['ZWARN','LOCATION','TILEID','TILELOCID','FIBERSTATUS','FIBERASSIGN_X','FIBERASSIGN_Y','PRIORITY','DELTA_X','DELTA_Y','EXPTIME','PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B','TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z','TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
    fgu = join(fgu,specf,keys=['LOCATION','TILEID','FIBER'],join_type='left')
    fgu.sort('TARGETID')
    outf = lspecdir+'/rancomb_'+str(rann)+tp+'wdupspec_zdone.fits'
    print(outf)
    fgu.write(outf,format='fits', overwrite=True)
    


def combran(tiles,rann,randir,ddir,tp,tmask,tc='SV3_DESI_TARGET',imask=False):

    s = 0
    td = 0
    #tiles.sort('ZDATE')
    print(len(tiles))
    delcols = ['DESI_TARGET','BGS_TARGET','MWS_TARGET','SUBPRIORITY','OBSCONDITIONS','PRIORITY_INIT',\
    'NUMOBS_INIT','SCND_TARGET','NUMOBS_MORE','NUMOBS','Z','ZWARN','TARGET_STATE','TIMESTAMP','VERSION','PRIORITY']
    for tile,zdate in zip(tiles['TILEID'],tiles['ZDATE']):
        tspec = combfibmap_and_scores(tile,zdate)
        pdict,gloc = goodlocdict(tspec)
        tspec.keep_columns(['LOCATION','FIBERSTATUS','DELTA_X','DELTA_Y','PSF_TO_FIBER_SPECFLUX','EXPTIME','OBJTYPE','TSNR2_ELG','TSNR2_LRG','TSNR2_QSO','TSNR2_BGS'])
        dt = ddir+'ALL'+str(tile)+'_full.dat.fits'
        ffa = randir+str(rann)+'/fba-'+str(tile).zfill(6)+'.fits'
        ffna = randir+str(rann)+'/tilenofa-'+str(tile)+'.fits'
        if os.path.isfile(ffa):
            fd = Table.read(dt)
           # print(np.sum(fd['LOCATION_ASSIGNED']),len(fd))
            #gloc = np.unique(fd['LOCATION_AVAIL']) #bad locations already removed from this files
            #print(np.sum(fd['LOCATION_ASSIGNED']),len(fd),len(gloc))
            if tp != 'dark' and tp != 'bright':
                wt = (fd[tc] & tmask[tp]) > 0
                fd = fd[wt]
            #print(np.sum(fd['LOCATION_ASSIGNED']),len(fd))
            nl,nla = countloc(fd)
            #commenting out zfailure stuff, not vetoing randoms based on that
            #wzf = fd['ZWARN'] != 0 
            #wzf &= fd['ZWARN'] != 999999
            #wzf &= fd['ZWARN']*0 == 0
            #loc_fail = np.unique(fd[wzf]['LOCATION'])
            #print('number of zfail locations',len(loc_fail))
            #
            #print(np.sum(fd['LOCATION_ASSIGNED']),len(np.unique(fd['LOCATION_AVAIL'])),np.sum(nla),np.sum(nl))
        # 
            #find the locations that were requested by type but not assigned
            fa = Table.read(ffa,hdu='FAVAIL')
            
            wg = np.isin(fa['LOCATION'],gloc)
            fa = fa[wg]
            fa = join(fa,tspec,keys=['LOCATION'],join_type='left')
            #fa['FIBER_GOOD'] = np.zeros(len(fa)).astype(int)
            #fa['FIBER_GOOD'][wg] = 1
            #fa['Z_NOTBAD'] = np.zeros(len(fa)).astype(int)
            #wnzf = ~np.isin(fa['LOCATION'],loc_fail)
            #fa['Z_NOTBAD'][wnzf] = 1
            fa['ZPOSS'] = np.zeros(len(fa)).astype(int)
            #fa['ZPOSSNOTBAD'] = np.zeros(len(fa)).astype(int)
            if tp != 'dark' and tp != 'bright':
                #fa['LOC_NOTBLOCK'] = np.zeros(len(fa)).astype(int)
                locsna = []
                for i in range(0,len(nla)):
                    if nla[i] == 0 and nl[i] > 0:
                        locsna.append(i)

                print('number of unassigned locations',len(locsna))
                ntloc = len(gloc)-len(locsna)#-len(loc_fail)
                print('total number of assignable positions',ntloc)
                was = ~np.isin(fa['LOCATION'],locsna)
                #fa['LOC_NOTBLOCK'][was] = 1
                #wg &= was
                fa['ZPOSS'][was] = 1
                #fa['ZPOSSNOTBAD'][was&wnzf] = 1
                #if maskzfail:
                #    wg &= wnzf
            
            #wzt = wpr & ~wzf & ~wna

            #fg = fa[wg]
            #print(len(fa),np.sum(fa['ZPOSSNOTBAD']))
            #fg = fa
            #print('before,after vetoing locations:')
            #print(len(fa),len(fg))
            #if tp != 'dark' and tp != 'bright':
            #    fa.sort('ZPOSS')
            #else:
            #    fg.sort('FIBER_GOOD') 
            fgun = unique(fa,keys=['TARGETID'],keep='last')
            ffna = Table.read(ffna)
            fgun = join(fgun,ffna,keys=['TARGETID'])
            fgun.remove_columns(delcols)
            if imask:
                wm = fgun['MASKBITS'] == 0
                fgun = fgun[wm]

            print(tile,td, len(tiles), str(len(fgun))+' unique new randoms')
            td += 1
            aa = np.chararray(len(fgun),unicode=True,itemsize=100)
            aa[:] = str(tile)
            fgun['TILE'] = int(tile)
            fgun['TILES'] = aa
            fgun['TILELOCID'] = 10000*tile +fgun['LOCATION']
            if s == 0:
                fgu = fgun
                s = 1
            else:   
                fv = vstack([fgu,fgun],metadata_conflicts='silent')
                fgo = fgu.copy()
                fgu = unique(fv,keys='TARGETID')#,keep='last') 
                
                dids = np.isin(fgun['TARGETID'],fgo['TARGETID']) #get the rows with target IDs that were duplicates in the new file
                didsc = np.isin(fgu['TARGETID'],fgun['TARGETID'][dids]) #get the row in the concatenated table that had dup IDs
                #print(len(fgu),len(fgo),len(fgun),len(fgu[didsc]),len(fgun[dids]))
                fgu['TILELOCID'][didsc] = fgun['TILELOCID'][dids] #give the repeats the new tilelocids, since those are the most likely to be available to low priority targets
                #if this works, can save vetoing until the end
                fgu['TSNR2_ELG'][didsc] = np.maximum(fgu['TSNR2_ELG'][didsc],fgun['TSNR2_ELG'][dids])
                fgu['TSNR2_QSO'][didsc] = np.maximum(fgu['TSNR2_QSO'][didsc],fgun['TSNR2_QSO'][dids]) 
                fgu['TSNR2_BGS'][didsc] = np.maximum(fgu['TSNR2_BGS'][didsc],fgun['TSNR2_BGS'][dids])
                fgu['TSNR2_LRG'][didsc] = np.maximum(fgu['TSNR2_LRG'][didsc],fgun['TSNR2_LRG'][dids]) 
                if tp != 'dark' and tp != 'bright':
                     #fgu['FIBER_GOOD'][didsc] = np.maximum(fgu['FIBER_GOOD'][didsc],fgun['FIBER_GOOD'][dids])
                     #fgu['LOC_NOTBLOCK'][didsc] = np.maximum(fgu['LOC_NOTBLOCK'][didsc],fgun['LOC_NOTBLOCK'][dids]) 
                     #fgu['Z_NOTBAD'][didsc] = np.maximum(fgu['Z_NOTBAD'][didsc],fgun['Z_NOTBAD'][dids])
                     fgu['ZPOSS'][didsc] = np.maximum(fgu['ZPOSS'][didsc],fgun['ZPOSS'][dids]) 
                     #fgu['ZPOSSNOTBAD'][didsc] = np.maximum(fgu['ZPOSSNOTBAD'][didsc],fgun['ZPOSSNOTBAD'][dids])

                aa = np.chararray(len(fgu['TILES']),unicode=True,itemsize=20)
                aa[:] = '-'+str(tile)
                #rint(aa)
                ms = np.core.defchararray.add(fgu['TILES'][didsc],aa[didsc])
                #print(ms)
                fgu['TILES'][didsc] = ms #add the tile info
                print(str(len(fgu))+' unique total randoms')
        else:
            print('did not find '+ffa)

    #fgu.sort('ZPOSS')
    #fgu['TILES'] = np.copy(fgu['TILE']).astype('<U100')
    #fu = unique(fgu,keys=['TARGETID'])#,keep='last')
    fu = fgu
    #fu.write(randir+str(rann)+'/rancomb_'+tp+'_Alltiles.fits',format='fits', overwrite=True)
    #return True
#     tiles = fgu['TILES']
#     tilesu = fu['TILES']
    #tlids = fgu['TILELOCIDS']
    #tlidsu = fu['TILELOCIDS']

#     for ii in range(0,len(tidsu)): #this takes a long time and something more efficient will be necessary
#         tid = tidsu[ii]#fu[ii]['TARGETID']
#         wt = tids == tid
#         ot = tilesu[ii]
#         #otl = tlidsu[ii]
#         tt = tiles[wt]
#         #tti = tlids[wt]
#         for tl in tt:
#             if tl != ot:
#                 tilesu[ii] += '-'+str(tl)
#         #for ti in tti:
#         #    if ti != otl:
#         #        tlidsu[ii] += '-'+str(ti)
#         if ii%1000 == 0:
#             print(ii)        
#     fu['TILES'] = tilesu
    #fu['TILELOCIDS'] = tlidsu
 
 
    fl = np.chararray(len(fu),unicode=True,itemsize=100)
    for ii in range(0,len(fu)):
        tl = fu['TILES'][ii]
        tls = tl.split('-')#.astype('int')
        tli = tls[0]
        if len(tls) > 1:
            #tls = tls.astype('int')
            tls.sort()
            tli = tls[0]
            for i in range(1,len(tls)):
                tli += '-'+tls[i]
        #else:
        #    tli = tls
        #print(tli)
        fl[ii] = tli   
    
    fu['TILES'] = fl
    print('number of unique tiles configurations '+str(len(np.unique(fu['TILES']))))


    NT = np.zeros(len(fgu))
    ros = np.zeros(len(fgu))
    print('counting tiles and finding rosette')
    for ii in range(0,len(fu['TILES'])): #not sure why, but this only works when using loop for Table.read but array option works for fitsio.read
        NT[ii] = np.char.count(fu['TILES'][ii],'-')+1
        ti = int(fu['TILES'][ii].split('-')[0])
        ros[ii] = tile2rosette(ti)
    fu['NTILE'] = NT    
            
    fu['rosette_number'] = ros
    print(np.unique(fu['rosette_number'],return_counts=True))

    fu.write(randir+str(rann)+'/rancomb_'+tp+'_Alltiles.fits',format='fits', overwrite=True)

def mkfullran(indir,rann,imbits,outf,tp,pd,bit,desitarg='SV3_DESI_TARGET',tsnr= 'TSNR2_ELG',notqso='',qsobit=4,fbcol='COADD_FIBERSTATUS'):

#     selz = dz['ZWARN'] != 999999
#     fs = dz[selz]
# 
#     #first, need to find locations to veto based data
#     nodata = fs["ZWARN_MTL"] & zwarn_mask["NODATA"] != 0
#     num_nod = np.sum(nodata)
#     print('number with no data '+str(num_nod))
#     badqa = fs["ZWARN_MTL"] & zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA") != 0
#     num_badqa = np.sum(badqa)
#     print('number with bad qa '+str(num_badqa))
#     nomtl = nodata & badqa
#     wfqa = ~nomtl
#     #wf = fs['FIBERSTATUS'] == 0
#     if specver == 'daily':
#         fbcol = 'FIBERSTATUS'
#     if specver == 'everest':
#         fbcol = 'COADD_FIBERSTATUS'
#     wf = fs[fbcol] == 0
#     print(len(fs[wf]),len(fs[wfqa]))
    zf = indir+'/datcomb_'+pd+'_tarspecwdup_zdone.fits'
    dz = Table.read(zf) 
    
    fs = get_specdat(indir,pd)
    stlid = 10000*fs['TILEID'] +fs['LOCATION']
    gtl = np.unique(stlid)

    wtype = ((dz[desitarg] & bit) > 0)
    if notqso == 'notqso':
        wtype &= ((dz[desitarg] & qsobit) == 0)

    wg = np.isin(dz['TILELOCID'],gtl)
    dz = dz[wtype&wg]
    print('length after selecting type and fiberstatus == 0 '+str(len(dz)))
    lznp = find_znotposs(dz)

    
    zf = indir+'/rancomb_'+str(rann)+pd+'wdupspec_zdone.fits'
    dz = Table.read(zf)
    #dz.remove_columns(['TILES','NTILE'])

    zfpd = indir+'/rancomb_'+str(rann)+pd+'_Alltilelocinfo.fits'
    dzpd = Table.read(zfpd)
    #dzpd.keep_columns(['TARGETID','TILES','NTILE'])
    dz = join(dz,dzpd,keys=['TARGETID'])
    #if maskzfail:
    #    wk = dz['ZPOSSNOTBAD'] == 1
    #else:
    #    wk = dz['ZPOSS'] == 1
    print('length before cutting to good positions '+str(len(dz)))
    wk = ~np.isin(dz['TILELOCID'],lznp)
    wk &= np.isin(dz['TILELOCID'],gtl)
    dz = dz[wk]    
    print('length after cutting to good positions '+str(len(dz)))
    dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/' 
    tcol = ['TARGETID','MASKBITS','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z'] #only including what are necessary for mask cuts for now
    #tcol = ['TARGETID','EBV','WISEMASK_W1','WISEMASK_W2','BRICKID','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G',\
    #'GALDEPTH_R','GALDEPTH_Z','PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','MASKBITS','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z']
    tarf = fitsio.read(dirrt+'/randoms-1-'+str(rann)+'.fits',columns=tcol)
    dz = join(dz,tarf,keys=['TARGETID'])
    del tarf
    dz = cutphotmask(dz,imbits)
    print('length after cutting to based on imaging veto mask '+str(len(dz)))
    dz.sort(tsnr) #should allow to later cut on tsnr for match to data
    dz = unique(dz,keys=['TARGETID'],keep='last')
    print('length after cutting to unique TARGETID '+str(len(dz)))
    print(np.unique(dz['NTILE']))
    
    dz.write(outf,format='fits', overwrite=True)
    del dz

def addcol_ran(fn,rann,dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/',ecol=['TARGETID','EBV','WISEMASK_W1','WISEMASK_W2','BRICKID','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']):
    dz = fitsio.read(fn)
    tarf = fitsio.read(dirrt+'/randoms-1-'+str(rann)+'.fits',columns=ecol)
    dz = join(dz,tarf,keys=['TARGETID'])
    dz.write(fn,format='fits', overwrite=True)
    del dz
        
    

def mkfulldat(zf,imbits,ftar,tp,bit,outf,ftiles,azf='',desitarg='DESI_TARGET',specver='daily',notqso='',qsobit=4):
    from scipy.special import erf
    #from desitarget.mtl import inflate_ledger
    if tp[:3] == 'BGS' or tp[:3] == 'MWS':
        pd = 'bright'        
        tscol = 'TSNR2_BGS'
    else:    
        pd = 'dark'
        tscol = 'TSNR2_ELG'
    #fs = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+specver+'/datcomb_'+pd+'_spec_zdone.fits')
#     dz = Table.read(zf) 
#     selz = dz['ZWARN_MTL'] != 999999
#     fs = dz[selz]
#     nodata = fs["ZWARN_MTL"] & zwarn_mask["NODATA"] != 0
#     num_nod = np.sum(nodata)
#     print('number with no data '+str(num_nod))
#     badqa = fs["ZWARN"] & zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA") != 0
#     num_badqa = np.sum(badqa)
#     print('number with bad qa '+str(num_badqa))
#     nomtl = nodata & badqa
#     wfqa = ~nomtl
#     #wf = fs['FIBERSTATUS'] == 0
#     if specver == 'daily':
#         fbcol = 'FIBERSTATUS'
#     if specver == 'everest':
#         fbcol = 'COADD_FIBERSTATUS'
#     wf = fs[fbcol] == 0
#     print(len(fs[wf]),len(fs[wfqa]))
    
    indir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+specver
    fs = get_specdat(indir,pd)
    stlid = 10000*fs['TILEID'] +fs['LOCATION']
    gtl = np.unique(stlid)

    dz = Table.read(zf) 
    
    
    wtype = ((dz[desitarg] & bit) > 0)
    if notqso == 'notqso':
        print('removing QSO targets')
        wtype &= ((dz[desitarg] & qsobit) == 0)

    wg = np.isin(dz['TILELOCID'],gtl)
    print(len(dz[wtype]))
    print(len(dz[wg]))
    dz = dz[wtype&wg]
    print('length after selecting type and good hardware '+str(len(dz)))
    lznp = find_znotposs(dz)
    wk = ~np.isin(dz['TILELOCID'],lznp)#dz['ZPOSS'] == 1
    dz = dz[wk]
    print('length after priority veto '+str(len(dz)))
    print('joining to full imaging')
    dz.remove_columns(['RA','DEC','DESI_TARGET','BGS_TARGET']) #these come back in with merge to full target file
    dz = join(dz,ftar,keys=['TARGETID'])
    #print('length after join to full targets (should be same) '+str(len(dz)))
    dz = cutphotmask(dz,imbits)
    dtl = Table.read(ftiles)
    dtl.keep_columns(['TARGETID','NTILE','TILES','TILELOCIDS'])
    dz = join(dz,dtl,keys='TARGETID')
    wz = dz['ZWARN'] != 999999 #this is what the null column becomes
    wz &= dz['ZWARN']*0 == 0 #just in case of nans
    dz['LOCATION_ASSIGNED'] = np.zeros(len(dz)).astype('bool')
    dz['LOCATION_ASSIGNED'][wz] = 1
    tlids = np.unique(dz['TILELOCID'][wz])
    wtl = np.isin(dz['TILELOCID'],tlids)
    dz['TILELOCID_ASSIGNED'] = 0
    dz['TILELOCID_ASSIGNED'][wtl] = 1
    print('number of unique targets at assigned tilelocid:')
    print(len(np.unique(dz[wtl]['TARGETID'])))

    if tp[:3] == 'ELG' and azf != '':# or tp == 'ELG_HIP':
        arz = fitsio.read(azf,columns=['TARGETID','LOCATION','TILEID','OII_FLUX','OII_FLUX_IVAR','SUBSET','DELTACHI2'])
        st = []
        for i in range(0,len(arz)):
            st.append(arz['SUBSET'][i][:4])
        st = np.array(st)
        #wg = arz[fbcol] == 0
        wg = st == "thru"
        arz = arz[wg]
        o2c = np.log10(arz['OII_FLUX'] * np.sqrt(arz['OII_FLUX_IVAR']))+0.2*np.log10(arz['DELTACHI2'])
        w = (o2c*0) != 0
        w |= arz['OII_FLUX'] < 0
        o2c[w] = -20
        #arz.keep_columns(['TARGETID','LOCATION','TILEID','o2c','OII_FLUX','OII_SIGMA'])#,'Z','ZWARN','TSNR2_ELG'])    
        arz = Table(arz)
        arz['o2c'] = o2c
        dz = join(dz,arz,keys=['TARGETID','LOCATION','TILEID'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['', '_OII'])
        
        dz.remove_columns(['SUBSET','DELTACHI2_OII'])#,fbcol+'_OII'])
        print('check length after merge with OII strength file:' +str(len(dz)))

    if tp[:3] == 'QSO' and azf != '':
        arz = Table.read(azf)
        arz.keep_columns(['TARGETID','LOCATION','TILEID','Z','ZERR','Z_QN'])
        print(arz.dtype.names)
        #arz['TILE'].name = 'TILEID'
        dz = join(dz,arz,keys=['TARGETID','TILEID','LOCATION'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
        dz['Z'].name = 'Z_RR' #rename the original redrock redshifts
        dz['Z_QF'].name = 'Z' #the redshifts from the quasar file should be used instead


    dz['sort'] = dz['LOCATION_ASSIGNED']*dz[tscol]+dz['TILELOCID_ASSIGNED']
    dz.sort('sort')
    dz = unique(dz,keys=['TARGETID'],keep='last')
    if tp == 'ELG' or tp == 'ELG_HIP':
        print('number of masked oII row (hopefully matches number not assigned) '+ str(np.sum(dz['o2c'].mask)))
    print('length after cutting to unique targetid '+str(len(dz)))
    print('LOCATION_ASSIGNED numbers')
    print(np.unique(dz['LOCATION_ASSIGNED'],return_counts=True))
   
    print('TILELOCID_ASSIGNED numbers')
    print(np.unique(dz['TILELOCID_ASSIGNED'],return_counts=True))
    #print('length after join to file with tiles info is '+str(len(dz)))
    #NT = np.zeros(len(dz))
    ros = np.zeros(len(dz))
    #ti = np.zeros(len(dz))

    probl = np.zeros(len(dz))
    #dr = fitsio.read(e2eout+ program+'/'+type+'_oneper_full.ran.fits')

    #get completeness based on unique sets of tiles
    compa = []
    tll = []
    ti = 0
    print('getting completenes')
    dz.sort('TILES')
    nts = len(np.unique(dz['TILES']))
    tlsl = dz['TILES']
    tlslu = np.unique(tlsl)
    laa = dz['LOCATION_ASSIGNED']
    
    #for tls in np.unique(dz['TILES']): #this is really slow now, need to figure out a better way
    i = 0
    while i < len(dz):
        tls  = []
        tlis = []
        nli = 0
        nai = 0
    
        while tlsl[i] == tlslu[ti]:
            nli += 1
            nai += laa[i]
            i += 1
            if i == len(dz):
                break
    
        if ti%1000 == 0:
            print('at tiles '+str(ti)+' of '+str(nts))

        #w = dz['TILES'] == tls
        #no = sum(dz[w]['LOCATION_ASSIGNED'])
        #nt = len(dz[w])
        cp = nai/nli#no/nt
        #print(tls,cp,no,nt)
        compa.append(cp)
        tll.append(tlslu[ti])
        ti += 1
    comp_dicta = dict(zip(tll, compa))
    fcompa = []
    for tl in dz['TILES']:
        fcompa.append(comp_dicta[tl]) 
    dz['COMP_TILE'] = np.array(fcompa)
    wc0 = dz['COMP_TILE'] == 0
    print('number of targets in 0 completeness regions '+str(len(dz[wc0])))       




    locl,nlocl = np.unique(dz['TILELOCID'],return_counts=True)
    #wa = dzz['LOCATION_ASSIGNED'] == 1
    #if len(dzz[wa]) != len(dzz):
     #   print('!found some zwarn = 0 without location_assigned = 1!')
    wz = dz['LOCATION_ASSIGNED'] == 1
    dzz = dz[wz]

    loclz,nloclz = np.unique(dzz['TILELOCID'],return_counts=True)
    print(np.max(nloclz),np.min(loclz))
    #print(np.histogram(nloclz))
    print(len(locl),len(nloclz),sum(nlocl),sum(nloclz))
    natloc = ~np.isin(dz['TILELOCID'],loclz)
    print('number of unique targets left around unassigned locations is '+str(np.sum(natloc)))
    locs = np.copy(dz['TILELOCID'])
# 
# 
    print('reassigning TILELOCID for duplicates and finding rosette')
    nch = 0
    nbl = 0
    tlids = dz['TILELOCIDS']
#     nf = 0
#     #dz.write('temp.fits',format='fits', overwrite=True)
#     #fdz = fitsio.read('temp.fits')
    for ii in range(0,len(dz['TILEID'])): #not sure why, but this only works when using loop for Table.read but array option works for fitsio.read
#         NT[ii] = np.char.count(dz['TILES'][ii],'-')+1
#         #ti[ii] = int(dz['TILE'][ii].split('-')[0])
#         tiles = dz['TILES'][ii].split('-')
#         ti = int(tiles[0])
        ti = dz[ii]['TILEID']
        
        if natloc[ii]:# == False:
            nbl += 1
            s = 0
            tids = tlids[ii].split('-')
            if s == 0:
                for tl in tids:
                    ttlocid  = int(tl)              
                    if np.isin(ttlocid,loclz):
                        #dz[ii]['TILELOCID'] = ttlocid
                        locs[ii] = ttlocid #use below instead and assign at end, maybe faster
                        nch += 1
                        s = 1
                        break
        if ii%10000 == 0:
            print(ii,len(dz['TILEID']),ti,ros[ii],nch,nbl)
     
    dz['TILELOCID'] = locs
    locl,nlocl = np.unique(dz['TILELOCID'],return_counts=True)
#     #wa = dzz['LOCATION_ASSIGNED'] == 1
#     #if len(dzz[wa]) != len(dzz):
#      #   print('!found some zwarn = 0 without location_assigned = 1!')
    loclz,nloclz = np.unique(dzz['TILELOCID'],return_counts=True)
#     print(np.max(nloclz),np.min(loclz))
#     #print(np.histogram(nloclz))
#     print(len(locl),len(nloclz),sum(nlocl),sum(nloclz))

    #NT = np.char.count(dz['TILE'],'-')
    #NT += 1
    print(np.unique(dz['NTILE']))

    #get tilelocid probs
    #wz = dz['ZWARN'] == 0
    print('getting fraction assigned for each tilelocid')
    nm = 0
    nmt =0
    pd = []
    nloclt = len(locl)
    for i in range(0,len(locl)):
        if i%10000 == 0:
            print('at row '+str(i)+' of '+str(nloclt))
        nt = nlocl[i]
        loc = locl[i]
        w = loclz == loc
        nz = 0
        if len(loclz[w]) == 1:
            nz = nloclz[w] #these are supposed all be 1...
            
        else:
            #print(loclz[w],nt) 
            nm += 1.
            nmt += nt
        if len(loclz[w]) > 1:
            print('why is len(loclz[w]) > 1?')
            #wa = dz['TILELOCID'] == loc
            #print(nz,nt,len(dz[wa]),len(loclz[w]),len(nloclz[w]),len(nz),nloclz[w])
            #probl[wa] = nz/nt
            #pd.append((loc,nz/nt)) 
        pd.append((loc,nz/nt))  
    pd = dict(pd)
    for i in range(0,len(dz)):
        probl[i] = pd[dz['TILELOCID'][i]]
    print('number of fibers with no observation, number targets on those fibers')
    print(nm,nmt)
    
    #print(np.min(probl),np.max(probl))
    #dz = Table.read(zf) #table is slow, so using fitsio above, Table here
    dz['FRACZ_TILELOCID'] = probl
    print('sum of 1/FRACZ_TILELOCID, 1/COMP_TILE, and length of input; should match')
    print(np.sum(1./dz[wz]['FRACZ_TILELOCID']),np.sum(1./dz[wz]['COMP_TILE']),len(dz))
    #print(np.unique(dz['TILE']))
    #dz['NTILE']  = NT
    dz['WEIGHT_ZFAIL'] = np.ones(len(dz))
            
    print(np.unique(dz['NTILE']))
    dz.write(outf,format='fits', overwrite=True)

def mkclusdat(fl,weighttileloc=True,zmask=False,tp='',dchi2=9,tsnrcut=80,rcut=None,ntilecut=0,ccut=None,ebits=None):
    '''
    fl is the root of the input/output file
    weighttileloc determines whether to include 1/FRACZ_TILELOCID as a completeness weight
    zmask determines whether to apply a mask at some given redshift
    tp is the target type
    dchi2 is the threshold for keeping as a good redshift
    tnsrcut determines where to mask based on the tsnr2 value (defined below per tracer)

    '''    
    ff = Table.read(fl+'full_noveto.dat.fits')
    
    
    if ebits is not None:
        print('number before imaging mask '+str(len(ff)))
        ff = cutphotmask(ff,ebits)
        print('number after imaging mask '+str(len(ff)))
    ff.write(fl+'full.dat.fits',overwrite=True,format='fits')
    wzm = ''
    if zmask:
        wzm = 'zmask_'
    if rcut is not None:
        wzm += 'rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
    if ntilecut > 0:
        wzm += 'ntileg'+str(ntilecut)+'_'    
    if ccut is not None:
        wzm += ccut+'_' #you could change this to however you want the file names to turn out
    outf = fl+wzm+'clustering.dat.fits'
    wz = ff['ZWARN'] == 0
    print('length before cutting to objects with redshifts '+str(len(ff)))
    print('length after cutting to zwarn == 0 '+str(len(ff[wz])))
    if tp == 'QSO':
        #good redshifts are currently just the ones that should have been defined in the QSO file when merged in full
        wz = ff['Z']*0 == 0
        wz &= ff['Z'] != 999999
        wz &= ff['Z'] != 1.e20
        wz &= ff['ZWARN'] != 999999
        wz &= ff['TSNR2_QSO'] > tsnrcut
    
    if tp[:3] == 'ELG':
        wz = ff['o2c'] > dchi2
        wz &= ff['ZWARN']*0 == 0
        wz &= ff['ZWARN'] != 999999
        print('length after oII cut '+str(len(ff[wz])))
        wz &= ff['LOCATION_ASSIGNED'] == 1
        print('length after also making sure location assigned '+str(len(ff[wz])))
        wz &= ff['TSNR2_ELG'] > tsnrcut
        print('length after tsnrcut '+str(len(ff[wz])))
    if tp == 'LRG':
        print('applying extra cut for LRGs')
        # Custom DELTACHI2 vs z cut from Rongpu
        drz = (10**(3 - 3.5*ff['Z']))
        mask_bad = (drz>30) & (ff['DELTACHI2']<30)
        mask_bad |= (drz<30) & (ff['DELTACHI2']<drz)
        mask_bad |= (ff['DELTACHI2']<10)
        wz &= ff['Z']<1.4
        wz &= (~mask_bad)

        #wz &= ff['DELTACHI2'] > dchi2
        print('length after Rongpu cut '+str(len(ff[wz])))
        wz &= ff['TSNR2_ELG'] > tsnrcut
        print('length after tsnrcut '+str(len(ff[wz])))

    if tp[:3] == 'BGS':
        print('applying extra cut for BGS')
        wz &= ff['DELTACHI2'] > dchi2
        print('length after dchi2 cut '+str(len(ff[wz])))
        wz &= ff['TSNR2_BGS'] > tsnrcut
        print('length after tsnrcut '+str(len(ff[wz])))

    
    ff = ff[wz]
    print('length after cutting to good z '+str(len(ff)))
    ff['WEIGHT'] = ff['WEIGHT_ZFAIL']
    if weighttileloc == True:
        ff['WEIGHT'] *= 1./ff['FRACZ_TILELOCID']

    #weights for imaging systematic go here
    ff['WEIGHT_SYS'] =  1.
#     if tp[:3] == 'ELG':
#         zmin = 0.8
#         zmax = 1.5
#         selz = ff['Z'] > zmin
#         selz &= ff['Z'] < zmax
#         ec = ff[selz]
#         hd = np.histogram(ec['EBV'],weights=1./ec['COMP_TILE'],range=(0,.15))
#         fer = fitsio.read(fl+'0_full.ran.fits')
#         hr = np.histogram(fer['EBV'],bins=hd[1])
#         norm = sum(hr[0])/sum(hd[0])
#         xl = hd[1][:-1]+(hd[1][1]-hd[1][0])/2.
#         yl = hd[0]/hr[0]*norm
#         el = np.sqrt(hd[0])/hr[0]*norm
#         m,b = np.polyfit(xl,yl,1,w=1/el)
#         print('linear fits coefficients to EBV are '+str(m)+' '+str(b))
#         ff['WEIGHT_SYS'] = 1./(m*ff['EBV']+b)
#         hd = np.histogram(np.log(ec['GALDEPTH_G']),weights=1./ec['COMP_TILE'],range=(5.5,8.))
#         hr = np.histogram(np.log(fer['GALDEPTH_G']),bins=hd[1])
#         norm = sum(hr[0])/sum(hd[0])
#         xl = hd[1][:-1]+(hd[1][1]-hd[1][0])/2.
#         yl = hd[0]/hr[0]*norm
#         el = np.sqrt(hd[0])/hr[0]*norm
#         m,b = np.polyfit(xl,yl,1,w=1/el)
#         print('linear fits coefficients to GALDEPTH_G are '+str(m)+' '+str(b))
#         ff['WEIGHT_SYS'] *= 1./(m*np.log(ff['GALDEPTH_G'])+b)

    ff['WEIGHT'] *= ff['WEIGHT_SYS']

    if zmask:
        whz = ff['Z'] < 1.6
        ff = ff[whz]

        fzm = fitsio.read('/global/homes/m/mjwilson/desi/DX2DROPOUT/radial_mask.fits')
        zma = []
        for z in ff['Z']:
            zind = int(z/1e-6)
            zma.append(fzm[zind]['RADIAL_MASK'])        
        zma = np.array(zma)
        wm = zma == 0
        ff = ff[wm]    
    #apply any cut on rosette radius
    if rcut is not None:
        wr = ff['rosette_r'] > rcut[0]
        wr &= ff['rosette_r'] <  rcut[1]
        print('length before rosette radius cut '+str(len(ff)))
        ff = ff[wr]
        print('length after rosette radius cut '+str(len(ff)))
    #apply cut on ntile
    if ntilecut > 0:
        print('length before ntile cut '+str(len(ff)))
        wt = ff['NTILE'] > ntilecut
        ff = ff[wt]
        print('length after ntile cut '+str(len(ff)))    
    if ccut == 'notQSO':
        wc = (ff['SV3_DESI_TARGET'] & sv3_targetmask.desi_mask['QSO']) ==  0
        print('length before cutting to not QSO '+str(len(ff)))
        ff = ff[wc]
        print('length after cutting to not QSO '+str(len(ff)))
    if ccut == 'zQSO':
        wc = ff['SPECTYPE'] ==  'QSO'
        print('length before cutting to spectype QSO '+str(len(ff)))
        ff = ff[wc]
        print('length after cutting to spectype QSO '+str(len(ff)))

    #select down to specific columns below and then also split N/S
    wn = ff['PHOTSYS'] == 'N'
    kl = ['RA','DEC','Z','WEIGHT','TARGETID','NTILE','TILES','WEIGHT_SYS']
    if tp[:3] == 'BGS':
        ff['flux_r_dered'] = ff['FLUX_R']/ff['MW_TRANSMISSION_R']
        kl.append('flux_r_dered')
        print(kl)
    ff.keep_columns(kl)
    print('minimum,maximum weight')
    print(np.min(ff['WEIGHT']),np.max(ff['WEIGHT']))
    ff.write(outf,format='fits', overwrite=True)
    outfn = fl+wzm+'N_clustering.dat.fits'
    ff[wn].write(outfn,format='fits', overwrite=True)
    outfn = fl+wzm+'S_clustering.dat.fits'
    ffs = ff[~wn]
    ffs.write(outfn,format='fits', overwrite=True)
    for reg in ['DS','DN']: #split DECaLS NGC/SGC
        outfn = fl+wzm+reg+'_clustering.dat.fits'
        sel = densvar.sel_reg(ffs['RA'],ffs['DEC'],reg)
        ffs[sel].write(outfn,format='fits', overwrite=True)

def mkclusran(fl,rann,rcols=['Z','WEIGHT'],zmask=False,tsnrcut=80,tsnrcol='TSNR2_ELG',ebits=None):
    #first find tilelocids where fiber was wanted, but none was assigned; should take care of all priority issues
    wzm = ''
    if zmask:
        wzm = 'zmask_'

    #ffd = Table.read(fl+'full.dat.fits')
    fcd = Table.read(fl+wzm+'clustering.dat.fits')
    ffr = Table.read(fl+str(rann)+'_full_noveto.ran.fits')
    if ebits is not None:
        print('number before imaging mask '+str(len(ffr)))
        ffr = cutphotmask(ffr,ebits)
        print('number after imaging mask '+str(len(ffr)))
    ffr.write(fl+str(rann)+'_full.ran.fits',overwrite=True,format='fits')

    #if type[:3] == 'ELG' or type == 'LRG':
    wz = ffr[tsnrcol] > tsnrcut
    #wif = np.isin(ffr['TILELOCID'],ffd['TILELOCID'])
    #wic = np.isin(ffr['TILELOCID'],fcd['TILELOCID'])
    #wb = wif & ~wic #these are the tilelocid in the full but not in clustering, should be masked
    #ffc = ffr[~wb]
    ffc = ffr[wz]
    print('length after,before tsnr cut:')
    print(len(ffc),len(ffr))
    inds = np.random.choice(len(fcd),len(ffc))
    dshuf = fcd[inds]

    for col in rcols: 
        ffc[col] = dshuf[col] 
    wn = ffc['PHOTSYS'] == 'N'
    ffc.keep_columns(['RA','DEC','Z','WEIGHT','TARGETID','NTILE','TILES'])  
    outf =  fl+wzm+str(rann)+'_clustering.ran.fits' 
    ffc.write(outf,format='fits', overwrite=True)

    outfn =  fl+wzm+'N_'+str(rann)+'_clustering.ran.fits' 
    fcdn = Table.read(fl+wzm+'N_clustering.dat.fits')
    ffcn = ffc[wn]
    inds = np.random.choice(len(fcdn),len(ffcn))
    dshuf = fcdn[inds]
    for col in rcols: 
        ffcn[col] = dshuf[col]     
    ffcn.write(outfn,format='fits', overwrite=True)

    outfs =  fl+wzm+'S_'+str(rann)+'_clustering.ran.fits' 
    fcds = Table.read(fl+wzm+'S_clustering.dat.fits')
    ffcs = ffc[~wn]
    inds = np.random.choice(len(fcds),len(ffcs))
    dshuf = fcds[inds]
    for col in rcols: 
        ffcs[col] = dshuf[col]     
    ffcs.write(outfs,format='fits', overwrite=True)

    for reg in ['DS','DN']: #split DECaLS NGC/SGC
        outfn = fl+wzm+reg+'_'+str(rann)+'_clustering.ran.fits'
        sel = densvar.sel_reg(ffcs['RA'],ffcs['DEC'],reg)
        fcd = Table.read(fl+wzm+reg+'_clustering.dat.fits')
        ffss = ffcs[sel]
        inds = np.random.choice(len(fcd),len(ffss))
        dshuf = fcd[inds]
        for col in rcols: 
            ffss[col] = dshuf[col]     

        ffss.write(outfn,format='fits', overwrite=True)





    

def randomtiles_allmain(tiles,dirout='/global/cfs/cdirs/desi/survey/catalogs/main/LSS/random',imin=0,imax=18,rann=1,dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/' ):
    '''
    tiles should be a table containing the relevant info
    '''
    trad = desimodel.focalplane.get_tile_radius_deg()*1.1 #make 10% greater just in case
    
    for ii in range(imin,imax):
        print(trad,ii)
        rt = fitsio.read(dirrt+'/randoms-'+str(rann)+'-'+str(ii)+'.fits')
        print('loaded random file '+str(ii)) 
    
        for i in range(0,len(tiles)):
            
            #print('length of tile file is (expected to be 1):'+str(len(tiles)))
            tile = tiles['TILEID'][i]
            fname = dirout+str(ii)+'/tilenofa-'+str(tile)+'.fits'
            if os.path.isfile(fname):
                print(fname +' already exists')
            else:
                tdec = tiles['DEC'][i]
                decmin = tdec - trad
                decmax = tdec + trad
                wdec = (rt['DEC'] > decmin) & (rt['DEC'] < decmax)
                print(len(rt[wdec]))
                inds = desimodel.footprint.find_points_radec(tiles['RA'][i], tdec,rt[wdec]['RA'], rt[wdec]['DEC'])
                print('got indexes')
                rtw = rt[wdec][inds]
                rmtl = Table(rtw)
                #rmtl['TARGETID'] = np.arange(len(rmtl))
                print(len(rmtl['TARGETID'])) #checking this column is there
                rmtl['DESI_TARGET'] = np.ones(len(rmtl),dtype=int)*2
                rmtl['NUMOBS_INIT'] = np.zeros(len(rmtl),dtype=int)
                rmtl['NUMOBS_MORE'] = np.ones(len(rmtl),dtype=int)
                rmtl['PRIORITY'] = np.ones(len(rmtl),dtype=int)*3400
                rmtl['OBSCONDITIONS'] = np.ones(len(rmtl),dtype=int)*516#tiles['OBSCONDITIONS'][i]
                rmtl['SUBPRIORITY'] = np.random.random(len(rmtl))
                rmtl.write(fname,format='fits', overwrite=True)
                print('added columns, wrote to '+fname)

def randomtiles_main_fromran(tiles,rt,rann,dirout='/global/cfs/cdirs/desi/survey/catalogs/main/LSS/random' ):
    '''
    tiles should be a table containing the relevant info
    take the input random, rt, as an argument so when doing in parallel only one copy in memory
    '''
    trad = desimodel.focalplane.get_tile_radius_deg()*1.1 #make 10% greater just in case
    print(trad)
    
    for i in range(0,len(tiles)):
        
        #print('length of tile file is (expected to be 1):'+str(len(tiles)))
        tile = tiles['TILEID'][i]
        fname = dirout+str(rann)+'/tilenofa-'+str(tile)+'.fits'
        if os.path.isfile(fname):
            print(fname +' already exists')
        else:
            tdec = tiles['DEC'][i]
            decmin = tdec - trad
            decmax = tdec + trad
            wdec = (rt['DEC'] > decmin) & (rt['DEC'] < decmax)
            print(len(rt[wdec]))
            inds = desimodel.footprint.find_points_radec(tiles['RA'][i], tdec,rt[wdec]['RA'], rt[wdec]['DEC'])
            print('got indexes')
            rtw = rt[wdec][inds]
            rmtl = Table(rtw)
            #rmtl['TARGETID'] = np.arange(len(rmtl))
            print(len(rmtl['TARGETID'])) #checking this column is there
            rmtl['DESI_TARGET'] = np.ones(len(rmtl),dtype=int)*2
            rmtl['NUMOBS_INIT'] = np.zeros(len(rmtl),dtype=int)
            rmtl['NUMOBS_MORE'] = np.ones(len(rmtl),dtype=int)
            rmtl['PRIORITY'] = np.ones(len(rmtl),dtype=int)*3400
            rmtl['OBSCONDITIONS'] = np.ones(len(rmtl),dtype=int)*516#tiles['OBSCONDITIONS'][i]
            rmtl['SUBPRIORITY'] = np.random.random(len(rmtl))
            rmtl.write(fname,format='fits', overwrite=True)
            print('added columns, wrote to '+fname)
