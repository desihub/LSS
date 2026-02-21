'''
python functions to do various useful date processing/manipulation
'''
import numpy as np
from scipy.special import erf
import fitsio
import glob
import os
import astropy.io.fits as fits
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt
import desimodel.footprint as foot
import desimodel.focalplane
from random import random
from LSS import common_tools as common
from LSS import ssr_tools
from desitarget.io import read_targets_in_tiles
from desitarget.sv3 import sv3_targetmask
import datetime

from LSS.Cosmo import distance

def tile2rosette(tile):
    if tile < 433:
        return (tile-1)//27
    else:
        if tile >= 433 and tile < 436:
            return 13
        if tile >= 436 and tile < 439:
            return 14
        if tile >= 439 and tile < 442:
            return 15
        if tile >= 442 and tile <=480:
            return (tile-442)//3
            
        if tile > 480:
            return tile//30    
    return 999999 #shouldn't be any more?

def calc_rosr(rosn,ra,dec):
    #given rosetter number and ra,dec, calculate distance from center 
    roscen = {0:(150.100,2.182),1:(179.6,0),2:(183.1,0),3:(189.9,61.8),4:(194.75,28.2)\
    ,5:(210.0,5.0),6:(215.5,52.5),7:(217.8,34.4),8:(216.3,-0.6),9:(219.8,-0.6)\
    ,10:(218.05,2.43),11:(242.75,54.98),12:(241.05,43.45),13:(245.88,43.45),14:(252.5,34.5)\
    ,15:(269.73,66.02),16:(194.75,24.7),17:(212.8,-0.6),18:(269.73,62.52),19:(236.1,43.45)}
    ra = ra*np.pi/180.
    dec = dec*np.pi/180.
    rac,decc = roscen[rosn]
    rac = rac*np.pi/180.
    decc = decc*np.pi/180.
    cd = np.sin(dec)*np.sin(decc)+np.cos(dec)*np.cos(decc)*np.cos(rac-ra)
    ad = np.arccos(cd)*180./np.pi
    if ad > 2.5:
        print(rosn,ra,dec,rac,decc)
    return ad
    
def combtile_spec(tiles,outf='',rel='daily'):
    s = 0
    n = 0
    if os.path.isfile(outf):
        specd = Table.read(outf)
        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)
    else:
        tmask = np.ones(len(tiles)).astype('bool')    

    for tile,zdate in zip(tiles[tmask]['TILEID'],tiles[tmask]['LASTNIGHT']):
        zdate = str(zdate)
        tspec = combspecdata(tile,zdate,rel=rel)
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
    if n > 0:
        specd.write(outf,format='fits', overwrite=True)       
    else:
        print('There was nothing to have done') 

def combspecdata(tile,zdate,specroot='/global/cfs/cdirs/desi/spectro/redux/',rel='daily' ):
    #put data from different spectrographs together, one table for fibermap, other for z
    coaddir=specroot+rel+'/tiles/cumulative/'
    specs = []
    #find out which spectrograph have data
    for si in range(0,10):
        
        try:
            ff = coaddir+str(tile)+'/'+zdate+'/zbest-'+str(si)+'-'+str(tile)+'-thru'+zdate+'.fits'
            fitsio.read(ff)

            specs.append(si)
        except:
            print('no spectrograph '+str(si)+ ' for tile '+str(tile))
            #print(ff)
    print('spectrographs with data:')
    print(specs)            
    if len(specs) == 0:
        return None
    tspec = Table.read(coaddir+str(tile)+'/'+zdate+'/zbest-'+str(specs[0])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='ZBEST')
    tf = Table.read(coaddir+str(tile)+'/'+zdate+'/zbest-'+str(specs[0])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='FIBERMAP')
    ts = Table.read(coaddir+str(tile)+'/'+zdate+'/coadd-'+str(specs[0])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='SCORES')
    for i in range(1,len(specs)):
        tn = Table.read(coaddir+str(tile)+'/'+zdate+'/zbest-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='ZBEST')
        tnf = Table.read(coaddir+str(tile)+'/'+zdate+'/zbest-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='FIBERMAP')
        try:
            tns = Table.read(coaddir+str(tile)+'/'+zdate+'/coadd-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='SCORES')
            ts = vstack([ts,tns],metadata_conflicts='silent')
        except:
            print('did not find '+coaddir+str(tile)+'/'+zdate+'/coadd-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits')
        tspec = vstack([tspec,tn],metadata_conflicts='silent')
        tf = vstack([tf,tnf],metadata_conflicts='silent')
        
    
    tf = unique(tf,keys=['TARGETID'])
    #tf.keep_columns(['FIBERASSIGN_X','FIBERASSIGN_Y','TARGETID','LOCATION','FIBER','FIBERSTATUS','PRIORITY','FA_TARGET','FA_TYPE',\
    #'OBJTYPE','DELTA_X','DELTA_Y','PSF_TO_FIBER_SPECFLUX','EXPTIME','OBJTYPE','NIGHT','EXPID','MJD','SV3_DESI_TARGET','SV3_BGS_TARGET'])
    tspec = join(tspec,tf,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')
    tspec = join(tspec,ts,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')
    print(len(tspec),len(tf))
    #tspec['LOCATION'] = tf['LOCATION']
    #tspec['FIBERSTATUS'] = tf['FIBERSTATUS']
    #tspec['PRIORITY'] = tf['PRIORITY']
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

# def cutphotmask(aa,bits):
#     print(str(len(aa)) +' before imaging veto' )
#     keep = (aa['NOBS_G']>0) & (aa['NOBS_R']>0) & (aa['NOBS_Z']>0)
#     for biti in bits:
#         keep &= ((aa['MASKBITS'] & 2**biti)==0)
#     aa = aa[keep]
#     print(str(len(aa)) +' after imaging veto' )
#     return aa

def combtiles_wdup(tiles,mdir='',fout='',tarcol=['RA','DEC','TARGETID','SV3_DESI_TARGET','SV3_BGS_TARGET','SV3_MWS_TARGET','SUBPRIORITY','PRIORITY_INIT','TARGET_STATE','TIMESTAMP','ZWARN','PRIORITY']):
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
        wt = tiles['TILEID'] == tile
        #tars = read_targets_in_tiles(mdir,tiles[wt],mtl=True,isodate=fht['MTLTIME'])
        tars = read_targets_in_tiles(mdir,tiles[wt],mtl=True,isodate=fht['MTLTIME'],columns=tarcol)
        #tars.keep_columns(tarcols)
        #tars = tars[[b for b in tarcol]]
        
        tt = Table.read(faf,hdu='POTENTIAL_ASSIGNMENTS')
        tars = join(tars,tt,keys=['TARGETID'])
        tars['TILEID'] = tile
        tars['ZWARN'].name = 'ZWARN_MTL'
        if s == 0:
            tarsn = tars
            s = 1
        else:
            tarsn = vstack([tarsn,tars],metadata_conflicts='silent')
        tarsn.sort('TARGETID')
        n += 1
        print(tile,n,len(tiles[tmask]),len(tarsn)) 
    if os.path.isfile('tmp.fits'):
        os.system('rm tmp.fits')
    fd = fitsio.FITS('tmp.fits', "rw")
    fd.write(np.array(tarsn),extname='POTENTIAL_ASSINGMENT')
    fd['POTENTIAL_ASSINGMENT'].write_comment("concatenation of SV3 POTENTIAL_ASSIGNMENT information, joined to columns in targe files")
    fd['POTENTIAL_ASSINGMENT'].write_history("updated on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    fd.close()    
    os.system('mv tmp.fits '+fout)
    #tarsn.write(fout,format='fits', overwrite=True)       

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

#moved to common_tools 
#def find_znotposs(dz):
# 
#     dz.sort('TARGETID')
#     tidnoz = []
#     tids = np.unique(dz['TARGETID'])
#     ti = 0
#     i = 0
#     
#     print('finding targetids that were not observed')
#     while i < len(dz):
#         za = 0
#     
#         while dz[i]['TARGETID'] == tids[ti]:
#             if dz[i]['ZWARN'] != 999999:
#                 za = 1
#                 #break
#             i += 1
#             if i == len(dz):
#                 break
#         if za == 0:
#             tidnoz.append(tids[ti])
#       
#         if ti%30000 == 0:
#             print(ti)
#         ti += 1 
# 
#     
#     selnoz = np.isin(dz['TARGETID'],tidnoz)
#     tidsb = np.unique(dz[selnoz]['TILELOCID'])
#     #dz = dz[selnoz]
#     dz.sort('TILELOCID')
#     tids = np.unique(dz['TILELOCID'])
#     print('number of targetids with no obs '+str(len(tidnoz)))
#     tlidnoz = []
#     lznposs = []
#     
#     ti = 0
#     i = 0
#     
#     while i < len(dz):
#         za = 0
#     
#         while dz[i]['TILELOCID'] == tids[ti]:
#             if dz[i]['ZWARN'] != 999999:
#                 za = 1
#                 #break
#             i += 1
#             if i == len(dz):
#                 break
#         if za == 0:
#             tlidnoz.append(tids[ti])
#             #if np.isin(tids[ti],tidsb):
#             #    lznposs.append(tids[ti])
#       
#         if ti%30000 == 0:
#             print(ti,len(tids))
#         ti += 1 
#     #the ones to veto are now the join of the two
#     wtbtlid = np.isin(tlidnoz,tidsb)
#     tlidnoz = np.array(tlidnoz)
#     lznposs = tlidnoz[wtbtlid]
#     print('number of locations where assignment was not possible because of priorities '+str(len(lznposs)))
#     return lznposs
    
def count_tiles_better(fs,dr,pd,rann=0,specrel='daily',fibcol='COADD_FIBERSTATUS'):
    '''
    from files with duplicates that have already been sorted by targetid, quickly go 
    through and get the multi-tile information
    dr is either 'dat' or 'ran'
    returns file with TARGETID,NTILE,TILES,TILELOCIDS
    '''
    
    #fs = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'+specrel+'/datcomb_'+pd+'_specwdup_Alltiles.fits')
    #wf = fs['FIBERSTATUS'] == 0
    wf = fs[fibcol] == 0
    stlid = 10000*fs['TILEID'] +fs['LOCATION']
    gtl = np.unique(stlid[wf])
    
    if dr == 'dat':
        fj = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'+specrel+'/datcomb_'+pd+'_tarspecwdup_Alltiles.fits')
        #outf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/datcomb_'+pd+'ntileinfo.fits' 
    if dr == 'ran':
        fj = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'+specrel+'/rancomb_'+str(rann)+pd+'wdupspec_Alltiles.fits')
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


def combran_wdup(tiles,rann,randir,tp,sv3dir,specf,keepcols=[]):

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
    #specf = Table.read(sv3dir+'datcomb_'+tp+'_specwdup_Alltiles.fits')
    specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    specf.keep_columns(keepcols)
    #specf.keep_columns(['ZWARN','LOCATION','TILEID','TILELOCID','FIBERSTATUS','FIBERASSIGN_X','FIBERASSIGN_Y','PRIORITY','DELTA_X','DELTA_Y','EXPTIME','PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B','TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z','TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
    fgu = join(fgu,specf,keys=['LOCATION','TILEID','FIBER'])
    fgu.sort('TARGETID')
    outf = sv3dir+'/rancomb_'+str(rann)+tp+'wdupspec_Alltiles.fits'
    print(outf)
    if os.path.isfile('tmp.fits'):
        os.system('rm tmp.fits')
    fd = fitsio.FITS('tmp.fits', "rw")
    fd.write(np.array(fgu),extname='POTENTIAL_ASSINGMENT')
    fd['POTENTIAL_ASSINGMENT'].write_comment("concatenation of SV3 POTENTIAL_ASSIGNMENT information for randoms, joined to columns in targe files")
    fd['POTENTIAL_ASSINGMENT'].write_history("updated on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
    fd.close()    
    os.system('mv tmp.fits '+fout)

    #fgu.write(outf,format='fits', overwrite=True)
    


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

def mkfullran(gtl,lznp,indir,rann,imbits,outf,tp,pd,notqso='',maxp=103400,min_tsnr2=0,tlid_full=None,badfib=None):

    '''
    indir is directory with inputs
    rann is the random file number (0-17)
    imbits are the maskbits for the imaging veto mask
    outf is the name (including full path) of the output file
    tp is the target type
    pd is the program, dark or bright
    bit is the bit to use to select to the target type
    randir doesn't get used anymore
    desitarg is the column to use to select the target type
    tsnr is the tsnr2 used for this sample
    '''
    
    #first, need to find locations to veto based on data
    #the same is done in mkfulldat
    #fs = fitsio.read(indir+'datcomb_'+pd+'_specwdup_Alltiles.fits')
#     wf = fs[fbcol] == 0
#     stlid = 10000*fs['TILEID'] +fs['LOCATION']
#     gtl = np.unique(stlid[wf])
#     #gtl now contains the list of good locations
#     #we now want to load in the bigger data file with all the target info
#     #we use it to find the locations where observations of the given type were not possible and then mask them
#     zf = indir+'datcomb_'+pd+'_tarspecwdup_Alltiles.fits'
#     dz = Table.read(zf) 
#     wtype = ((dz[desitarg] & bit) > 0)
#     if notqso == 'notqso':
#         wtype &= ((dz[desitarg] & qsobit) == 0)
# 
#     
#     dz = dz[wtype]#&wg]

    #print('length after selecting type and fiberstatus == 0 '+str(len(dz)))
    #lznp = common.find_znotposs(dz)
    #lznp,lfull = common.find_znotposs_tloc(dz)
    
    #aloc = dz['ZWARN'] != 999999
    #aloc &= dz['ZWARN']*0 != 0
    #alocid = np.unique(dz[aloc]['TILELOCID'])

    #lznp will later be used to veto
    #load in random file
    if pd == 'bright':
        tscol = 'TSNR2_BGS'
    else:
        tscol = 'TSNR2_ELG'


    print('in mkfullran')
    zf = indir+'/rancomb_'+str(rann)+pd+'wdupspec_Alltiles.fits'
    dz = Table(fitsio.read(zf))
    print('loaded wdup file')
    wg = np.isin(dz['TILELOCID'],gtl)
    if badfib is not None:
        bad = np.isin(dz['FIBER'],badfib)
        print('number at bad fibers '+str(sum(bad)))
        wg &= ~bad

    dz['GOODHARDLOC'] = np.zeros(len(dz)).astype('bool')
    dz['GOODHARDLOC'][wg] = 1

    #wa = np.isin(dz['TILELOCID'],alocid)
    #dz['LOC_ASSIGNED'] = np.zeros(len(dz)).astype('bool')
    #dz['LOC_ASSIGNED'][wa] = 1


    wk = ~np.isin(dz['TILELOCID'],lznp)
    dz['ZPOSSLOC'] = np.zeros(len(dz)).astype('bool')
    dz['ZPOSSLOC'][wk] = 1

    dz['LOCFULL'] = np.zeros(len(dz)).astype('bool')
    if tlid_full is not None:
        wf = np.isin(dz['TILELOCID'],tlid_full)
    
        dz['LOCFULL'][wf] = 1


    #load in tileloc info for this random file and join it
    zfpd = indir+'/rancomb_'+str(rann)+pd+'_Alltilelocinfo.fits'
    dzpd = fitsio.read(zfpd)
    dz = join(dz,dzpd,keys=['TARGETID'])
    print('length before cutting to good positions '+str(len(dz)))
    #cut to good and possible locations
    #wk = ~np.isin(dz['TILELOCID'],lznp)
    #wk &= np.isin(dz['TILELOCID'],gtl)
    #dz = dz[wk]    
    print('length after cutting to good positions '+str(len(dz)))
    #get all the additional columns desired from original random files through join
    tarf = Table(fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random'+str(rann)+'/alltilesnofa.fits'))
    delcols = ['RA','DEC','DESI_TARGET','BGS_TARGET','MWS_TARGET','SUBPRIORITY','OBSCONDITIONS','PRIORITY_INIT','NUMOBS_INIT','SCND_TARGET',\
    'NUMOBS_MORE','NUMOBS','Z','ZWARN','TARGET_STATE','TIMESTAMP','VERSION','PRIORITY']
    tarf.remove_columns(delcols)
    dz = join(dz,tarf,keys=['TARGETID'])
    
    #apply imaging vetos
    dz = common.cutphotmask(dz,imbits)
    print('length after cutting to based on imaging veto mask '+str(len(dz)))
    #pl = np.copy(dz['PRIORITY']).astype(float)#dz['PRIORITY']
    #sp = pl <= 0
    #pl[sp] = .1
    dz['GOODPRI'] = np.zeros(len(dz)).astype('bool')
    sel = dz['PRIORITY'] <= maxp
    dz['GOODPRI'][sel] = 1
    dz = np.array(dz)
    np.random.shuffle(dz)
    dz = Table(dz)
    t0 = dz[tscol]*0 != 0
    t0 |= dz[tscol] == 999999
    t0 |= dz[tscol] == 1.e20
    dz[tscol][t0] = 0
    dz['GOODTSNR'] = np.zeros(len(dz)).astype('bool')
    sel = dz[tscol] > min_tsnr2
    dz['GOODTSNR'][sel] = 1
    dz['sort'] =  dz['GOODPRI']*dz['GOODHARDLOC']*dz['ZPOSSLOC']*dz['GOODTSNR']*1-0.5*dz['LOCFULL']#+0.5*dz['LOC_ASSIGNED']#*(1+dz[tsnr])
    #dz[tsnr]*dz['GOODHARDLOC']*dz['ZPOSSLOC']+dz['GOODHARDLOC']*dz['ZPOSSLOC']+dz['GOODHARDLOC']*dz['ZPOSSLOC']/pl
    #sort by tsnr, like done for data, so that the highest tsnr are kept
    dz.sort('sort') 
    dz = unique(dz,keys=['TARGETID'],keep='last')
    print('length after cutting to unique TARGETID '+str(len(dz)))
    sel = dz['GOODPRI']
    sel &= dz['ZPOSSLOC']
    sel &= dz['GOODTSNR']
    sel &= dz['GOODHARDLOC']
    print('number passing vetos (not cut) '+str(len(dz[sel])))
    #tids,cts = np.unique(dz['TILEID'],return_counts=True)
    #plt.plot(tids,cts/np.sum(cts))
    #plt.xlabel('TILEID')
    #plt.ylabel('fraction of randoms')
    #plt.show()
    dz['rosette_number'] = 0
    dz['rosette_r'] = np.zeros(len(dz))
    for ii in range(0,len(dz)):
        rosn = tile2rosette(dz[ii]['TILEID'])
        rosd = calc_rosr(rosn,dz[ii]['RA'],dz[ii]['DEC']) #calculates distance in degrees from the rosette center
        dz[ii]['rosette_number'] = rosn
        dz[ii]['rosette_r'] = rosd
    print(np.unique(dz['NTILE']))
    if int(rann) < 10:
        cof = fitsio.read(outf[:-23]+'_comp_tile.fits')
    else:
        cof = fitsio.read(outf[:-24]+'_comp_tile.fits')  
    comp_dicta = dict(zip(cof['TILES'], cof['COMP_TILE']))
    fcompa = []
    tls = dz['TILES']
    ctls = cof['TILES']
    ctiles = np.zeros(len(dz))
    tlsd = np.isin(tls,cof['TILES'])
    print('number of tiles groups in randoms not in data '+str(len(np.unique(tls[~tlsd]))))
    for i in range(0,len(tls)):
        if tlsd[i]:#np.isin(tl,ctls):
            ctiles[i] = comp_dicta[tls[i]]
        #    fcompa.append(comp_dicta[tl]) 
        #else:
        #    fcompa.append(0)
    dz['COMP_TILE'] = ctiles#np.array(fcompa)
    wc0 = dz['COMP_TILE'] == 0
    print('number of randoms in 0 completeness regions '+str(len(dz[wc0])))   
    
#     cof = fitsio.read(outf[:-23]+'_comp_tileloc.fits')
#     pd = dict(zip(cof['TILELOCID'],cof['FRACZ_TILELOCID']))
#     probl = np.zeros(len(dz))
#     no = 0
#     tidsd = np.isin(dz['TILELOCID'],cof['TILELOCID'])
#     for i in range(0,len(dz)):
#         if tidsd[i]:#np.isin(dz['TILELOCID'][i],cof['TILELOCID']):
#             probl[i] = pd[dz['TILELOCID'][i]]
#         else:
#             no += 1
#     print('number of tilelocid in randoms not in data '+str(no))    
#     dz['FRACZ_TILELOCID'] = probl
    
    comments = ["SV3 'full' LSS catalog for random # "+str(rann)+" without any vetos applied","entries are for all targetid that showed up in POTENTIAL_ASSIGNMENTS"]
    common.write_LSS(dz,outf,comments)

#     tmpfn = outf+'.tmp'
#     if os.path.isfile(tmpfn):
#         os.system('rm '+tmpfn)
#     fd = fitsio.FITS(tmpfn, "rw")
#     fd.write(np.array(dz),extname='LSS')
#     fd['LSS'].write_comment("SV3 'full' LSS catalog for random # "+str(rann)+" without any vetos applied")
#     fd['LSS'].write_comment("entries are for all targetid that showed up in POTENTIAL_ASSIGNMENTS")
#     fd['LSS'].write_history("updated on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
#     fd.close()    
#     os.system('mv '+tmpfn+' '+outf)

    #dz.write(outf,format='fits', overwrite=True)
    #print('moved output to '+outf)
    


def mkfulldat(zf,imbits,tdir,tp,bit,outf,ftiles,azf='',azfm='cumul',desitarg='SV3_DESI_TARGET',specver='fuji',notqso='',qsobit=4,bitweightfile=None,min_tsnr2=0,badfib=None):
    '''
    zf is the name of the file containing all of the combined spec and target info compiled already
    imbits is the list of imaging mask bits to mask out
    tdir is the directory for the targets
    tp is the target type
    bit is the SV3_{type}_MASK bit to use for select the correct target type
    outf is the full path + name for the output file
    ftiles is the name of the file containing information on, e.g., how many tiles each target was available on
    azf is the file name for OII flux info (relevant for ELGs only)
    desitarg is the column to use for the target type cut (all use SV3_DESI_TARGET except BGS_BRIGHT)
    specver is the version of the pipeline used for the redshift info; only 'daily' exists for now
    '''
    
    
    #from desitarget.mtl import inflate_ledger
    if tp[:3] == 'BGS' or tp[:3] == 'MWS':
        pd = 'bright'        
        tscol = 'TSNR2_BGS'
    else:    
        pd = 'dark'
        tscol = 'TSNR2_ELG'
    #load in the appropriate dark/bright combined spec file and use to denote the tileid + location that had good observations:
    #fs = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'+specver+'/datcomb_'+pd+'_specwdup_Alltiles.fits')
    if specver == 'daily':
        #fbcol = 'FIBERSTATUS'
        print('no longer supported')
        return False
    #if specver == 'everest' or specver == 'fuji':
    #else:
    #    fbcol = 'COADD_FIBERSTATUS'
    #wf = fs[fbcol] == 0
    #stlid = 10000*fs['TILEID'] +fs['LOCATION']
    #gtl = np.unique(stlid[wf])
    #gtl now contains the list of 'good' tilelocid

    #read in the big combined data file
    dz = Table.read(zf) 
    #find the rows that satisfy the target type
    wtype = ((dz[desitarg] & bit) > 0)
    if notqso == 'notqso':
        print('removing QSO targets')
        wtype &= ((dz[desitarg] & qsobit) == 0)
    #find the rows that are 'good' tilelocid
    
    print(len(dz[wtype]))
    print('length after selecting type '+str(len(dz)))
    #print(len(dz[wg]))
    #down-select to target type of interest and good tilelocid
    dz = dz[wtype]#&wg]
    
    wz = dz['ZWARN'] != 999999 #this is what the null column becomes
    wz &= dz['ZWARN']*0 == 0 #just in case of nans
    wz &= dz['COADD_FIBERSTATUS'] == 0
    
    if badfib is not None:
        bad = np.isin(dz['FIBER'],badfib)
        print('number at bad fibers '+str(sum(bad)))
        wz &= ~bad
    fs = dz[wz]

    print('number of good obs '+str(len(fs)))
    #fs = common.cut_specdat(dz)
    gtl = np.unique(fs['TILELOCID'])
    wg = np.isin(dz['TILELOCID'],gtl)
    dz['GOODHARDLOC'] = np.zeros(len(dz)).astype('bool')
    dz['GOODHARDLOC'][wg] = 1

    #print('length after selecting type and fiberstatus == 0 '+str(len(dz)))
    #print('length of unique targetid after selecting type and fiberstatus == 0 '+str(len(np.unique(dz['TARGETID']))))
    
    #find targets that were never available at the same location as a target of the same type that got assigned to a good location
    #those that were never available are assumed to have 0 probability of assignment so we want to veto this location
    #lznp = find_znotposs(dz)
    #wk = ~np.isin(dz['TILELOCID'],lznp)#dz['ZPOSS'] == 1
    #dz = dz[wk] #0 probability locations now vetoed
    #print('length after priority veto '+str(len(dz)))
    print('joining to full imaging')
    ftar = Table.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/'+pd+'_targets.fits')
    ftar.keep_columns(['TARGETID','EBV','FLUX_G','FLUX_R','FLUX_Z','FLUX_IVAR_G','FLUX_IVAR_R','FLUX_IVAR_Z','MW_TRANSMISSION_G','MW_TRANSMISSION_R',\
            'MW_TRANSMISSION_Z','FRACFLUX_G','FRACFLUX_R','FRACFLUX_Z','FRACMASKED_G','FRACMASKED_R','FRACMASKED_Z','FRACIN_G','FRACIN_R',\
            'FRACIN_Z','NOBS_G','NOBS_R','NOBS_Z','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','FLUX_W1',\
            'FLUX_W2','FLUX_IVAR_W1','FLUX_IVAR_W2','MW_TRANSMISSION_W1','MW_TRANSMISSION_W2','ALLMASK_G','ALLMASK_R','ALLMASK_Z','FIBERFLUX_G',\
            'FIBERFLUX_R','FIBERFLUX_Z','FIBERTOTFLUX_G','FIBERTOTFLUX_R','FIBERTOTFLUX_Z','WISEMASK_W1','WISEMASK_W2','MASKBITS',\
            'RELEASE','BRICKID','BRICKNAME','BRICK_OBJID','MORPHTYPE','PHOTSYS','SHAPE_R'])
    dz = join(dz,ftar,keys=['TARGETID'])
    print('length after join to full targets (should be same) '+str(len(dz)))
    
    #apply imaging veto mask
    dz = common.cutphotmask(dz,imbits)
    
    #load in file with information about where repeats occurred and join it
    dtl = Table.read(ftiles)
    dtl.keep_columns(['TARGETID','NTILE','TILES','TILELOCIDS'])
    dz = join(dz,dtl,keys='TARGETID')
    
    #find the rows where we have spectroscopic observations
    wz = dz['ZWARN'] != 999999 #this is what the null column becomes
    wz &= dz['ZWARN']*0 == 0 #just in case of nans
    
    
    #mark them as having LOCATION_ASSIGNED
    dz['LOCATION_ASSIGNED'] = np.zeros(len(dz)).astype('bool')
    dz['LOCATION_ASSIGNED'][wz] = 1
    #find the TILELOCID that were assigned and mark them as so
    tlids = np.unique(dz['TILELOCID'][wz])
    #test that all with goodhardloc have z
    gin = np.isin(gtl,tlids)
    print('gtl in tlids, should be all',np.sum(gin),len(gtl))
    wtl = np.isin(dz['TILELOCID'],tlids)
    dz['TILELOCID_ASSIGNED'] = 0
    dz['TILELOCID_ASSIGNED'][wtl] = 1
    print('number of unique targets at assigned tilelocid:')
    print(len(np.unique(dz[wtl]['TARGETID'])))


    
    #sort and then cut to unique targetid; sort prioritizes observed targets and then TSNR2
    wnts = dz[tscol]*0 != 0
    wnts |= dz[tscol] == 999999
    dz[tscol][wnts] = 0
    print(np.max(dz[tscol]))
    dz['GOODTSNR'] = np.zeros(len(dz)).astype('bool')
    sel = dz[tscol] > min_tsnr2
    dz['GOODTSNR'][sel] = 1

    #dz['sort'] = dz['LOCATION_ASSIGNED']*np.clip(dz[tscol],0,200)*dz['GOODHARDLOC']+dz['TILELOCID_ASSIGNED']*dz['GOODHARDLOC']+dz['GOODHARDLOC']
    dz['sort'] = dz['LOCATION_ASSIGNED']*dz['GOODTSNR']*dz['GOODHARDLOC']*(1+np.clip(dz[tscol],0,200))+dz['TILELOCID_ASSIGNED']*dz['GOODHARDLOC']*1+dz['GOODHARDLOC']*1

    print('sort min/max',np.min(dz['sort']),np.max(dz['sort']))
    

    #get OII flux info for ELGs
    if tp[:3] == 'ELG' and azfm == 'cumul':
        if azf != '':
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
   
            dz.remove_columns(['SUBSET','DELTACHI2_OII'])
            print('check length after merge with OII strength file:' +str(len(dz)))
            #join changes order, so get wz again
            wz = dz['ZWARN'] != 999999 #this is what the null column becomes
            wz &= dz['ZWARN']*0 == 0 #just in case of nans

    if tp[:3] == 'QSO' and azfm == 'cumul':
        if azf != '':
            arz = Table(fitsio.read(azf))
            arz.keep_columns(['TARGETID','LOCATION','TILEID','Z','ZERR','Z_QN'])
            arz['TILEID'] = arz['TILEID'].astype(int)
            print(arz.dtype.names)
            #arz['TILE'].name = 'TILEID'
            dz = join(dz,arz,keys=['TARGETID','TILEID','LOCATION'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
            dz['Z'].name = 'Z_RR' #rename the original redrock redshifts
            dz['Z_QF'].name = 'Z' #the redshifts from the quasar file should be used instead
            #join changes order, so get wz again
            wz = dz['ZWARN'] != 999999 #this is what the null column becomes
            wz &= dz['ZWARN']*0 == 0 #just in case of nans

    dz.sort('sort')
    dz = unique(dz,keys=['TARGETID'],keep='last')

    if azfm == 'hp':
        arz = fitsio.read('/global/cfs/cdirs/desi/spectro/redux/'+specver+'/zcatalog/zpix-sv3-'+pd+'.fits',columns=['TARGETID','Z','DELTACHI2','TSNR2_ELG','TSNR2_LRG','TSNR2_QSO','TSNR2_BGS']  )
        dz = join(dz,arz,keys=['TARGETID'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['', '_HP'])


    if tp[:3] == 'ELG' and azfm == 'hp':
        if azf != '':
            arz = fitsio.read(azf,columns=['TARGETID','OII_FLUX','OII_FLUX_IVAR','DELTACHI2'])
            o2c = np.log10(arz['OII_FLUX'] * np.sqrt(arz['OII_FLUX_IVAR']))+0.2*np.log10(arz['DELTACHI2'])
            w = (o2c*0) != 0
            w |= arz['OII_FLUX'] < 0
            o2c[w] = -20
            arz = Table(arz)
            arz['o2c'] = o2c
            dz = join(dz,arz,keys=['TARGETID'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['', '_OII'])
   
            dz.remove_columns(['DELTACHI2_OII'])
            print('check length after merge with OII strength file:' +str(len(dz)))

    if tp[:3] == 'QSO' and azfm == 'hp':
        if azf != '':
            arz = Table(fitsio.read(azf))
            arz.keep_columns(['TARGETID','Z','ZERR','Z_QN'])
            print(arz.dtype.names)
            dz = join(dz,arz,keys=['TARGETID'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
            dz['Z_HP'].name = 'Z_RR' #rename the original redrock redshifts
            dz['Z_QF'].name = 'Z_HP' #the redshifts from the quasar file should be used instead

    if tp == 'ELG' or tp == 'ELG_HIP':
        print('number of masked oII row (hopefully matches number not assigned) '+ str(np.sum(dz['o2c'].mask)))
    if tp == 'QSO':
        if azfm == 'hp':
            print('number of good z according to qso file '+str(len(dz)-np.sum(dz['Z_HP'].mask)))
            dz['Z_HP'] = dz['Z_HP'].filled(999999)
        else:
            print('number of good z according to qso file '+str(len(dz)-np.sum(dz['Z'].mask)))
    dz['Z'] = dz['Z'].filled(999999)
    selm = dz['Z'] == 999999
    print('999999s for Z',len(dz[selm]))
    print('length after cutting to unique targetid '+str(len(dz)))
    print('LOCATION_ASSIGNED numbers')
    print(np.unique(dz['LOCATION_ASSIGNED'],return_counts=True))
   
    print('TILELOCID_ASSIGNED numbers')
    print(np.unique(dz['TILELOCID_ASSIGNED'],return_counts=True))

    
    #get completeness based on unique sets of tiles "comp_tile"
    tll,compa = common.comp_tile(dz)
    comp_dicta = dict(zip(tll, compa))
    fcompa = []
    for tl in dz['TILES']:
        fcompa.append(comp_dicta[tl]) 
    dz['COMP_TILE'] = np.array(fcompa)
    wc0 = dz['COMP_TILE'] == 0
    print('number of targets in 0 completeness regions '+str(len(dz[wc0])))   
    
    #write out comp_tile info
    tll = np.array(tll).astype(dz['TILES'].dtype)
    co = Table()
    co['TILES'] = tll
    co['COMP_TILE'] = compa
    cof = outf.replace('_full_noveto.dat.fits','_comp_tile.fits')
    print('writing comp_tile completeness to '+cof)
    co.write(cof,overwrite=True,format='fits')    

    #get counts at unique TILELOCID
    locl,nlocl = np.unique(dz['TILELOCID'],return_counts=True)
    #do same after cutting to only the data with location_assigned
    wz = dz['LOCATION_ASSIGNED'] == 1
    dzz = dz[wz]
    loclz,nloclz = np.unique(dzz['TILELOCID'],return_counts=True)
    natloc = ~np.isin(dz['TILELOCID'],loclz)
    print('number of unique targets left around unassigned locations is '+str(np.sum(natloc)))
#    locs = np.copy(dz['TILELOCID'])
# 
# 
    print('reassigning TILELOCID for duplicates and finding rosette')
    #re-assigning "naked" targets; if we gave a targetid a tilelocid that was not assigned 
    #by the same target was available at a location that was assigned, we re-assign its tilelocid
    nch = 0
    nbl = 0
    tlids = dz['TILELOCIDS']
    ros = np.zeros(len(dz))
    rosr = np.zeros(len(dz))
    for ii in range(0,len(dz['TILEID'])): #not sure why, but this only works when using loop for Table.read but array option works for fitsio.read
        ti = dz[ii]['TILEID']
        rosn = tile2rosette(ti) #get rosette id
        rosr[ii] = calc_rosr(rosn,dz[ii]['RA'],dz[ii]['DEC']) #calculates distance in degrees from rosette center
        ros[ii] = rosn
#         if natloc[ii]:# == False:
#             nbl += 1
#             s = 0
#             tids = tlids[ii].split('-')
#             if s == 0:
#                 for tl in tids:
#                     ttlocid  = int(tl)              
#                     if np.isin(ttlocid,loclz):
#                         locs[ii] = ttlocid 
#                         nch += 1
#                         s = 1
#                         break
#         if ii%10000 == 0:
#             print(ii,len(dz['TILEID']),ti,ros[ii],nch,nbl)
     

    dz['rosette_number'] = ros
    dz['rosette_r'] = rosr
    print('rosette number and the number on each rosette')
    print(np.unique(dz['rosette_number'],return_counts=True))

#    dz['TILELOCID'] = locs
    #get numbers again after the re-assignment
#     locl,nlocl = np.unique(dz['TILELOCID'],return_counts=True)
#     loclz,nloclz = np.unique(dzz['TILELOCID'],return_counts=True)
# 
#     print('getting fraction assigned for each tilelocid')
#     #should be one (sometimes zero, though) assigned target at each tilelocid and we are now counting how many targets there are per tilelocid
#     #probability of assignment is then estimated as 1/n_tilelocid
#     nm = 0
#     nmt =0
#     pd = []
#     for i in range(0,len(locl)):
#         if i%10000 == 0:
#             print('at row '+str(i))
#         nt = nlocl[i]
#         loc = locl[i]
#         w = loclz == loc
#         nz = 0
#         if len(loclz[w]) == 1:
#             nz = nloclz[w] #these are supposed all be 1...            
#         else:            
#             nm += 1.
#             nmt += nt
#         if len(loclz[w]) > 1:
#             print('why is len(loclz[w]) > 1?') #this should never happen
#         pd.append((loc,nz/nt))  
    loco,fzo = common.comp_tileloc(dz)
    pd = dict(zip(loco,fzo))
    probl = np.zeros(len(dz))
    for i in range(0,len(dz)):
        probl[i] = pd[dz['TILELOCID'][i]]
    dz['FRACZ_TILELOCID'] = probl

    #write out FRACZ_TILELOCID info
    #loco = np.array(loco).astype(dz['TILELOCID'].dtype)
    co = Table()
    co['TILELOCID'] = loco
    co['FRACZ_TILELOCID'] = fzo
    cof = outf.strip('_full_noveto.dat.fits')+'_comp_tileloc.fits'
    print('writing comp_tileloc completeness to '+cof)
    co.write(cof,overwrite=True,format='fits')    


    print('sum of 1/FRACZ_TILELOCID, 1/COMP_TILE, length of input, number of inputs with obs; no longer rejecting unobserved loc, so wont match')
    print(np.sum(1./dz[wz]['FRACZ_TILELOCID']),np.sum(1./dz[wz]['COMP_TILE']),len(dz),len(dz[wz]))
    
    oct = np.copy(dz['COMP_TILE'])
    if bitweightfile is not None:
        fb = fitsio.read(bitweightfile)
        dz = join(dz,fb,keys=['TARGETID'])
    wz = dz['LOCATION_ASSIGNED'] == 1 #join re-ordered array, reget mask for assigned locations and check comp_tile
    print('length after join with bitweight file and sum of 1/comp_tile',len(dz),np.sum(1./dz[wz]['COMP_TILE']),len(dz[wz]))
    #print('check comp_tile array',np.array_equal(oct,dz['COMP_TILE']))

    
    #for debugging writeout
#     for col in dz.dtype.names:
#         to = Table()
#         to[col] = dz[col]
#         #print(col)
#         try:
#             to.write('temp.fits',format='fits', overwrite=True)
#         except:
#             print(col+' failed!')            

    comments = ["SV3 'full' LSS catalog for data without any vetos applied","entries are for all targetid that showed up in POTENTIAL_ASSIGNMENTS"]
    common.write_LSS(dz,outf,comments)

#     tmpfn = outf +'.tmp'
#     if os.path.isfile(tmpfn):
#         os.system('rm '+tmpfn)
#     fd = fitsio.FITS(tmpfn, "rw")
#     fd.write(np.array(dz),extname='LSS')
#     fd['LSS'].write_comment("'full' LSS catalog for data without any vetos applied")
#     fd['LSS'].write_comment("entries are for all targetid that showed up in POTENTIAL_ASSIGNMENTS")
#     fd['LSS'].write_history("updated on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
#     fd.close()    
#     os.system('mv '+tmpfn+' '+outf)
    
    #dz.write(outf,format='fits', overwrite=True)

def mkclusdat(fl,weightmd='tileloc',zmask=False,tp='',dchi2=9,tsnrcut=80,rcut=None,ntilecut=0,ccut=None,ebits=None,nreal=128,zmin=0.01,zmax=6,hp='hp'):
    '''
    fl is the root of the input/output file
    weighttileloc determines whether to include 1/FRACZ_TILELOCID as a completeness weight
    zmask determines whether to apply a mask at some given redshift
    tp is the target type
    dchi2 is the threshold for keeping as a good redshift
    tnsrcut determines where to mask based on the tsnr2 value (defined below per tracer)

    '''    
    ff = Table.read(fl+'full.dat.fits')
#     if ebits is not None:
#         print('number before imaging mask '+str(len(ff)))
#         if ebits == 'lrg_mask':
#             sel = ff['lrg_mask'] == 0
#             ff = ff[sel]
#         else:
#             ff = cutphotmask(ff,ebits)
#         print('number after imaging mask '+str(len(ff)))
#    ff.write(fl+'full.dat.fits',overwrite=True,format='fits')
    wzm = ''
    if zmask:
        wzm = 'zmask_'
    if rcut is not None:
        wzm += 'rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
    if ntilecut > 0:
        wzm += 'ntileg'+str(ntilecut)+'_'    
    if ccut is not None:
        wzm += ccut+'_' #you could change this to however you want the file names to turn out

    if ccut == 'main':
        if tp != 'LRG':
            print('this is only defined for LRGs!' )
        else:
            lrgmaintar = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/LRGtargetsDR9v1.1.1.fits',columns=['TARGETID'])
            sel = np.isin(ff['TARGETID'],lrgmaintar['TARGETID'])
            print('numbers before/after cut:')
            print(len(ff),len(ff[sel]))
            ff = ff[sel]   
            ff.write(fl+wzm+'full.dat.fits',format='fits',overwrite='True')
    zfcol = 'Z_not4clus'
    if hp == 'hp':
        zfcol = 'Z_HP'
    ff[zfcol].name = 'Z'    
    '''
    Not doing systematic weights for SV3, just setting them to 1
    '''

    ff['WEIGHT_ZFAIL'] = np.ones(len(ff))
    #The LRGs just have this fairly ad hoc model that AJR fit in the notebook, definitely needs refinement/automation
#     if tp == 'LRG':
#         fibfluxz = ff['FIBERFLUX_Z']/ff['MW_TRANSMISSION_Z']
#         coeff = [117.46,-60.91,11.49,-0.513] #from polyfit, 3rd to zeroth order in 1/fiberflu
#         efs = coeff[-1]+coeff[-2]*(1/fibfluxz)+coeff[-3]*(1/fibfluxz)**2.+coeff[-4]*(1/fibfluxz)**3.
#         ems = erf((ff['TSNR2_LRG']-13.2)/39.7)*.9855
#         ff['WEIGHT_ZFAIL'] = 1./(1. -(1.-ems)*efs)



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
    
    if tp == 'ELG' or tp == 'ELG_HIP':
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
        #drz = (10**(3 - 3.5*ff['Z']))
        #mask_bad = (drz>30) & (ff['DELTACHI2']<30)
        #mask_bad |= (drz<30) & (ff['DELTACHI2']<drz)
        #mask_bad |= (ff['DELTACHI2']<10)
        #wz &= ff['Z']<1.4
        #wz &= (~mask_bad)
        wz &= ff['ZWARN']*0 == 0
        wz &= ff['ZWARN'] != 999999
        wz &= ff['ZWARN'] != 1.e20

        selg = ssr_tools.LRG_goodz(ff)
        wz &= selg

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
    print('minimum,maximum Z',min(ff['Z']),max(ff['Z']))
    ff['WEIGHT'] = ff['WEIGHT_ZFAIL']
    ff['WEIGHT_COMP'] = np.ones(len(ff))
    if weightmd == 'tileloc':
        ff['WEIGHT_COMP'] = 1./ff['FRACZ_TILELOCID']
    
    if weightmd == 'probobs' :         
        nassign = nreal*ff['PROB_OBS']+1 #assignment in actual observation counts
        ff['WEIGHT_COMP'] *= (nreal+1)/nassign#1./ff['PROB_OBS']
        print(np.min(ff['WEIGHT']),np.max(ff['WEIGHT']))
        #wzer = ff['PROB_OBS'] == 0
        #ff['WEIGHT'][wzer] = 0
        #print(str(len(ff[wzer]))+' galaxies with PROB_OBS 0 getting assigned weight of 0 (should not happen, at minimum adjust weights to reflect 1 real realization happened)')
    
    ff['WEIGHT'] *= ff['WEIGHT_COMP']
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

    selz = ff['Z'] > zmin
    selz &= ff['Z'] < zmax
    ff = ff[selz]


    #select down to specific columns below and then also split N/S
    wn = ff['PHOTSYS'] == 'N'
    #if tp[:3] == 'BGS':
    #    kl = ['RA','DEC','Z','WEIGHT','TARGETID','NTILE','rosette_number','rosette_r','TILES','WEIGHT_ZFAIL','FRACZ_TILELOCID']
    #else:
    kl = ['RA','DEC','Z','WEIGHT','TARGETID','NTILE','COMP_TILE','rosette_number','rosette_r','TILES','WEIGHT_ZFAIL','FRACZ_TILELOCID','PROB_OBS','BITWEIGHTS']    
    if tp[:3] == 'BGS':
        #ff['flux_r_dered'] = ff['FLUX_R']/ff['MW_TRANSMISSION_R']
        fcols = ['G','R','Z','W1','W2']
        ff = common.add_dered_flux(ff,fcols)
        for col in fcols:
            kl.append('flux_'+col.lower()+'_dered')
        print(kl)

    ff.keep_columns(kl)#,'PROB_OBS'
    print('minimum,maximum weight')
    print(np.min(ff['WEIGHT']),np.max(ff['WEIGHT']))

    #comments = ["SV3 'clustering' LSS catalog for data, all regions","entries are only for data with good redshifts"]
    #common.write_LSS(ff,outf,comments)

    outfn = fl+wzm+'N_clustering.dat.fits'
    comments = ["SV3 'clustering' LSS catalog for data, just in BASS/MzLS region","entries are only for data with good redshifts"]
    common.write_LSS(ff[wn],outfn,comments)

    outfn = fl+wzm+'S_clustering.dat.fits'
    comments = ["SV3 'clustering' LSS catalog for data, just in DECaLS region","entries are only for data with good redshifts"]
    common.write_LSS(ff[~wn],outfn,comments)

#     tmpfn = outf+'.tmp'
#     if os.path.isfile(tmpfn):
#         os.system('rm '+tmpfn)
#     fd = fitsio.FITS(tmpfn, "rw")
#     fd.write(np.array(ff),extname='LSS')
#     fd['LSS'].write_comment("'cluster' LSS catalog for data, all regions")
#     fd['LSS'].write_comment("entries are only for data with good redshifts")
#     fd['LSS'].write_history("updated on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
#     fd.close()    
#     os.system('mv '+tmpfn+' '+outf)

    #ff.write(outf,format='fits', overwrite=True)
#     outfn = fl+wzm+'N_clustering.dat.fits'
#     tmpfn = outfn+'.tmp'
#     if os.path.isfile(tmpfn):
#         os.system('rm '+tmpfn)
#     fd = fitsio.FITS(tmpfn, "rw")
#     fd.write(np.array(ff[wn]),extname='LSS')
#     fd['LSS'].write_comment("'cluster' LSS catalog for data, just in BASS/MzLS region")
#     fd['LSS'].write_comment("entries are only for data with good redshifts")
#     fd['LSS'].write_history("updated on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
#     fd.close()    
#     os.system('mv '+tmpfn+' '+outf)
    #ff[wn].write(outfn,format='fits', overwrite=True)
#     outfn = fl+wzm+'S_clustering.dat.fits'
#     tmpfn = outfn+'.tmp'
#     if os.path.isfile(tmpfn):
#         os.system('rm '+tmpfn)
#     fd = fitsio.FITS(tmpfn, "rw")
#     fd.write(np.array(ff[~wn]),extname='LSS')
#     fd['LSS'].write_comment("'cluster' LSS catalog for data, just in DECaLS region")
#     fd['LSS'].write_comment("entries are only for data with good redshifts")
#     fd['LSS'].write_history("updated on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
#     fd.close()    
#     os.system('mv '+tmpfn+' '+outfn)
    #ff[~wn].write(outfn,format='fits', overwrite=True)

def mkclusran(fl,rann,reg='_S',rcols=['Z','WEIGHT'],zmask=False,tsnrcut=80,tsnrcol='TSNR2_ELG',rcut=None,ntilecut=0,ccut=None,ebits=None):
    '''
    fl is the root of our catalog file names
    rann is the random number
    rcols are the columns that we randomly select from the data file
    zmask is whether or not we mask out certain redshift
    tsnrcut is the tsnr2 value below which we discard data
    tsnrcol is the specific column used for the tsnrcut
    '''
    
    wzm = ''
    if zmask:
        wzm = 'zmask_'
    if rcut is not None:
        wzm += 'rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
    if ntilecut > 0:
        wzm += 'ntileg'+str(ntilecut)+'_'    
    if ccut is not None:
        wzm += '_'+ccut  #you could change this to however you want the file names to turn out

    #load in data clustering catalog
    fcd = Table.read(fl+wzm+reg+'_clustering.dat.fits')
    #load in full random file
    ffr = Table.read(fl+'_'+str(rann)+'_full.ran.fits')
    selr = ffr['PHOTSYS'] == reg.strip('_')
    ffr = ffr[selr]
#     if ebits is not None:
#         #print(ebits)
#         print('number before imaging mask '+str(len(ffr)))
#         if ebits == 'lrg_mask':
#             sel = ffr['lrg_mask'] == 0
#             ffr = ffr[sel]
#         else:    
#             ffr = cutphotmask(ffr,ebits)
#         print('number after imaging mask '+str(len(ffr)))
# 
#     ffr.write(fl+str(rann)+'_full.ran.fits',overwrite=True,format='fits')

    #mask mask on tsnr
    wz = ffr[tsnrcol] > tsnrcut
    ffc = ffr[wz]
    print('length after,before tsnr cut:')
    print(len(ffc),len(ffr))
    #apply any cut on rosette radius
    if rcut is not None:
        wr = ffc['rosette_r'] > rcut[0]
        wr &= ffc['rosette_r'] <  rcut[1]
        print('length before rosette radius cut '+str(len(ffc)))
        ffc = ffc[wr]
        print('length after rosette radius cut '+str(len(ffc)))
    #apply cut on ntile
    if ntilecut > 0:
        print('length before ntile cut '+str(len(ffc)))
        wt = ffc['NTILE'] > ntilecut
        ffc = ffc[wt]
        print('length after ntile cut '+str(len(ffc)))    


    #randomly sample data rows to apply redshifts, weights, etc. to randoms
    inds = np.random.choice(len(fcd),len(ffc))
    dshuf = fcd[inds]
    kl =  ['RA','DEC','TARGETID','NTILE','COMP_TILE','rosette_number','rosette_r','TILES'] + rcols 

    for col in rcols: 
        ffc[col] = dshuf[col] 
    #cut to desired small set of columns and write out files, splitting N/S as well
    #wn = ffc['PHOTSYS'] == 'N'
    
    ffc.keep_columns(kl)  
    outf =  fl+wzm+reg+'_'+str(rann)+'_clustering.ran.fits' 

    comments = ["SV3 'clustering' LSS catalog for random #"+str(rann)+reg+" regions","columns that are not ra,dec are sampled from data with good redshifts"]
    common.write_LSS(ffc,outf,comments)


#     tmpfn = outf+'.tmp'
#     if os.path.isfile(tmpfn):
#         os.system('rm '+tmpfn)
#     fd = fitsio.FITS(tmpfn, "rw")
#     fd.write(np.array(ffc),extname='LSS')
#     fd['LSS'].write_comment("'cluster' LSS catalog for random #"+str(rann)+", all regions")
#     fd['LSS'].write_comment("columns that are not ra,dec are sampled from data with good redshifts")
#     fd['LSS'].write_history("updated on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
#     fd.close()    
#     os.system('mv '+tmpfn+' '+outf)
    #ffc.write(outf,format='fits', overwrite=True)
#     outfn =  fl+wzm+'N_'+str(rann)+'_clustering.ran.fits' 
#     fcdn = Table.read(fl+wzm+'N_clustering.dat.fits')
#     ffcn = ffc[wn]
#     inds = np.random.choice(len(fcdn),len(ffcn))
#     dshuf = fcdn[inds]
#     for col in rcols: 
#         ffcn[col] = dshuf[col]     
# 
#     comments = ["SV3 'clustering' LSS catalog for random #"+str(rann)+", BASS/MzLS region","columns that are not ra,dec are sampled from data with good redshifts"]
#     common.write_LSS(ffcn,outfn,comments)


#     tmpfn = outfn+'.tmp'
#     if os.path.isfile(tmpfn):
#         os.system('rm '+tmpfn)
#     fd = fitsio.FITS(tmpfn, "rw")
#     fd.write(np.array(ffcn),extname='LSS')
#     fd['LSS'].write_comment("'cluster' LSS catalog for random #"+str(rann)+", BASS/MzLS region")
#     fd['LSS'].write_comment("columns that are not ra,dec are sampled from data with good redshifts")
#     fd['LSS'].write_history("updated on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
#     fd.close()    
#     os.system('mv '+tmpfn+' '+outfn)
    #ffcn.write(outfn,format='fits', overwrite=True)

#     outfs =  fl+wzm+'S_'+str(rann)+'_clustering.ran.fits' 
#     fcds = Table.read(fl+wzm+'S_clustering.dat.fits')
#     ffcs = ffc[~wn]
#     inds = np.random.choice(len(fcds),len(ffcs))
#     dshuf = fcds[inds]
#     for col in rcols: 
#         ffcs[col] = dshuf[col]     
# 
#     comments = ["SV3 'clustering' LSS catalog for random #"+str(rann)+", DECaLS region","columns that are not ra,dec are sampled from data with good redshifts"]
#     common.write_LSS(ffcs,outfs,comments)


#     tmpfn = outfs+'.tmp'
#     if os.path.isfile(tmpfn):
#         os.system('rm '+tmpfn)
#     fd = fitsio.FITS(tmpfn, "rw")
#     fd.write(np.array(ffcs),extname='LSS')
#     fd['LSS'].write_comment("'cluster' LSS catalog for random #"+str(rann)+", DECaLS region")
#     fd['LSS'].write_comment("columns that are not ra,dec are sampled from data with good redshifts")
#     fd['LSS'].write_history("updated on "+datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
#     fd.close()    
#     os.system('mv '+tmpfn+' '+outfs)
    #ffcs.write(outfs,format='fits', overwrite=True)

# def addnbar(fb,nran=18,bs=0.01,zmin=0.01,zmax=1.6):
#     nzd = np.loadtxt(fb+'_nz.dat').transpose()[3] #column with nbar values
#     fn = fb+'_clustering.dat.fits'
#     fd = fitsio.read(fn) #reading in data with fitsio because it is much faster to loop through than table
#     zl = fd['Z']
#     nl = np.zeros(len(zl))
#     for ii in range(0,len(zl)):
#         z = zl[ii]
#         zind = int((z-zmin)/bs)
#         if z > zmin and z < zmax:
#             nl[ii] = nzd[zind]
#     del fd
#     ft = Table.read(fn)
#     ft['NZ'] = nl
#     ft.write(fn,format='fits',overwrite=True)        
#     print('done with data')
#     for rann in range(0,nran):
#         fn = fb+'_'+str(rann)+'_clustering.ran.fits'
#         fd = fitsio.read(fn) #reading in data with fitsio because it is much faster to loop through than table
#         zl = fd['Z']
#         nl = np.zeros(len(zl))
#         for ii in range(0,len(zl)):
#             z = zl[ii]
#             zind = int((z-zmin)/bs)
#             if z > zmin and z < zmax:
#                 nl[ii] = nzd[zind]
#         del fd
#         ft = Table.read(fn)
#         ft['NZ'] = nl
#         ft.write(fn,format='fits',overwrite=True)      
#         print('done with random number '+str(rann))  
#     return True        
#     
# 
# def mknz(fcd,fcr,fout,bs=0.01,zmin=0.01,zmax=1.6,om=0.31519):
#     
#     cd = distance(om,1-om)
#     ranf = fitsio.read(fcr) #should have originally had 2500/deg2 density, so can convert to area
#     area = len(ranf)/2500.
#     print('area is '+str(area))
#     
#     df = fitsio.read(fcd)
#     
#     nbin = int((zmax-zmin)/bs)
#     zhist = np.histogram(df['Z'],bins=nbin,range=(zmin,zmax),weights=df['WEIGHT'])
#     outf = open(fout,'w')
#     outf.write('#area is '+str(area)+'square degrees\n')
#     outf.write('#zmid zlow zhigh n(z) Nbin Vol_bin\n')
#     for i in range(0,nbin):
#         zl = zhist[1][i]
#         zh = zhist[1][i+1]
#         zm = (zh+zl)/2.
#         voli = area/(360.*360./np.pi)*4.*np.pi/3.*(cd.dc(zh)**3.-cd.dc(zl)**3.)
#         nbarz =  zhist[0][i]/voli
#         outf.write(str(zm)+' '+str(zl)+' '+str(zh)+' '+str(nbarz)+' '+str(zhist[0][i])+' '+str(voli)+'\n')
#     outf.close()

    

def randomtiles_allSV3(tiles,dirout='/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random',imin=0,imax=18):
    '''
    tiles should be a table containing the relevant info
    '''
    trad = desimodel.focalplane.get_tile_radius_deg()*1.1 #make 10% greater just in case
    print(trad)
    for ii in range(imin,imax):
        rt = fitsio.read(dirout+str(ii)+'/alltilesnofa.fits')
        #rt = fitsio.read(minisvdir+'random/random_mtl.fits')
        print('loaded random file') 
    
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
                rmtl['SV3_DESI_TARGET'] = np.ones(len(rmtl),dtype=int)*2
                rmtl['NUMOBS_INIT'] = np.zeros(len(rmtl),dtype=int)
                rmtl['NUMOBS_MORE'] = np.ones(len(rmtl),dtype=int)
                rmtl['PRIORITY'] = np.ones(len(rmtl),dtype=int)*3400
                rmtl['OBSCONDITIONS'] = np.ones(len(rmtl),dtype=int)*516#tiles['OBSCONDITIONS'][i]
                rmtl['SUBPRIORITY'] = np.random.random(len(rmtl))
                rmtl.write(fname,format='fits', overwrite=True)
                print('added columns, wrote to '+fname)
