'''
python functions to do various useful date processing/manipulation
'''
import numpy as np
import glob
import os
from random import random

import astropy.io.fits as fits
from astropy.table import Table,join,unique,vstack,setdiff

import fitsio

from matplotlib import pyplot as plt


import healpy as hp

#from LSS.Cosmo import distance
from LSS.imaging import densvar

#from LSS.common_tools import find_znotposs



import logging
logger = logging.getLogger('cattools')
logger.setLevel(level=logging.INFO)

#logging.getLogger("QSO_CAT_UTILS").setLevel(logging.ERROR)



def combtile_qso(tiles,outf='',restart=False,release='guadalupe'):
    s = 0
    n = 0
    nfail = 0
    kl = ['TARGET_RA','TARGET_DEC','DESI_TARGET','TARGETID', 'Z', 'LOCATION',  'TSNR2_LYA', 'TSNR2_QSO', 'DELTA_CHI2_MGII', 'A_MGII', 'SIGMA_MGII', 'B_MGII', 'VAR_A_MGII', 'VAR_SIGMA_MGII', 'VAR_B_MGII', 'Z_RR', 'Z_QN', 'C_LYA', 'C_CIV', 'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha', 'Z_LYA', 'Z_CIV', 'Z_CIII', 'Z_MgII', 'Z_Hbeta', 'Z_Halpha', 'QSO_MASKBITS', 'TILEID']
    #
    if os.path.isfile(outf) and restart == False:
        #specd = Table.read(outf)
        specd = fitsio.read(outf)

        #dt = specd.dtype
        #specd = np.empty(len(specio),dtype=dt)
        #cols = fw.dtype.names
        #for colname in cols:
        #    specd[colname][...] = specio[colname][...]
        #del specio

        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)
    else:
        infl = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/QSO_catalog_'+release+'.fits'
        specd = Table(fitsio.read(infl))
        specd.keep_columns(kl)
        specd = np.array(specd)
        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)
    print('will add QSO info for '+str(len(tiles[tmask]))+' tiles')
    #kl = list(specd.dtype.names)
    for tile,zdate,tdate in zip(tiles[tmask]['TILEID'],tiles[tmask]['ZDATE'],tiles[tmask]['THRUDATE']):
        tdate = str(tdate)
        tspec = combQSOdata(tile,zdate,tdate,cols=kl)
        if tspec:
            #this is stupid but should speed up concatenation
            #tspec.write('temp.fits',format='fits', overwrite=True)
            #tspec = fitsio.read('temp.fits')
            tspec = np.array(tspec)
            #tspec = np.empty(len(tspecio),dtype=dt)

            if s == 0:
                specd = tspec
                s = 1
            else:
                #specd = vstack([specd,tspec],metadata_conflicts='silent')
                #column order got mixed up
                new = np.empty(len(tspec),dtype=specd.dtype)
                cols = specd.dtype.names
                for colname in cols:
                    new[colname][...] = tspec[colname][...]

                #specd = np.hstack((specd,tspec))
                specd = np.hstack((specd,new))
            #specd.sort('TARGETID')
            kp = (specd['TARGETID'] > 0)
            specd = specd[kp]

            n += 1
            print(tile,n,len(tiles[tmask]),len(specd))
        else:
            print(str(tile)+' failed')
            nfail += 1
    print('total number of failures was '+str(nfail))
    if n > 0:
        #specd.write(outf,format='fits', overwrite=True)
        fitsio.write(outf,specd,clobber=True)
        return True
    else:
        return False

def combtile_qso_alt(tiles,outf='',coaddir=''):
    #to be used for rerun
    s = 0
    n = 0
    nfail = 0
    kl = ['TARGET_RA','TARGET_DEC','DESI_TARGET','TARGETID', 'Z', 'LOCATION',  'TSNR2_LYA', 'TSNR2_QSO', 'DELTA_CHI2_MGII', 'A_MGII', 'SIGMA_MGII', 'B_MGII', 'VAR_A_MGII', 'VAR_SIGMA_MGII', 'VAR_B_MGII', 'Z_RR', 'Z_QN', 'C_LYA', 'C_CIV', 'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha', 'Z_LYA', 'Z_CIV', 'Z_CIII', 'Z_MgII', 'Z_Hbeta', 'Z_Halpha', 'QSO_MASKBITS', 'TILEID']
    print('will add QSO info for '+str(len(tiles))+' tiles')
    #kl = list(specd.dtype.names)
    for tile,zdate in zip(tiles['TILEID'],tiles['THRUDATE']):
        zdate = str(zdate)
        tspec = combQSOdata_alt(tile,zdate,cols=kl,coaddir=coaddir)
        if tspec:
            #better would be to concatenate all at once than one at a time
            tspec = np.array(tspec)

            if s == 0:
                specd = tspec
                s = 1
            else:
                new = np.empty(len(tspec),dtype=specd.dtype)
                cols = specd.dtype.names
                for colname in cols:
                    new[colname][...] = tspec[colname][...]

                
                specd = np.hstack((specd,new))
            
            kp = (specd['TARGETID'] > 0)
            specd = specd[kp]

            n += 1
            print(tile,n,len(tiles),len(specd))
        else:
            print(str(tile)+' failed')
            nfail += 1
    print('total number of failures was '+str(nfail))
    if n > 0:
        
        fitsio.write(outf,specd,clobber=True)
        return True
    else:
        return False



def combtile_spec(tiles,outf='',md='',specver='daily',redo='n',specrel='guadalupe',prog='dark',par='n'):
    import LSS.common_tools as common
    logger = logging.getLogger('comb_inputs')
    s = 0
    n = 0
    nfail = 0
    tl = []
    if os.path.isfile(outf) and redo == 'n':
        #specd = Table.read(outf)
        specd = fitsio.read(outf)
        tl.append(specd)
        #dt = specd.dtype
        #specd = np.empty(len(specio),dtype=dt)
        #cols = fw.dtype.names
        #for colname in cols:
        #    specd[colname][...] = specio[colname][...]
        #del specio
        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)

    elif redo == 'y':
        specd = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/'+specrel+'/datcomb_'+prog+'_spec_zdone.fits')
        tl.append(specd)
        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)
    else:
        tmask = np.ones(len(tiles)).astype('bool')

    if par == 'y':
        print('defunct')
    else:
        for tile,zdate,tdate in zip(tiles[tmask]['TILEID'],tiles[tmask]['ZDATE'],tiles[tmask]['THRUDATE']):
            tdate = str(tdate)
            if md == 'zmtl':
                #if specver ==
                #tspec = combzmtl(tile,zdate,tdate)
                tspec = combzmtl(tile,tdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/'+specver+'/tiles/cumulative/')
            else:
                tspec = combspecdata(tile,zdate,tdate)
            if tspec:
                tspec['TILEID'] = tile
                tspec = np.array(tspec)
                #this is stupid but should speed up concatenation
                #tspec.write('temp.fits',format='fits', overwrite=True)
                #tspec = fitsio.read('temp.fits')
                #tspec = np.empty(len(tspecio),dtype=dt)
    
                if s == 0:
                    specd = tspec
                    s = 1
                #else:
                #specd = vstack([specd,tspec],metadata_conflicts='silent')
                #column order got mixed up
                new = np.empty(len(tspec),dtype=specd.dtype)
                cols = specd.dtype.names
                for colname in cols:
                    new[colname][...] = tspec[colname][...]
    
                    #specd = np.hstack((specd,tspec))
                    #specd = np.hstack((specd,new))
                tl.append(new)
                #specd.sort('TARGETID')
                #kp = (specd['TARGETID'] > 0)
                #specd = specd[kp]
    
                n += 1
                common.printlog(str(tile)+','+str(n)+','+str(len(tiles[tmask])),logger)#,len(specd))
            else:
                print(str(tile)+' failed')
                nfail += 1
    specd = np.hstack(tl)
    kp = (specd['TARGETID'] > 0)
    specd = specd[kp]
    print('total number of failures was '+str(nfail))
    if n > 0:
        #specd.write(outf,format='fits', overwrite=True)
        #fitsio.write(outf,specd,clobber=True)
        common.write_LSS_scratchcp(specd,outf,logger=logger)
        
        return True
    else:
        return False

def combtile_spec_alt(tiles,outf='',md='',redo='n',prog='dark',coaddir=''):
    s = 0
    n = 0
    nfail = 0
    tl = []
    if os.path.isfile(outf) and redo == 'n':
        #specd = Table.read(outf)
        specd = fitsio.read(outf)
        tl.append(specd)
        #dt = specd.dtype
        #specd = np.empty(len(specio),dtype=dt)
        #cols = fw.dtype.names
        #for colname in cols:
        #    specd[colname][...] = specio[colname][...]
        #del specio
        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)

    else:
        tmask = np.ones(len(tiles)).astype('bool')

    for tile,tdate in zip(tiles[tmask]['TILEID'],tiles[tmask]['THRUDATE']):
        tdate = str(tdate)
        
        tspec = combspecdata_alt(tile,tdate,coaddir=coaddir)
        if tspec:
            tspec['TILEID'] = tile
            tspec = np.array(tspec)
            
            if s == 0:
                specd = tspec
                s = 1
            new = np.empty(len(tspec),dtype=specd.dtype)
            cols = specd.dtype.names
            for colname in cols:
                new[colname][...] = tspec[colname][...]
            tl.append(new)

            n += 1
            print(tile,n,len(tiles[tmask]))#,len(specd))
        else:
            print(str(tile)+' failed')
            nfail += 1
    specd = np.hstack(tl)
    kp = (specd['TARGETID'] > 0)
    specd = specd[kp]
    print('total number of failures was '+str(nfail))
    if n > 0:
        #specd.write(outf,format='fits', overwrite=True)
        fitsio.write(outf,specd,clobber=True)
        return True
    else:
        return False


def combspecdata(tile,zdate,tdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/archive/',md='' ):
    #put data from different spectrographs together, one table for fibermap, other for z
    zdate = str(zdate)
    specs = []
    #find out which spectrograph have data
    zfn = 'zbest'
    zhdu = 'ZBEST'
    shdu = 'SCORES'
    if int(tdate) >  20210730:
        zfn = 'redrock'
        zhdu = 'REDSHIFTS'
        #shdu = 'TSNR2'


    for si in range(0,10):
        ff = coaddir+str(tile)+'/'+zdate+'/'+zfn+'-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
        if os.path.isfile(ff):
            fq = coaddir+str(tile)+'/'+zdate+'/zmtl-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
            if os.path.isfile(fq):

                specs.append(si)
            else:
                print('did not find '+fq)
        elif zfn == 'zbest':
            zfnt = 'redrock'
            ff = coaddir+str(tile)+'/'+zdate+'/'+zfnt+'-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
            if os.path.isfile(ff):
                fq = coaddir+str(tile)+'/'+zdate+'/zmtl-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
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
    tsl = []
    tql = []
    tspecl = []
    tfl = []
    for i in range(0,len(specs)):
        tn = Table.read(coaddir+str(tile)+'/'+zdate+'/'+zfn+'-'+str(specs[i])+'-'+str(tile)+'-thru'+tdate+'.fits',hdu=zhdu)
        tspecl.append(tn)
        tnq = Table.read(coaddir+str(tile)+'/'+zdate+'/zmtl-'+str(specs[i])+'-'+str(tile)+'-thru'+tdate+'.fits')
        tql.append(tnq)
        tnf = Table.read(coaddir+str(tile)+'/'+zdate+'/'+zfn+'-'+str(specs[i])+'-'+str(tile)+'-thru'+tdate+'.fits',hdu='FIBERMAP')
        tfl.append(tnf)
        tns = Table.read(coaddir+str(tile)+'/'+zdate+'/coadd-'+str(specs[i])+'-'+str(tile)+'-thru'+tdate+'.fits',hdu=shdu)
        tsl.append(tns)
#         if i == 0:
#            tspec = tn
#            tq = tnq
#            tf = tnf
#            ts = tns
#         else:
#             ts = vstack([ts,tns],metadata_conflicts='silent')
#             tq = vstack([tq,tnq],metadata_conflicts='silent')
#             tspec = vstack([tspec,tn],metadata_conflicts='silent')
#             tf = vstack([tf,tnf],metadata_conflicts='silent')

    ts = vstack(tsl,metadata_conflicts='silent')
    tq = vstack(tql,metadata_conflicts='silent')
    tspec = vstack(tspecl,metadata_conflicts='silent')
    tf = vstack(tfl,metadata_conflicts='silent')


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

def combspecdata_alt(tile,tdate,coaddir='',md='' ):
    #put data from different spectrographs together, one table for fibermap, other for z
    tdate = str(tdate)
    specs = []
    #find out which spectrograph have data
    #zfn = 'zbest'
    #zhdu = 'ZBEST'
    shdu = 'SCORES'
    #if int(tdate) >  20210730:
    zfn = 'redrock'
    zhdu = 'REDSHIFTS'
    #shdu = 'TSNR2'


    for si in range(0,10):
        ff = coaddir+str(tile)+'/'+tdate+'/'+zfn+'-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
        if os.path.isfile(ff):
            specs.append(si)
        else:
            print('did not find '+ff)
    print('spectrographs with data on tile '+str(tile)+':')
    print(specs)
    if len(specs) == 0:
        return None
    #tsl = []
    #tql = []
    tspecl = []
    tfl = []
    for i in range(0,len(specs)):
        tn = Table.read(coaddir+str(tile)+'/'+tdate+'/'+zfn+'-'+str(specs[i])+'-'+str(tile)+'-thru'+tdate+'.fits',hdu=zhdu)
        tspecl.append(tn)
        tnf = Table.read(coaddir+str(tile)+'/'+tdate+'/'+zfn+'-'+str(specs[i])+'-'+str(tile)+'-thru'+tdate+'.fits',hdu='FIBERMAP')
        tfl.append(tnf)
#         if i == 0:
#            tspec = tn
#            tq = tnq
#            tf = tnf
#            ts = tns
#         else:
#             ts = vstack([ts,tns],metadata_conflicts='silent')
#             tq = vstack([tq,tnq],metadata_conflicts='silent')
#             tspec = vstack([tspec,tn],metadata_conflicts='silent')
#             tf = vstack([tf,tnf],metadata_conflicts='silent')

    #ts = vstack(tsl,metadata_conflicts='silent')
    #tq = vstack(tql,metadata_conflicts='silent')
    tspec = vstack(tspecl,metadata_conflicts='silent')
    tf = vstack(tfl,metadata_conflicts='silent')


    tf = unique(tf,keys=['TARGETID'])
    if md == '4combtar': #target files should contain the rest of the info
        tf.keep_columns(['FIBERASSIGN_X','FIBERASSIGN_Y','TARGETID','LOCATION','FIBERSTATUS','PRIORITY','DELTA_X','DELTA_Y','PSF_TO_FIBER_SPECFLUX','EXPTIME','OBJTYPE'])
    #tq.keep_columns(['TARGETID','Z_QN','Z_QN_CONF','IS_QSO_QN','ZWARN'])
    #tq['ZWARN'].name = 'ZWARN_MTL'
    tspec = join(tspec,tf,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')
    #tspec = join(tspec,ts,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')
    #tspec = join(tspec,tq,keys=['TARGETID'],join_type='left',metadata_conflicts='silent')

    print(len(tspec),len(tf))
    #tspec['LOCATION'] = tf['LOCATION']
    #tspec['FIBERSTATUS'] = tf['FIBERSTATUS']
    #tspec['PRIORITY'] = tf['PRIORITY']
    return tspec


def combtile_em(tiles,outf='',md='',prog='dark',redo='n'):
    s = 0
    n = 0
    nfail = 0
    guadtid = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/datcomb_'+prog+'_spec_zdone.fits',columns=['TILEID'])
    #print(len(guadtid))
    guadtid = np.unique(guadtid['TILEID'])

    if os.path.isfile(outf) and redo == 'n':
        #specd = Table.read(outf)
        specd = fitsio.read(outf)
        #dt = specd.dtype
        #specd = np.empty(len(specio),dtype=dt)
        #cols = fw.dtype.names
        #for colname in cols:
        #    specd[colname][...] = specio[colname][...]
        #del specio
        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)

    
        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)
    else:
        tmask = np.ones(len(tiles)).astype('bool')

    newl = []
    for tile,zdate,tdate in zip(tiles[tmask]['TILEID'],tiles[tmask]['ZDATE'],tiles[tmask]['THRUDATE']):
        tdate = str(tdate)
        tspec = None
        if np.isin(tile,guadtid):
            #if specver ==
            #tspec = combzmtl(tile,zdate,tdate)
            tspec = combEMdata_guad(tile,tdate)
        else:
            tnm = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/emtiles/emline-'+str(tile)+'.fits'
            if os.path.isfile(tnm):
                tspec = fitsio.read(tnm)
        if tspec is not None:
            tspec = np.array(tspec)

            if s == 0:
                specd = tspec
                s = 1
            else:
                #specd = vstack([specd,tspec],metadata_conflicts='silent')
                #column order got mixed up
                new = np.empty(len(tspec),dtype=specd.dtype)
                cols = specd.dtype.names
                for colname in cols:
                    new[colname][...] = tspec[colname][...]

                newl.append(new)
                #specd = np.hstack((specd,new))
            #specd.sort('TARGETID')
            #kp = (specd['TARGETID'] > 0)
            #specd = specd[kp]

            n += 1
            print(tile,n,len(tiles[tmask]),len(specd))
        else:
            print(str(tile)+' failed')
            nfail += 1
    print('total number of failures was '+str(nfail))
    newtot = np.hstack(newl)
    specd = np.hstack((specd,newtot))
    kp = (specd['TARGETID'] > 0)
    specd = specd[kp]
    if n > 0:
        #specd.write(outf,format='fits', overwrite=True)
        fitsio.write(outf,specd,clobber=True)
        return True
    else:
        return False

def combtile_em_alt(tiles,outf='',md='',prog='dark',coaddir=''):

    s = 0
    n = 0
    nfail = 0

    if os.path.isfile(outf) and redo == 'n':
        #specd = Table.read(outf)
        specd = fitsio.read(outf)
        #dt = specd.dtype
        #specd = np.empty(len(specio),dtype=dt)
        #cols = fw.dtype.names
        #for colname in cols:
        #    specd[colname][...] = specio[colname][...]
        #del specio
        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)

    
        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)
    else:
        tmask = np.ones(len(tiles)).astype('bool')

    for tile,tdate in zip(tiles[tmask]['TILEID'],tiles[tmask]['THRUDATE']):
        tdate = str(tdate)
        tspec = None
        tspec = combEMdata_rel(tile,tdate,coaddir=coaddir)
        if tspec is not None:
            tspec = np.array(tspec)

            if s == 0:
                specd = tspec
                s = 1
            else:
                #specd = vstack([specd,tspec],metadata_conflicts='silent')
                #column order got mixed up
                new = np.empty(len(tspec),dtype=specd.dtype)
                cols = specd.dtype.names
                for colname in cols:
                    new[colname][...] = tspec[colname][...]

                #specd = np.hstack((specd,tspec))
                specd = np.hstack((specd,new))
            #specd.sort('TARGETID')
            kp = (specd['TARGETID'] > 0)
            specd = specd[kp]

            n += 1
            print(tile,n,len(tiles[tmask]),len(specd))
        else:
            print(str(tile)+' '+coaddir+' failed')
            nfail += 1
    print('total number of failures was '+str(nfail))
    if n > 0:
        #specd.write(outf,format='fits', overwrite=True)
        fitsio.write(outf,specd,clobber=True)
        return True
    else:
        return False

def combtile_em_daily(tiles,outf='',md='',prog='dark',coaddir='',redo='n'):

    s = 0
    n = 0
    nfail = 0

    if os.path.isfile(outf) and redo == 'n':
        #specd = Table.read(outf)
        specd = fitsio.read(outf)
        #dt = specd.dtype
        #specd = np.empty(len(specio),dtype=dt)
        #cols = fw.dtype.names
        #for colname in cols:
        #    specd[colname][...] = specio[colname][...]
        #del specio
        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)

    
        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)
    else:
        tmask = np.ones(len(tiles)).astype('bool')

    for tile,zdate,tdate in zip(tiles[tmask]['TILEID'],tiles[tmask]['ZDATE'],tiles[tmask]['THRUDATE']):
        tdate = str(tdate)
        tspec = None
        tspec = combEMdata_daily(tile,str(zdate),tdate)
        if tspec is not None:
            tspec = np.array(tspec)

            if s == 0:
                specd = tspec
                s = 1
            else:
                #specd = vstack([specd,tspec],metadata_conflicts='silent')
                #column order got mixed up
                new = np.empty(len(tspec),dtype=specd.dtype)
                cols = specd.dtype.names
                for colname in cols:
                    new[colname][...] = tspec[colname][...]

                #specd = np.hstack((specd,tspec))
                specd = np.hstack((specd,new))
            #specd.sort('TARGETID')
            kp = (specd['TARGETID'] > 0)
            specd = specd[kp]

            n += 1
            print(tile,n,len(tiles[tmask]),len(specd))
        else:
            print(str(tile)+' failed')
            nfail += 1
    print('total number of failures was '+str(nfail))
    if n > 0:
        #specd.write(outf,format='fits', overwrite=True)
        fitsio.write(outf,specd,clobber=True)
        return True
    else:
        return False



def combEMdata_rel(tile,tdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/guadalupe/tiles/cumulative/'):
    remcol = ['Z', 'ZWARN', 'SPECTYPE', 'DELTACHI2', 'TARGET_RA', 'TARGET_DEC', 'OBJTYPE']
    zfn = 'emline'
    dl = []
    for si in range(0,10):
        ff = coaddir+str(tile)+'/'+tdate+'/'+zfn+'-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
        fz = coaddir+str(tile)+'/'+tdate+'/redrock'+'-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
        if os.path.isfile(ff):
            d = Table(fitsio.read(ff))
            fm = fitsio.read(fz,ext='FIBERMAP',columns=['LOCATION'])
            d['LOCATION'] = fm['LOCATION']
            dl.append(d)
    if len(dl) > 0:
        dt = vstack(dl,metadata_conflicts='silent')
        dt['TILEID'] = int(tile)
        dt.remove_columns(remcol)
        return dt
    else:
        return None

def combEMdata_daily(tile,zdate,tdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/archive/'):
    remcol = ['Z', 'ZWARN', 'SPECTYPE', 'DELTACHI2', 'TARGET_RA', 'TARGET_DEC', 'OBJTYPE']
    zfn = 'emline'
    dl = []
    for si in range(0,10):
        ff = coaddir+str(tile)+'/'+zdate+'/'+zfn+'-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
        fz = coaddir+str(tile)+'/'+zdate+'/redrock'+'-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
        if os.path.isfile(ff):
            d = Table(fitsio.read(ff))
            fm = fitsio.read(fz,ext='FIBERMAP',columns=['LOCATION'])
            d['LOCATION'] = fm['LOCATION']
            dl.append(d)
    if len(dl) > 0:
        dt = vstack(dl,metadata_conflicts='silent')
        dt['TILEID'] = tile
        dt.remove_columns(remcol)
        return dt
    else:
        return None


def combEMdata_daily_old(tile,zdate,tdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/archive/',outf='temp.fits'):
    from desispec.io.emlinefit import read_emlines_inputs
    from desispec.emlinefit import get_emlines

    allems = ['OII','HDELTA','HGAMMA','HBETA','OIII','HALPHA']
    props = ['FLUX','FLUX_IVAR','SIGMA','SIGMA_IVAR','CONT','CONT_IVAR','SHARE','SHARE_IVAR','EW','EW_IVAR','CHI2','NDOF']
    zdate = str(zdate)
    specs = []
    #find out which spectrograph have data
    zfn = 'zbest'
    zhdu = 'ZBEST'
    shdu = 'SCORES'
    if int(tdate) >  20210730:
        zfn = 'redrock'
    dl = []
    for si in range(0,10):
        ff = coaddir+str(tile)+'/'+zdate+'/'+zfn+'-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
        cf = coaddir+str(tile)+'/'+zdate+'/coadd-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
        if os.path.isfile(ff) and os.path.isfile(cf):
            d = Table.read(ff, "REDSHIFTS")
            df = Table.read(ff, "FIBERMAP")
            tids, zs = d["TARGETID"], d["Z"]            
            _, _, waves, fluxes, ivars = read_emlines_inputs(ff, cf, targetids=tids)
            emdict = get_emlines(zs, waves, fluxes, ivars, emnames=allems)
            t = Table()
            for em in allems:
                for prop in props:
                    t[em+'_'+prop] = emdict[em][prop]
            t['TARGETID'] = d['TARGETID']
            t['LOCATION'] = df['LOCATION']
            dl.append(t)
    if len(dl) > 0:
        dt = vstack(dl,metadata_conflicts='silent')
        dt['TILEID'] = tile
        dt.write(outf,format='fits',overwrite=True)
    else:
        print('no data to combine for tile '+str(tile))    
    #outf = outdir+'emline-'+str(tile)+'.fits'
    #'/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/emtiles/'
    
    #return dt

def combQSOdata_alt(tile,zdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/archive/',cols=None ):
    #combine for some rerun, won't be from archive so dates are the same
    from LSS.qso_cat_utils import qso_catalog_maker
    #put data from different spectrographs together, one table for fibermap, other for z
    zdate = str(zdate)
    tdate = zdate
    specs = []
    #find out which spectrograph have data
    #zfn = 'zbest'
    #zhdu = 'ZBEST'
    #shdu = 'SCORES'
    #if int(tdate) >  20210730:
    zfn = 'redrock'
    zhdu = 'REDSHIFTS'
        #shdu = 'TSNR2'


    for si in range(0,10):
        ff = coaddir+str(tile)+'/'+zdate+'/'+zfn+'-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
        if os.path.isfile(ff):
            fq = coaddir+str(tile)+'/'+zdate+'/qso_mgii-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
            fqn = coaddir+str(tile)+'/'+zdate+'/qso_qn-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
            if os.path.isfile(fq) and os.path.isfile(fqn):

                specs.append(si)
            else:
                print('did not find '+fq+' or '+fqn)
        elif zfn == 'zbest':
            zfnt = 'redrock'
            ff = coaddir+str(tile)+'/'+zdate+'/'+zfnt+'-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
            if os.path.isfile(ff):
                fq = coaddir+str(tile)+'/'+zdate+'/qso_mgii-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
                fqn = coaddir+str(tile)+'/'+zdate+'/qso_qn-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
                zfn = zfnt
                zhdu = 'REDSHIFTS'
                if os.path.isfile(fq):

                    specs.append(si)
                else:
                    print('did not find '+fq+' '+fqn)
            else:
                print('did not find '+ff)
        else:
            print('did not find '+ff)
    print('spectrographs with data on tile '+str(tile)+':')
    print(specs)
    if len(specs) == 0:
        return None
    qsocats = []    
    for i in range(0,len(specs)):
        rr = coaddir+str(tile)+'/'+zdate+'/'+zfn+'-'+str(specs[i])+'-'+str(tile)+'-thru'+tdate+'.fits'
        mgii = coaddir+str(tile)+'/'+zdate+'/'+'qso_mgii'+'-'+str(specs[i])+'-'+str(tile)+'-thru'+tdate+'.fits'
        qn = coaddir+str(tile)+'/'+zdate+'/'+'qso_qn'+'-'+str(specs[i])+'-'+str(tile)+'-thru'+tdate+'.fits'
        old_extname_redrock = True if zhdu == 'ZBEST' else False
        old_extname_for_qn = False #if int(tdate) >= 20220118 else True
        qso_cati = Table.from_pandas(qso_catalog_maker(rr, mgii, qn, old_extname_redrock, old_extname_for_qn))
        #qso_cati = Table(qso_catalog_maker(rr, mgii, qn, old_extname_redrock, old_extname_for_qn))
        qsocats.append(qso_cati)
        #if i == 0:
        #    qso_cat = qso_cati
        #else:
        #    qso_cat = vstack([qso_cat,qso_cati],metadata_conflicts='silent')

    qso_cat = vstack(qsocats,metadata_conflicts='silent')
    qso_cat['TILEID'] = tile
    #print(qso_cat.dtype.names)
    if cols is not None:
        qso_cat = Table(qso_cat)
        qso_cat.keep_columns(cols)

    return qso_cat


def combQSOdata(tile,zdate,tdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/archive/',cols=None ):
    from LSS.qso_cat_utils import qso_catalog_maker
    #put data from different spectrographs together, one table for fibermap, other for z
    zdate = str(zdate)
    specs = []
    #find out which spectrograph have data
    zfn = 'zbest'
    zhdu = 'ZBEST'
    shdu = 'SCORES'
    if int(tdate) >  20210730:
        zfn = 'redrock'
        zhdu = 'REDSHIFTS'
        #shdu = 'TSNR2'


    for si in range(0,10):
        ff = coaddir+str(tile)+'/'+zdate+'/'+zfn+'-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
        if os.path.isfile(ff):
            fq = coaddir+str(tile)+'/'+zdate+'/qso_mgii-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
            fqn = coaddir+str(tile)+'/'+zdate+'/qso_qn-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
            if os.path.isfile(fq) and os.path.isfile(fqn):

                specs.append(si)
            else:
                print('did not find '+fq+' or '+fqn)
        elif zfn == 'zbest':
            zfnt = 'redrock'
            ff = coaddir+str(tile)+'/'+zdate+'/'+zfnt+'-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
            if os.path.isfile(ff):
                fq = coaddir+str(tile)+'/'+zdate+'/qso_mgii-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
                fqn = coaddir+str(tile)+'/'+zdate+'/qso_qn-'+str(si)+'-'+str(tile)+'-thru'+tdate+'.fits'
                zfn = zfnt
                zhdu = 'REDSHIFTS'
                if os.path.isfile(fq):

                    specs.append(si)
                else:
                    print('did not find '+fq+' '+fqn)
            else:
                print('did not find '+ff)
        else:
            print('did not find '+ff)
    print('spectrographs with data on tile '+str(tile)+':')
    print(specs)
    if len(specs) == 0:
        return None
    qsocats = []    
    for i in range(0,len(specs)):
        rr = coaddir+str(tile)+'/'+zdate+'/'+zfn+'-'+str(specs[i])+'-'+str(tile)+'-thru'+tdate+'.fits'
        mgii = coaddir+str(tile)+'/'+zdate+'/'+'qso_mgii'+'-'+str(specs[i])+'-'+str(tile)+'-thru'+tdate+'.fits'
        qn = coaddir+str(tile)+'/'+zdate+'/'+'qso_qn'+'-'+str(specs[i])+'-'+str(tile)+'-thru'+tdate+'.fits'
        old_extname_redrock = True if zhdu == 'ZBEST' else False
        old_extname_for_qn = False #if int(tdate) >= 20220118 else True
        qso_cati = Table.from_pandas(qso_catalog_maker(rr, mgii, qn, old_extname_redrock, old_extname_for_qn, update_qn_zwarn = False))
        #qso_cati = Table(qso_catalog_maker(rr, mgii, qn, old_extname_redrock, old_extname_for_qn))
        qsocats.append(qso_cati)
        #if i == 0:
        #    qso_cat = qso_cati
        #else:
        #    qso_cat = vstack([qso_cat,qso_cati],metadata_conflicts='silent')

    qso_cat = vstack(qsocats,metadata_conflicts='silent')
    qso_cat['TILEID'] = tile
    #print(qso_cat.dtype.names)
    if cols is not None:
        qso_cat = Table(qso_cat)
        qso_cat.keep_columns(cols)

    return qso_cat


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

def get_skystd(tile,zdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/',clip=3. ):
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
    #print('spectrographs with data:')
    #print(specs)
    if len(specs) == 0:
        return None
    
    cam_stats = {}
    for camera in ["B", "R", "Z"]:
        cam_stats[camera] = []
    
    for i in range(0,len(specs)):
        fn = coaddir+str(tile)+'/'+zdate+'/coadd-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits'
        h = fits.open(fn)
        fm = Table(h["FIBERMAP"].data)
        sel = fm["OBJTYPE"] == "SKY"
        #nsky = sel.sum()
        #snrs = np.zeros((sel.sum(), 0))
        #tab_petal = Table()
        #tab_petal['TILEID'] = int(tile)
        
        for camera in ["B", "R", "Z"]:
            fluxes = h["{}_FLUX".format(camera)].data[sel, :]
            ivars = h["{}_IVAR".format(camera)].data[sel, :]
            snrs = fluxes * np.sqrt(ivars) #array of snr per pixel, one row per spectrum
            stds = np.std(snrs,axis=1) #array of the standard deviations of snr per spectrum
            disp_std = (stds-np.std(stds))/np.mean(stds) #how many standard deviations is each spectrum away for petal/camera
            snrs = snrs[disp_std < clip] #remove any spectra that are more than 3 sigma
            snrs = snrs.flatten()
            snrs = snrs[snrs != 0] #should remove any masked pixels?
            disp = snrs.std()
            if disp > 2 or disp == 0:
                print('tile '+str(tile)+' petal '+str(i)+' dispersion is '+str(disp) )
            cam_stats[camera].append(disp)
    tab = Table()
    tab['TILEID'] = np.ones(len(specs))*int(tile)
    tab['PETAL_LOC'] = np.array(specs,dtype=int)
    for camera in ["B", "R", "Z"]:
        tab['SKY_DISP_'+camera] = np.array(cam_stats[camera])
    return tab

def combtile_skystd(tiles,outf='',md='',specver='daily',redo='n',specrel='guadalupe',clip=3):
    s = 0
    n = 0
    nfail = 0
    tl = []
    if os.path.isfile(outf) and redo == 'n':
        specd = Table.read(outf)
        tl.append(specd)
        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)

    else:
        tmask = np.ones(len(tiles)).astype('bool')

    newtabs = []
    for tile,zdate,tdate in zip(tiles[tmask]['TILEID'],tiles[tmask]['ZDATE'],tiles[tmask]['THRUDATE']):
        tdate = str(tdate)
        coaddir='/global/cfs/cdirs/desi/spectro/redux/'+specver+'/tiles/cumulative/'
        tile_data = get_skystd(tile,tdate,coaddir,clip)
        newtabs.append(tile_data)
    newtabs = vstack(newtabs)
    if s == 1:
        specd.write(outf,format='fits',overwrite=True)
    else:
        newtabs.write(outf,format='fits',overwrite=True)

    
    return True

def combtile_petalqa(tiles,outf='',md='',specver='daily',redo='y'):
    s = 0
    n = 0
    nfail = 0
    tl = []
    if os.path.isfile(outf) and redo == 'n':
        specd = Table.read(outf)
        tl.append(specd)
        s = 1
        tdone = np.unique(specd['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)

    else:
        tmask = np.ones(len(tiles)).astype('bool')

    newtabs = []
    for tile,zdate,tdate in zip(tiles[tmask]['TILEID'],tiles[tmask]['ZDATE'],tiles[tmask]['THRUDATE']):
        tdate = str(tdate)
        coaddir='/global/cfs/cdirs/desi/spectro/redux/'+specver+'/tiles/cumulative/'
        fn = coaddir  +str(tile)+'/'+tdate+'/tile-qa-'+str(tile)+'-thru'+tdate+'.fits'
        tile_data = Table(fitsio.read(fn,ext= "PETALQA"))
        tile_data['TILEID'] = np.ones(len(tile_data),dtype=int)*tile
        newtabs.append(tile_data)
    newtabs = vstack(newtabs)
    if s == 1:
        specd = vstack([specd,newtabs]) 
        specd.write(outf,format='fits',overwrite=True)
    else:
        newtabs.write(outf,format='fits',overwrite=True)
    return True



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

def get_tiletab(tile_row,tarcol=['RA','DEC','TARGETID','DESI_TARGET','BGS_TARGET','MWS_TARGET','SUBPRIORITY','PRIORITY_INIT','TARGET_STATE','TIMESTAMP','ZWARN','PRIORITY']):
    from desitarget.io import read_targets_in_tiles
    
    tile = tile_row['TILEID'][0]
    ts = str(tile).zfill(6)
    faf = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
    fht = fitsio.read_header(faf)
    mdir = '/global/cfs/cdirs/desi'+fht['MTL'][8:]+'/'
    if mdir == '/global/cfs/cdirs/desi/survey/ops/staging/mtl/main/dark/':
        mdir = '/global/cfs/cdirs/desi/target/catalogs/mtl/1.0.0/mtl/main/dark/'
    if mdir == '/global/cfs/cdirs/desi/survey/ops/staging/mtl/main/bright/':
        mdir = '/global/cfs/cdirs/desi/target/catalogs/mtl/1.0.0/mtl/main/bright/'
    #wt = tiles['TILEID'] == tile
    tars = read_targets_in_tiles(mdir,tile_row,mtl=True,isodate=fht['MTLTIME'])
    #tars.keep_columns(tarcols)
    tars = tars[[b for b in tarcol]]

    tt = Table.read(faf,hdu='POTENTIAL_ASSIGNMENTS')
    tars = join(tars,tt,keys=['TARGETID'])
    tars['TILEID'] = tile
    tars.remove_columns(['ZWARN'])
    return tars


def combtiles_wdup(tiles,fout='',tarcol=['RA','DEC','TARGETID','DESI_TARGET','BGS_TARGET','MWS_TARGET','SUBPRIORITY','PRIORITY_INIT','TARGET_STATE','TIMESTAMP','ZWARN','PRIORITY']):
    from desitarget.io import read_targets_in_tiles
    s = 0
    n = 0
    tl = []
    if os.path.isfile(fout):
        tars = Table.read(fout)
        tl.append(tars)
        #s = 1
        tdone = np.unique(tars['TILEID'])
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
        tl.append(tars)
        #if s == 0:
        #    tarsn = tars
        #    s = 1
        #else:
        #    tarsn = vstack([tarsn,tars],metadata_conflicts='silent')
        #tarsn.sort('TARGETID')
        n += 1
        print(tile,n,len(tiles[tmask]))#,len(tarsn))
    if np.sum(tmask) > 0:
        print('about to stack')
        tarsn = vstack(tl)
        tarsn.sort('TARGETID')
        common.write_LSS(tarsn,fout)
        #tarsn.write(fout,format='fits', overwrite=True)
    else:
        print('nothing to update, done')

def combtiles_wdup_hp(hpx,tiles,fout='',tarcol=['RA','DEC','TARGETID','DESI_TARGET','BGS_TARGET','MWS_TARGET','SUBPRIORITY','PRIORITY_INIT','TARGET_STATE','TIMESTAMP','ZWARN','PRIORITY'],logger=None):
    import desimodel.footprint as foot
    from desitarget.io import read_targets_in_tiles
    import LSS.common_tools as common
    s = 0
    n = 0

    tarsn = None
    tls = foot.pix2tiles(8,[hpx],tiles)
    if os.path.isfile(fout):
        tarsn = Table.read(fout)
        s = 1
        tdone = np.unique(tarsn['TILEID'])
        tmask = ~np.isin(tls['TILEID'],tdone)
    else:
        tmask = np.ones(len(tls)).astype('bool')
    common.printlog('there are potentially '+str(len(tls[tmask]))+' to get updates from, out of a possible '+str(len(tls))+' overlapping this pixel',logger)
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
        tars = read_targets_in_tiles(mdir,tls[wt],mtl=True,isodate=fht['MTLTIME'])
        #tars.keep_columns(tarcols)
        tars = tars[[b for b in tarcol]]
        if 'MTL2' in fht.keys():
            mdir2 = '/global/cfs/cdirs/desi'+fht['MTL2'][8:]+'/'
            tars2 = read_targets_in_tiles(mdir2,tls[wt],mtl=True,isodate=fht['MTLTIME'])
            tars2 = tars2[[b for b in tarcol]]
            #common.printlog(str(tars.dtype.names),logger)
            #common.printlog(str(tars2.dtype.names),logger)
            common.printlog('adding '+str(len(tars2))+ ' entries from MTL2 to tile '+str(tile),logger)
            #tars = vstack([tars,tars2])
            tars = np.concatenate([tars,tars2])
        theta, phi = np.radians(90-tars['DEC']), np.radians(tars['RA'])
        tpix = hp.ang2pix(8,theta,phi,nest=True)
        sel = tpix == hpx
        tars = tars[sel]
        tt = Table.read(faf,hdu='POTENTIAL_ASSIGNMENTS')
        if np.sum(np.isin(tt['TARGETID'],tars['TARGETID'])) > 0:
            tars = join(tars,tt,keys=['TARGETID'])
            tars['TILEID'] = tile
            tars.remove_columns(['ZWARN'])
            if s == 0:
                tarsn = tars
                s = 1
            else:
                tarsn = vstack([tarsn,tars],metadata_conflicts='silent')
            tarsn.sort('TARGETID')

            #print(tile,n,len(tls[tmask]),len(tarsn))

        else:
            common.printlog('no overlapping targetid for tile '+str(tile),logger)
        n += 1
    if tarsn is not None and n > 0:
        tarsn.write(fout,format='fits', overwrite=True)
    else:
        if tarsn == None:
            print('did not find any targets actually in this pixel '+str(hpx))
        else:
            print('no tiles to update for this pixel '+str(hpx))

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


def get_specdat(indir,pd,ver='daily',badfib=None):
    from desitarget.targetmask import zwarn_mask
    #indir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+specrel
    #zf = indir+'/datcomb_'+pd+'_tarspecwdup.fits'
    zf = indir+'/datcomb_'+pd+'_spec_zdone.fits'
    if ver == 'everest' or ver == 'guadalupe':
        zf = indir+'/datcomb_'+pd+'_tarspecwdup_zdone.fits'
    #if ver == 'daily' or ver == 'himalayas':
        
    print('getting spec data from '+zf)
    dz = Table.read(zf)
    #dz = fitsio.read(zf)
    selz = dz['ZWARN'] != 999999
    selz &= dz['ZWARN']*0 == 0 #just in case of nans
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
    if badfib is not None:
        bad = np.isin(fs['FIBER'],badfib)
        print('number at bad fibers '+str(sum(bad)))
        wfqa &= ~bad

    return fs[wfqa]

def cut_specdat(dz):
    import LSS.common_tools as common
    print('moved to common_tools')
    return common.cut_specdat(dz)
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
#     nomtl = nodata | badqa
#     wfqa = ~nomtl
#     return fs[wfqa]

def count_tiles_input_alt(fjg,logger=None):
    import LSS.common_tools as common
    '''
    take input array with require columns TARGETID TILEID
    return table with unique TARGETID and the number of tiles it showed up on (NTILE), the TILES 
    '''
    #fjg = Table(fjg)
    #fjg.keep_columns(['TARGETID','TILEID'])
    #fjg = np.array(fjg)
    #print(fjg.dtype.names)
    tidstot = fjg['TARGETID']
    tileidstot = fjg['TILEID']
    del fjg
    indsort = np.argsort(tidstot)
    tidstot = tidstot[indsort]
    tileidstot = tileidstot[indsort]

    tids,cnts = np.unique(tidstot,return_counts=True)
    common.printlog('counting tiles, going through '+str(len(tidstot))+' rows with '+str(len(tids))+' unique targetid',logger)
    cntold = 0
    tl = []
    nt = []
    for i in range(0,len(tids)):
        #get indices, given array was sorted
        cntstot = cntold + cnts[i]
        inds = cntold,cntstot 
        cntold = cntstot
        #get tileids for given targetid
        tls = tileidstot[inds[0]:inds[1]]
        tlsu = np.unique(tls)    
        nt.append(len(tlsu))
        tl.append("-".join(tlsu.astype(str)))
       

        if i%1000000 == 0:
            common.printlog(str(i),logger)
        
    
    tc = Table()
    tc['TARGETID'] = tids
    tc['NTILE'] = nt
    tc['TILES'] = tl
    #tc['TILELOCIDS'] = tli
    
    return tc


def count_tiles_input(fjg,logger=None):
    import LSS.common_tools as common
    '''
    take input array with require columns TARGETID TILEID TILELOCID
    return table with unique TARGETID and the number of tiles it showed up on (NTILE), the TILES and the TILELOCIDS
    '''
    fjg = Table(fjg)
    fjg.keep_columns(['TARGETID','TILEID','TILELOCID'])
    fjg = np.array(fjg)
    #print(fjg.dtype.names)
    fjg = fjg[np.argsort(fjg['TARGETID'])]

    tids = np.unique(fjg['TARGETID'])
    common.printlog('counting tiles, going through '+str(len(fjg))+' rows with '+str(len(tids))+' unique targetid',logger)
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
            #tlis.append(fjg[i]['TILELOCID'])
            i += 1
            if i == len(fjg):
                break
        nloc.append(nli)
        tlsu = np.unique(tls)
        #tlisu = np.unique(tlis)
        nt.append(len(tlsu))
        tl.append("-".join(tlsu.astype(str)))
        #tli.append("-".join(tlisu.astype(str)))

        if ti%1000000 == 0:
            common.printlog(str(ti),logger)
        ti += 1
    del fjg
    tc = Table()
    tc['TARGETID'] = tids
    tc['NTILE'] = nt
    tc['TILES'] = tl
    #tc['TILELOCIDS'] = tli
    
    return tc

def count_tiles_better(dr,pd,rann=0,specrel='daily',fibcol='COADD_FIBERSTATUS',px=False,survey='main',indir=None,gtl=None,badfib=None, prog_ = 'dark'):
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

    if indir is None:
        indir = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/'+specrel
    ps = pd
    if pd[:3] == 'LRG' or pd[:3] == 'ELG' or pd[:3] =='QSO':
        ps = 'dark'
    if pd[:3] == 'BGS' or pd[:3] == 'MWS_ANY':
        ps = 'bright'
    if gtl is None:
        print('getting good tileloc')
        fs = get_specdat(indir,ps,specrel,badfib=badfib)

        stlid = 10000*fs['TILEID'] +fs['LOCATION']
        gtl = np.unique(stlid)

    if dr == 'dat':
        fj = fitsio.read(indir+'/datcomb_'+pd+'_tarspecwdup_zdone.fits',columns=['TARGETID','TILEID','TILELOCID'])
        #outf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/datcomb_'+pd+'ntileinfo.fits'
    if dr == 'ran':
        if px:
            fj = fitsio.read(indir+'/healpix/rancomb_'+str(rann)+pd+'_'+str(px)+'_wdupspec_zdone.fits',columns=['TARGETID','TILEID','TILELOCID'])
        else:
            fj = fitsio.read(indir+'/rancomb_'+str(rann)+pd+'wdupspec_zdone.fits',columns=['TARGETID','TILEID','TILELOCID'])
    
    if (dr == 'mock'):
        fj = fitsio.read(indir+ '/comb' + prog_ +'_wdupspec_zdone_' + pd + '.fits' , columns = ['TARGETID', 'TILEID', 'TILELOCID'])
        
        #outf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random'+str(rann)+'/rancomb_'+pd+'ntileinfo.fits'
    wg = np.isin(fj['TILELOCID'],gtl)
    fjg = fj[wg]
    del fj
    fjg = fjg[np.argsort(fjg['TARGETID'])]

    tids = np.unique(fjg['TARGETID'])
    print('going through '+str(len(fjg))+' rows with '+str(len(tids))+' unique targetid')
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



def count_tiles_better_px(dr,pd,gtl,rann=0,specrel='daily',fibcol='COADD_FIBERSTATUS',px=None,survey='main'):
    '''
    from files with duplicates that have already been sorted by targetid, quickly go
    through and get the multi-tile information
    dr is either 'dat' or 'ran'
    returns file with TARGETID,NTILE,TILES,TILELOCIDS
    '''

    if dr == 'dat':
        fj = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/'+specrel+'/datcomb_'+pd+'_tarspecwdup_zdone.fits')
        #outf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/datcomb_'+pd+'ntileinfo.fits'
    if dr == 'ran':
        if px is not None:
            fj = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/'+specrel+'/healpix/rancomb_'+str(rann)+pd+'_'+str(px)+'_wdupspec_zdone.fits')
        else:
            fj = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/'+specrel+'/rancomb_'+str(rann)+pd+'wdupspec_zdone.fits')

        #outf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random'+str(rann)+'/rancomb_'+pd+'ntileinfo.fits'
    wg = np.isin(fj['TILELOCID'],gtl)
    fjg = fj[wg]
    fjg = fjg[np.argsort(fjg['TARGETID'])]

    tids = np.unique(fjg['TARGETID'])
    print('going through '+str(len(fjg))+' rows with '+str(len(tids))+' unique targetid')
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


def combran_wdup(tiles,rann,randir,outf,keepcols=[],redo=True):

    s = 0
    td = 0
    #tiles.sort('ZDATE')
    print(len(tiles))
    #delcols = ['DESI_TARGET','BGS_TARGET','MWS_TARGET','SUBPRIORITY','OBSCONDITIONS','PRIORITY_INIT',\
    #'NUMOBS_INIT','SCND_TARGET','NUMOBS_MORE','NUMOBS','Z','ZWARN','TARGET_STATE','TIMESTAMP','VERSION','PRIORITY']
    #outf = randir+str(rann)+'/rancomb_'+tp+'wdup_Alltiles.fits'
    tldata = []
    if os.path.isfile(outf) and redo == False:
        fgu = Table(fitsio.read(outf))
        #tarsn.keep_columns(['RA','DEC','TARGETID''LOCATION','FIBER','TILEID'])
        s = 1
        tdone = np.unique(fgu['TILEID'])
        tmask = ~np.isin(tiles['TILEID'],tdone)
        tldata.append(fgu)
    else:
        tmask = np.ones(len(tiles)).astype('bool')
    
    tot = 0
    for tile in tiles[tmask]['TILEID']:
        ffa = randir+str(rann)+'/fba-'+str(tile).zfill(6)+'.fits'
        ffna = randir+str(rann)+'/tilenofa-'+str(tile)+'.fits'
        if os.path.isfile(ffa):
            fa = Table(fitsio.read(ffa,ext='FAVAIL'))

            ffna = Table(fitsio.read(ffna))
            fgun = join(fa,ffna,keys=['TARGETID'],join_type='left')
            #fgun.remove_columns(delcols)

            td += 1
            fgun['TILEID'] = int(tile)
            fgun.keep_columns(['RA','DEC','TARGETID','LOCATION','FIBER','TILEID'])
            tldata.append(fgun)
            #if s == 0:
            #    fgu = fgun
            #    s = 1
            #else:
            #    fgu = vstack([fgu,fgun],metadata_conflicts='silent')
            #fgu.sort('TARGETID')
            tot += len(fgun)
            print(tile,td, len(tiles),len(fa),len(fgun),tot)#, len(fgun),len(fgu))
            
        else:
            print('did not find '+ffa)

    fgu = vstack(tldata)
    print(len(fgu),tot)
    if len(tiles[tmask]['TILEID']) > 0:
        fgu.write(outf,format='fits', overwrite=True)
        rv = True
    else:
        rv = False
    return rv

def combran_wdupspec(rann,tp,lspecdir,specf,infile,keepcols=[],mask_coll=True,collf='', alt_out = None, mock_priority_mask = 'n', mock_tr = 'LRG',logger=None):
    from LSS.common_tools import write_LSS,write_LSS_scratchcp,printlog
    fgu = Table(fitsio.read(infile.replace('global','dvs_ro')))
    if mask_coll:
        printlog('length before masking collisions '+str(len(fgu)),logger)
        if 'COLLISION' in list(fgu.dtype.names):
            sel = fgu['COLLISION'] == 0
            fgu = fgu[sel]
        else:
            coll = Table(fitsio.read(collf.replace('global','dvs_ro')))
            fgu = setdiff(fgu,coll,keys=['TARGETID','LOCATION','TILEID'])
        printlog('length after masking collisions '+str(len(fgu)),logger)
    specf.keep_columns(keepcols)
    #specf.keep_columns(['ZWARN','LOCATION','TILEID','TILELOCID','FIBERSTATUS','FIBERASSIGN_X','FIBERASSIGN_Y','PRIORITY','DELTA_X','DELTA_Y','EXPTIME','PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B','TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z','TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
    printlog('joining to spec data',logger)
    fgu = join(fgu,specf,keys=['LOCATION','TILEID','FIBER'],join_type='left')
    #fgu.sort('TARGETID')
    if alt_out != None:
        #outf = alt_out + '/comb' + tp + '_' + mock_tr + '_'+ 'wdupspec_zdone.fits'
        outf = alt_out + '/comb' + tp + '_' + 'wdupspec_zdone.fits'
        if mock_priority_mask == 'y':
            if mock_tr == 'QSO':
                maxp = 3400
            else:
                maxp = 3200
            pr_mask = fgu["PRIORITY"] <= maxp
            fgu = fgu[pr_mask]
    else:
        outf = lspecdir+'/rancomb_'+str(rann)+tp+'wdupspec_zdone.fits'
    printlog('writing to '+outf,logger)
    write_LSS_scratchcp(fgu,outf,logger=logger)
    #fgu.write(outf,format='fits', overwrite=True)
    


def combran_wdup_hp(hpx,tiles,rann,randir,tp,lspecdir,specf,keepcols=[],outf='',redos=False):
    import desimodel.footprint as foot
    s = 0

    #tiles.sort('ZDATE')
    print(len(tiles))
    #delcols = ['DESI_TARGET','BGS_TARGET','MWS_TARGET','SUBPRIORITY','OBSCONDITIONS','PRIORITY_INIT',\
    #'NUMOBS_INIT','SCND_TARGET','NUMOBS_MORE','NUMOBS','Z','ZWARN','TARGET_STATE','TIMESTAMP','VERSION','PRIORITY']
    outf = randir+str(rann)+'/healpix/rancomb_'+tp+'_'+str(hpx)+'_wdup_Alltiles.fits'
    outfs = lspecdir+'healpix/rancomb_'+str(rann)+tp+'_'+str(hpx)+'_wdupspec_zdone.fits'
    tarsn = None
    tls = foot.pix2tiles(8,[hpx],tiles)
    if os.path.isfile(outf):
    #try:
        fgu = Table.read(outf)
        s = 1
        tdone = np.unique(fgu['TILEID'])
        tmask = ~np.isin(tls['TILEID'],tdone)
    else:
    #except:
        tmask = np.ones(len(tls)).astype('bool')

    td = len(tls[~tmask])
    #if td > 0:
    #if hpx == 79:
    #    print(outf)
    #    print(tls)
    #    print(tls[tmask])
    if len(tls[tmask]) > 0:

        for tile in tls[tmask]['TILEID']:
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
                print(tile,td, len(tls), len(fgun),len(fgu))
            else:
                print('did not find '+ffa)

        if len(tls[tmask]['TILEID']) > 0:
            fgu.write(outf,format='fits', overwrite=True)
        #specf = Table.read(lspecdir+'datcomb_'+tp+'_spec_zdone.fits')
        specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
        specf.keep_columns(keepcols)
        #specf.keep_columns(['ZWARN','LOCATION','TILEID','TILELOCID','FIBERSTATUS','FIBERASSIGN_X','FIBERASSIGN_Y','PRIORITY','DELTA_X','DELTA_Y','EXPTIME','PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B','TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z','TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
        fgu = join(fgu,specf,keys=['LOCATION','TILEID','FIBER'],join_type='left')
        fgu.sort('TARGETID')

        print(outfs)
        fgu.write(outfs,format='fits', overwrite=True)
        return True
    else:
        if os.path.isfile(outfs) == False or redos:
            fgu = fitsio.read(outf)
            sel = np.isin(fgu['TILEID'],tls['TILEID'])
            fgu = fgu[sel]
            specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
            specf.keep_columns(keepcols)
            #specf.keep_columns(['ZWARN','LOCATION','TILEID','TILELOCID','FIBERSTATUS','FIBERASSIGN_X','FIBERASSIGN_Y','PRIORITY','DELTA_X','DELTA_Y','EXPTIME','PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B','TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z','TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
            fgu = join(fgu,specf,keys=['LOCATION','TILEID','FIBER'],join_type='left')
            fgu.sort('TARGETID')

            print(outfs)
            fgu.write(outfs,format='fits', overwrite=True)
            if os.path.isfile(outfs) == False:
                return True
            else:
                return False
        else:
            print('no new data to add')
            return False



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

def mkfullran_prog(gtl,indir,rann,imbits,outf,pd,tlid_full=None,badfib=None,ftiles=None):
    import LSS.common_tools as common
    #import logging
    logger = logging.getLogger('LSSran')



    zf = indir.replace('global','dvs_ro')+'/rancomb_'+str(rann)+pd+'wdupspec_zdone.fits'
    logger.info('about to load '+zf)
    in_cols = ['LOCATION', 'FIBER', 'TARGETID', 'RA', 'DEC', 'TILEID', 'PRIORITY']#, 'TILELOCID']
    dz = Table(fitsio.read(zf,columns=in_cols))
    logger.info(dz.dtype.names)

    cols = list(dz.dtype.names)
 

    dz['TILELOCID'] = 10000*dz['TILEID'] +dz['LOCATION'] #reset it here in case was set by specdat and some matches were missing

    wg = np.isin(dz['TILELOCID'],gtl)
    if badfib is not None:
        bad = np.isin(dz['FIBER'],badfib)
        logger.info('number at bad fibers '+str(sum(bad)))
        wg &= ~bad

    dz['GOODHARDLOC'] = np.zeros(len(dz)).astype('bool')
    dz['GOODHARDLOC'][wg] = 1
    if ftiles is None:
        logger.info('counting tiles from dz with columns '+str(dz.dtype.names))
        dzpd = count_tiles_input(dz[wg],logger=logger)#.keep_columns(['TARGETID','TILEID','TILELOCID']))
    else:
        dzpd = Table.read(ftiles)


    logger.info(str(dz.dtype.names))
    p4sort = np.copy(dz['PRIORITY'] )
    sel =  p4sort*0 != 0
    p4sort[sel] = 999999
    logger.info(str(np.sum(sel))+' had nan priority')
    sel = p4sort <= 0
    p4sort[sel] = 1
    logger.info(str(np.sum(sel))+' had priority <= 0 set to 1 for sort')
    dz['sort'] =  dz['GOODHARDLOC'] + 1/p4sort
    #dz['sort'] =  dz['GOODPRI']*dz['GOODHARDLOC']*dz['ZPOSSLOC']#*(1+dz[tsnr])
    logger.info(dz.dtype.names)
    logger.info(str(rann)+' about to do sort')

    dz.sort('sort') #should allow to later cut on tsnr for match to data
    dz = unique(dz,keys=['TARGETID'],keep='last')
    logger.info(str(rann)+' length after cutting to unique TARGETID '+str(len(dz)))
    dz = join(dz,dzpd,keys=['TARGETID'],join_type='left')
    
    tin = np.isin(dz['TARGETID'],dzpd['TARGETID'])
    dz['NTILE'][~tin] = 0
    del dzpd
    logger.info(str(rann)+' length after joining to tiles info '+str(len(dz)))
    logger.info(str(rann)+' '+str(np.unique(dz['NTILE'])))

    if len(imbits) > 0:
        logger.info(str(rann)+' joining with original randoms to get mask properties')
        dirrt='/dvs_ro/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
        tcol = ['TARGETID','MASKBITS','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z'] #only including what are necessary for mask cuts for now
        #tcol = ['TARGETID','EBV','WISEMASK_W1','WISEMASK_W2','BRICKID','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G',\
        #'GALDEPTH_R','GALDEPTH_Z','PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','MASKBITS','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z']
        tarf = fitsio.read(dirrt+'/randoms-1-'+str(rann)+'.fits',columns=tcol)
        dz = join(dz,tarf,keys=['TARGETID'])
        logger.info(str(rann)+' completed join with original randoms to get mask properties')
        del tarf
        dz = common.cutphotmask(dz,imbits,logger=logger)
        logger.info(str(rann)+' length after cutting to based on imaging veto mask '+str(len(dz)))


    if 'PHOTSYS' not in cols:
        dz['PHOTSYS'] = 'N'
        sel = dz['DEC'] < 32.375
        wra = (dz['RA'] > 100-dz['DEC'])
        wra &= (dz['RA'] < 280 +dz['DEC'])
        sel |= ~wra
        dz['PHOTSYS'][sel] = 'S'


    common.write_LSS_scratchcp(dz,outf,logger=logger)
    #dz.write(outf,format='fits', overwrite=True)
    logger.info('wrote to '+outf)
    del dz


def mk_maskedran_wdup(gtl,indir,rann,imbits,outf,pd,ebits,notqso='',hpmapcut='_HPmapcut',ftiles=None,mapn=None,maps=None,mapcuts=None,reccircmasks=None):
    #apply the masks to the duplicated randoms associated with the data version
    #this will make mock processing more efficient, as it should only need to be done once and then used by all mock realizations
    #gtl should contain all hardware masking
    import LSS.common_tools as common
    #import logging
    logger = logging.getLogger('LSSran')

    zf = indir.replace('global','dvs_ro')+'/rancomb_'+str(rann)+pd+'wdupspec_zdone.fits'
    logger.info('about to load '+zf)
    dz = Table(fitsio.read(zf,columns=['LOCATION', 'TARGETID', 'RA', 'DEC', 'TILEID']))
    #logger.info('unique collision values '+str(np.unique(dz['COLLISION'])))
    logger.info('loaded '+zf+ ' '+str(len(dz))+' rows with columns '+str(dz.dtype.names))

    dz['TILELOCID'] = 10000*dz['TILEID'] +dz['LOCATION'] #reset it here in case was set by specdat and some matches were missing
    
    wg = np.isin(dz['TILELOCID'],gtl)
    dz = dz[wg]
    logger.info('now has '+str(len(dz))+' rows after masking bad hardware')
    
    #imaging veto mask
        
    logger.info(str(rann)+' getting mask info from original randoms ')
    dirrt='/dvs_ro/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
    tcol = ['TARGETID','MASKBITS','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z'] 
    tarf = fitsio.read(dirrt+'/randoms-1-'+str(rann)+'.fits',columns=tcol)
    keep = (tarf['NOBS_G']>0) & (tarf['NOBS_R']>0) & (tarf['NOBS_Z']>0)

    for biti in imbits:
        keep &= ((tarf['MASKBITS'] & 2**biti)==0)
    logger.info('from parent randoms for initial mask, '+str(np.sum(keep))+' kept out of '+str(len(tarf)))



    if isinstance(ebits, str):
        gtids = tarf['TARGETID'][keep]
        del tarf
        if 'lrg' in ebits:
            tracer_mask = 'lrg'
        mask_fn = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/main/LSS/randoms-1-'+str(rann)+tracer_mask+'imask.fits'
        tarf = fitsio.read(mask_fn)
        glrg = tarf['lrg_mask'] == 0
        glrgtids = tarf['TARGETID'][glrg]
        sel_lrg = np.isin(gtids,glrgtids)
        gtids = gtids[sel_lrg]
        logger.info('applying '+ebits+' leaves '+str(len(gtids)))
    else:
        for biti in ebits:
            keep &= ((tarf['MASKBITS'] & 2**biti)==0)
        gtids = tarf['TARGETID'][keep]
        del tarf
        logger.info('applying extra imaging mask leaves '+str(len(gtids)))

    sel_immask = np.isin(dz['TARGETID'],gtids)
    dz = dz[sel_immask]
    logger.info(str(len(dz))+ ' left in dup random file')

    if reccircmasks is not None:
        for maskfn in reccircmasks:
            mask = common.maskcircandrec(dz,maskfn,logger=logger)
            dz = dz[~mask]
            logger.info(str(len(dz))+ ' left in dup random file after reccircmasks')

    #healpix map veto
    if hpmapcut == '_HPmapcut':
        dz = common.apply_map_veto_arrays(dz,mapn,maps,mapcuts,nside=256,logger=logger)
    
    #add tile info
    if ftiles is None:
        logger.info('counting tiles from dz with columns '+str(dz.dtype.names))
        dzpd = count_tiles_input(np.copy(dz),logger=logger)#.keep_columns(['TARGETID','TILEID','TILELOCID']))
    else:
        dzpd = fitsio.read(ftiles)
    dz = join(dz,dzpd,keys=['TARGETID'],join_type='left')
    tin = np.isin(dz['TARGETID'],dzpd['TARGETID'])
    dz['NTILE'][~tin] = 0
    common.write_LSS_scratchcp(dz,outf,logger=logger)
    del dz
    return True    


def mkfullran(gtl,lznp,indir,rann,imbits,outf,tp,pd,notqso='',maxp=3400,min_tsnr2=0,tlid_full=None,badfib=None,ftiles=None,badfib_status=None):
    import LSS.common_tools as common
    #import logging
    logger = logging.getLogger('LSSran')

    if pd == 'bright':
        tscol = 'TSNR2_BGS'
    else:
        tscol = 'TSNR2_ELG'




    zf = indir.replace('global','dvs_ro')+'/rancomb_'+str(rann)+pd+'wdupspec_zdone.fits'
    logger.info('about to load '+zf)
    dz = Table.read(zf)
    logger.info('loaded '+zf+ ' with columns '+str(dz.dtype.names))
    #logger.info(dz.dtype.names)


    
    #dz = join(dz,dzpd,keys=['TARGETID'])
    #print('length including duplicates '+str(len(dz)))

    cols = list(dz.dtype.names)
    if tscol not in cols:
        dz[tscol] = np.ones(len(dz))


    #if 'TILELOCID' not in list(dz.dtype.names):
    dz['TILELOCID'] = 10000*dz['TILEID'] +dz['LOCATION'] #reset it here in case was set by specdat and some matches were missing
    #print(dz.dtype.names)
    
    
    if lznp is None:
        dz['ZPOSSLOC'] = np.ones(len(dz)).astype('bool')
    else:    
        dz['ZPOSSLOC'] = np.zeros(len(dz)).astype('bool')
        wk = ~np.isin(dz['TILELOCID'],lznp)
        dz['ZPOSSLOC'][wk] = 1

    wg = np.isin(dz['TILELOCID'],gtl)
    if badfib is not None:
        bad = np.isin(dz['FIBER'],badfib)
        logger.info('number at bad fibers '+str(sum(bad)))
        wg &= ~bad


    dz['GOODHARDLOC'] = np.zeros(len(dz)).astype('bool')
    dz['GOODHARDLOC'][wg] = 1
    if ftiles is None:
        logger.info('counting tiles from dz with columns '+str(dz.dtype.names))
        dzpd = count_tiles_input(dz[wg],logger=logger)#.keep_columns(['TARGETID','TILEID','TILELOCID']))
    else:
        dzpd = Table.read(ftiles)




    #zfpd = indir.replace('global','dvs_ro')+'/rancomb_'+str(rann)+pd+'_Alltilelocinfo.fits'
    #dzpd = Table.read(zfpd)


    #dzpd = count_tiles_input(np.array(dz[wg].keep_columns(['TARGETID','TILEID','TILELOCID']))

    #dz['LOCFULL'] = np.zeros(len(dz)).astype('bool')
    #if tlid_full is not None:
    #    wf = np.isin(dz['TILELOCID'],tlid_full)
    #    dz['LOCFULL'][wf] = 1

    logger.info(str(dz.dtype.names))
    dz['GOODPRI'] = np.zeros(len(dz)).astype('bool')
    sel = dz['PRIORITY'] <= maxp
    dz['GOODPRI'][sel] = 1
    #t0 = dz[tscol]*0 != 0
    #t0 |= dz[tscol] == 999999
    #t0 |= dz[tscol] == 1.e20
    #dz[tscol][t0] = 0
    
    #dz['GOODTSNR'] = np.zeros(len(dz)).astype('bool')
    #sel = dz[tscol] > min_tsnr2
    #dz['GOODTSNR'][sel] = 1
    #dz['sort'] =  dz['GOODPRI']*dz['GOODHARDLOC']*dz['ZPOSSLOC']*dz['GOODTSNR']*1+dz['GOODPRI']*dz['GOODHARDLOC']*dz['GOODTSNR']*1#-0.5*dz['LOCFULL']#*(1+dz[tsnr])
    dz['sort'] =  dz['GOODPRI']*dz['GOODHARDLOC']*dz['ZPOSSLOC']*1+dz['GOODPRI']*dz['GOODHARDLOC']*1
    #dz['sort'] =  dz['GOODPRI']*dz['GOODHARDLOC']*dz['ZPOSSLOC']#*(1+dz[tsnr])
    logger.info(dz.dtype.names)
    logger.info(str(rann)+' about to do sort')

    dz.sort('sort') #should allow to later cut on tsnr for match to data
    dz = unique(dz,keys=['TARGETID'],keep='last')
    logger.info(str(rann)+' length after cutting to unique TARGETID '+str(len(dz)))
    dz = join(dz,dzpd,keys=['TARGETID'],join_type='left')
    tin = np.isin(dz['TARGETID'],dzpd['TARGETID'])
    dz['NTILE'][~tin] = 0

    logger.info(str(rann)+' length after joining to tiles info '+str(len(dz)))
    logger.info(str(rann)+' '+str(np.unique(dz['NTILE'])))

    if len(imbits) > 0:
        logger.info(str(rann)+' joining with original randoms to get mask properties')
        dirrt='/dvs_ro/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
        tcol = ['TARGETID','MASKBITS','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z'] #only including what are necessary for mask cuts for now
        #tcol = ['TARGETID','EBV','WISEMASK_W1','WISEMASK_W2','BRICKID','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G',\
        #'GALDEPTH_R','GALDEPTH_Z','PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','MASKBITS','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z']
        tarf = fitsio.read(dirrt+'/randoms-1-'+str(rann)+'.fits',columns=tcol)
        dz = join(dz,tarf,keys=['TARGETID'])
        logger.info(str(rann)+' completed join with original randoms to get mask properties')
        del tarf
        dz = common.cutphotmask(dz,imbits)
        logger.info(str(rann)+' length after cutting to based on imaging veto mask '+str(len(dz)))


    if 'PHOTSYS' not in cols:
        dz['PHOTSYS'] = 'N'
        sel = dz['DEC'] < 32.375
        wra = (dz['RA'] > 100-dz['DEC'])
        wra &= (dz['RA'] < 280 +dz['DEC'])
        sel |= ~wra
        dz['PHOTSYS'][sel] = 'S'


    common.write_LSS_scratchcp(dz,outf,logger=logger)
    #dz.write(outf,format='fits', overwrite=True)
    logger.info('wrote to '+outf)
    del dz

def mkfullran_px(indir,rann,imbits,outf,tp,pd,gtl,lznp,px,dirrt,maxp=3400,min_tsnr2=0,tlid_full=None):
    import LSS.common_tools as common
    if pd == 'bright':
        tscol = 'TSNR2_BGS'
    else:
        tscol = 'TSNR2_ELG'

    zf = indir+'/rancomb_'+str(rann)+pd+'_'+str(px)+'_wdupspec_zdone.fits'
    #fe = False
    dz = []
    if os.path.isfile(zf):
        dz = Table(fitsio.read(zf))
        #dz.remove_columns(['TILES','NTILE'])
        wg = np.isin(dz['TILELOCID'],gtl)
        dz['GOODHARDLOC'] = np.zeros(len(dz)).astype('bool')
        dz['GOODHARDLOC'][wg] = 1
        #fe = True
        zfpd = indir+'/rancomb_'+str(rann)+pd+'_'+str(px)+'__Alltilelocinfo.fits'
        try:
            dzpd = Table(fitsio.read(zfpd))
        except:
            print('failed to load '+zfpd)
            return False

    if len(dz) > 0 and len(dzpd) > 0:# and fe:
        #dzpd.keep_columns(['TARGETID','TILES','NTILE'])
        dz = join(dz,dzpd,keys=['TARGETID'],join_type='left')
        #if maskzfail:
        #    wk = dz['ZPOSSNOTBAD'] == 1
        #else:
        #    wk = dz['ZPOSS'] == 1
        #print('length before cutting to good positions '+str(len(dz)))
        wk = ~np.isin(dz['TILELOCID'],lznp)
        dz['ZPOSSLOC'] = np.zeros(len(dz)).astype('bool')

        dz['ZPOSSLOC'][wk] = 1#dz[wk]
        #print('length after cutting to good positions '+str(len(dz)))
        dz['LOCFULL'] = np.zeros(len(dz)).astype('bool')
        if tlid_full is not None:
            wf = np.isin(dz['TILELOCID'],tlid_full)
            dz['LOCFULL'][wf] = 1


        if len(dz) > 0:

            tcol = ['TARGETID','MASKBITS','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z'] #only including what are necessary for mask cuts for now
            #tcol = ['TARGETID','EBV','WISEMASK_W1','WISEMASK_W2','BRICKID','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G',\
            #'GALDEPTH_R','GALDEPTH_Z','PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','MASKBITS','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z']
            tarfn = dirrt+'/randoms-1-hp-'+str(px)+'.fits'
            
            if os.path.isfile(tarfn):
                tarf = fitsio.read(tarfn,columns=tcol)
                dz = join(dz,tarf,keys=['TARGETID'])
                del tarf
                ranisf = True

                dz = common.cutphotmask(dz,imbits)
            else:
                print('did not find '+tarfn)
                ranisf = False
            #print('length after cutting to based on imaging veto mask '+str(len(dz)))
            if len(dz) > 0 and ranisf:
                #pl = np.copy(dz['PRIORITY']).astype(float)#dz['PRIORITY']
                #sp = pl <= 0
                #pl[sp] = .1
                #dz['sort'] = dz[tsnr]*dz['GOODHARDLOC']*dz['ZPOSSLOC']+dz['GOODHARDLOC']*dz['ZPOSSLOC']+dz['GOODHARDLOC']*dz['ZPOSSLOC']/pl#/dz['PRIORITY']
                dz['GOODPRI'] = np.zeros(len(dz)).astype('bool')
                sel = dz['PRIORITY'] <= maxp
                dz['GOODPRI'][sel] = 1
                t0 = dz[tscol]*0 != 0
                t0 |= dz[tscol] == 999999
                t0 |= dz[tscol] == 1.e20
                dz[tscol][t0] = 0
                dz['GOODTSNR'] = np.zeros(len(dz)).astype('bool')
                sel = dz[tscol] > min_tsnr2
                dz['GOODTSNR'][sel] = 1
                dz['sort'] =  dz['GOODPRI']*dz['GOODHARDLOC']*dz['ZPOSSLOC']*dz['GOODTSNR']*1#-0.5*dz['LOCFULL']#*(1+dz[tsnr])

                #dz['sort'] =  dz['GOODPRI']*dz['GOODHARDLOC']*dz['ZPOSSLOC']#*(1+dz[tsnr])


                dz.sort('sort') #should allow to later cut on tsnr for match to data
                dz = unique(dz,keys=['TARGETID'],keep='last')
                dz.remove_columns(['sort'])
                #print('length after cutting to unique TARGETID '+str(len(dz)))
                #print(np.unique(dz['NTILE']))
                dz.write(outf,format='fits', overwrite=True)
            else:
                print('0 rows left after imaging veto for '+outf+' so nothing got written')

        else:
            print('0 rows left for '+outf+' so nothing got written')
    else:
        print('no input file or redshift data before or after cutting to good obs for '+outf+' so nothing got written')
        if len(dz)>0:
            print('no entries in the tileloc file: #of targets, #goodloc, #in tileloc file')
            print(len(dz),np.sum(dz['GOODHARDLOC']),len(dzpd))
    del dz


def addcol_ran(fn,rann,dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/',ecol=['TARGETID','EBV','WISEMASK_W1','WISEMASK_W2','BRICKID','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G','GALDEPTH_R','GALDEPTH_Z','PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z']):
    dz = fitsio.read(fn)
    tarf = fitsio.read(dirrt+'/randoms-1-'+str(rann)+'.fits',columns=ecol)
    dz = join(dz,tarf,keys=['TARGETID'])
    dz.write(fn,format='fits', overwrite=True)
    del dz

def mkfulldat_mock(zf,imbits,ftar,tp,bit,outf,ftiles,maxp=3400,azf='',azfm='cumul',desitarg='DESI_TARGET',survey='Y1',specver='daily',notqso='',qsobit=4,min_tsnr2=0,badfib=None,gtl_all=None, mockz=None, mask_coll=False, mocknum=None, mockassigndir=None,logger=None):
    import LSS.common_tools as common
    """Make 'full' data catalog, contains all targets that were reachable, with columns denoted various vetos to apply
    ----------
    zf : :class:`str` path to the file containing merged potential targets and redshift 
        info
    imbits : :class:`list`, the list of imaging bits to mask against; ignored in None
        is passed
    ftar : :class:`~numpy.array` or`~astropy.table.Table`, contains extra target info
        to merge to. Ignored if None is passed.
    tp : :class:`str`, the target class
    bit : :class:`int`, the targeting bit corresponding to the targeting class and desitarg
        argument.
    outf : :class:`str`, path to write output to
    ftiles : :class:`str`, path to file containing information on how and where each target
    azf : :class:`str`, path to where to find extra redshift info for ELG/QSO catalogs
    azfm : :class:`str`, whether to use per tile ('cumul') or healpix redshifts ('hp')
    desitarg : :class:`str`, column to use when selecting on targeting bit
    specver : :class:`str`, version of spectroscopic reductions
    notqso : :class:`str`, if 'notqso', quasar targets are rejected
    qsobit : :class:`int`, targeting bit to select quasar targets
    min_tsnr2 : :class: `float`, minimum TSNR2_ value to cut on
    badfib : :class: `str`, path to list of bad fibers to cut agains
    Returns
    -------
    nothing
    Notes
    -----
    """


    if tp[:3] == 'BGS' or tp[:3] == 'MWS':
        pd = 'bright'
        tscol = 'TSNR2_BGS'
        #CHANGE TO HANDLE MOCK PATHS PROPERLY
        collf = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/collisions-BRIGHT.fits'
    else:
        pd = 'dark'
        tscol = 'TSNR2_ELG'
        collf = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/collisions-DARK.fits'

    if mockz and mask_coll:
        collf = mask_coll

    dz = Table(fitsio.read(zf))
    wtype = ((dz[desitarg] & bit) > 0)
    if notqso == 'notqso':
        common.printlog('removing QSO targets',logger)
        wtype &= ((dz[desitarg] & qsobit) == 0)

    #print(len(dz[wtype]))
    dz = dz[wtype]

    if mask_coll:
        coll = Table(fitsio.read(collf.replace('global','dvs_ro')))
        common.printlog('length before masking collisions '+str(len(dz)),logger)
        dz = setdiff(dz,coll,keys=['TARGETID','LOCATION','TILEID'])
        common.printlog('length after masking collisions '+str(len(dz)),logger)

    #instead of full spec data, we are going to get type specific data and cut to unique entries
    #in the end, we can only use the data associated with an observation
    #NOTE, this is not what we want to do for randoms, where instead we want to keep all of the
    #locations where it was possible a target could have been assigned
    #changing behavior back, load file that is spec info including zmtl
    specdir = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/'+specver+'/'
    prog = 'dark'
    if tp[:3] == 'BGS':
        prog = 'bright'

    specf = specdir+'datcomb_'+prog+'_spec_zdone.fits'
    #print(specf)
    fs = fitsio.read(specf.replace('global', 'dvs_ro'))
    fs = common.cut_specdat(fs,badfib,tsnr_min=min_tsnr2,tsnr_col=tscol)
    fs = common.mask_bad_petal_nights(fs, prog)
    fs = Table(fs)
    fs['TILELOCID'] = 10000*fs['TILEID'] +fs['LOCATION']
    gtl = np.unique(fs['TILELOCID'])
    #print('size of gtl', len(gtl))

    ''' FOR MOCKS with fiberassign, PUT IN SOMETHING TO READ FROM MOCK FIBERASSIGN INFO'''
    if mockz:
        assignf = os.path.join(mockassigndir, 'datcomb_{PROG}assignwdup.fits').format(PROG=prog)
        fs = fitsio.read(assignf.replace('global', 'dvs_ro'))
        fs = Table(fs)
        fs['TILELOCID'] = 10000*fs['TILEID'] +fs['LOCATION']


    fs.keep_columns(['TILELOCID','PRIORITY'])
    dz = join(dz, fs, keys=['TILELOCID'],join_type='left', uniq_col_name='{col_name}{table_name}',table_names=['','_ASSIGNED'])
    del fs
    dz['PRIORITY_ASSIGNED'] = dz['PRIORITY_ASSIGNED'].filled(999999)
    dz['GOODPRI'] = np.zeros(len(dz)).astype('bool')
    selp = dz['PRIORITY_ASSIGNED'] <= maxp
    selp |=  dz['PRIORITY_ASSIGNED'] == 999999
    dz['GOODPRI'][selp] = 1
    
    wg = np.isin(dz['TILELOCID'],gtl)
    common.printlog('Size of sample after cutting to gtl from data '+str(len(dz[wg])),logger)

    if gtl_all is not None:
        wg &= np.isin(dz['TILELOCID'],gtl_all)
    #print(len(dz[wg]))
    #print(len(dz[wg]))
    dz['GOODHARDLOC'] = np.zeros(len(dz)).astype('bool')
    dz['GOODHARDLOC'][wg] = 1
    common.printlog('length after selecting type '+str(len(dz)),logger)
    wz = dz['ZWARN'] != 999999 #this is what the null column becomes
    wz &= dz['ZWARN']*0 == 0 #just in case of nans
    dz['LOCATION_ASSIGNED'] = np.zeros(len(dz)).astype('bool')
    dz['LOCATION_ASSIGNED'][wz] = 1
    common.printlog('number assigned '+str(np.sum(dz['LOCATION_ASSIGNED'])),logger)
    common.printlog('number assigned at good priority '+str(np.sum(dz['LOCATION_ASSIGNED']*dz['GOODPRI']*1.)),logger)
    common.printlog('number assigned at good priority and good hardware '+str(np.sum(dz['LOCATION_ASSIGNED']*dz['GOODPRI']*dz['GOODHARDLOC']*1.)),logger)
    
    tlids = np.unique(dz['TILELOCID'][wz])
    wtl = np.isin(dz['TILELOCID'],tlids)
    dz['TILELOCID_ASSIGNED'] = np.zeros(len(dz)).astype('bool')
    dz['TILELOCID_ASSIGNED'][wtl] = 1
    common.printlog('number of unique targets at assigned tilelocid: '+str(len(np.unique(dz[wtl]['TARGETID']))),logger)
    #print(len(np.unique(dz[wtl]['TARGETID'])))

    cols = list(dz.dtype.names)
    #if tscol not in cols:
    #    dz[tscol] = np.ones(len(dz))
    #    print('added '+tscol+' and set all to 1') 


    #wnts = dz[tscol]*0 != 0
    #wnts |= dz[tscol] == 999999
    #dz[tscol][wnts] = 0
    #print(np.max(dz[tscol]))
    #dz['GOODTSNR'] = np.ones(len(dz)).astype('bool')
    #if min_tsnr2 > 0:
    #    sel = dz[tscol] > min_tsnr2
    #    dz['GOODTSNR'][sel] = 1
    
    if ftiles is None:
        dtl = count_tiles_input(dz[wg],logger=logger)
    else:
        dtl = Table.read(ftiles)

    #print('ECHO ',dz['TARGETID'][0])
    #df = dz.to_pandas()
    #df = df.sample(frac = 1).reset_index(drop=True)
    #dz = Table.from_pandas(df)

    #print('ECHO ',dz['TARGETID'][0])
    #if tp[:3] != 'QSO':
    
    ##AURE MAYBE KEEP OR NOT? SHOULD NOT HARM
    dz['ransort'] = np.random.random(len(dz))
    dz.sort('ransort')
    common.printlog('randomly sorted',logger)
    dz.remove_column('ransort')
    if tp[:3] == 'QSO':
        selnp = dz['LOCATION_ASSIGNED'] == 0
        pv = dz['PRIORITY'] #we will multiply by priority in order to keep priority 3400 over lya follow-up
        pv[selnp] = 0
        dz['sort'] = dz['LOCATION_ASSIGNED']*dz['GOODHARDLOC']*dz['GOODPRI']*pv+dz['TILELOCID_ASSIGNED']*dz['GOODHARDLOC']*dz['GOODPRI']*1  + dz['GOODHARDLOC']*1 + dz['GOODPRI']*1

    else:
        dz['sort'] = dz['LOCATION_ASSIGNED']*dz['GOODHARDLOC']*dz['GOODPRI']*1+dz['TILELOCID_ASSIGNED']*dz['GOODHARDLOC']*dz['GOODPRI']*1  + dz['GOODHARDLOC']*1 + dz['GOODPRI']*1

    dz.sort('sort')
    common.printlog('sorted',logger)
    
    dz = unique(dz,keys=['TARGETID'],keep='last')
    dz.remove_column('sort')
    common.printlog('cut number assigned '+str(np.sum(dz['LOCATION_ASSIGNED'])),logger)
    common.printlog('cut number assigned at good priority '+str(np.sum(dz['LOCATION_ASSIGNED']*dz['GOODPRI'])),logger)
    common.printlog('cut number assigned at good priority and good hardware '+str(np.sum(dz['LOCATION_ASSIGNED']*dz['GOODPRI']*dz['GOODHARDLOC'])),logger)

    #print('cut number assigned',np.sum(dz['LOCATION_ASSIGNED']))
    #print('cut number assigned at good priority',np.sum(dz['LOCATION_ASSIGNED']*dz['GOODPRI']))
    #print('cut number assigned at good priority and good hardwared',np.sum(dz['LOCATION_ASSIGNED']*dz['GOODPRI']*dz['GOODHARDLOC']))


    common.printlog('length after cutting to unique targets '+str(len(dz)),logger)
    #dtl = Table.read(ftiles)

    dtl.keep_columns(['TARGETID','NTILE','TILES'])#,'TILELOCIDS'])
    dz = join(dz,dtl,keys='TARGETID',join_type='left')
    tin = np.isin(dz['TARGETID'],dtl['TARGETID'])
    dz['NTILE'][~tin] = 0
    #print(np.unique(dz['NTILE']))
    if ftar is not None:
        print('joining to full imaging')
        remcol = ['RA','DEC','DESI_TARGET','BGS_TARGET']
    
        for col in remcol:
            if col in cols:
                dz.remove_columns([col]) #these come back in with merge to full target file
        dz = join(dz,ftar,keys=['TARGETID'])
    
    if specver == 'daily':
        spec_cols = ['TARGETID','TILEID','LOCATION','Z','ZERR','SPECTYPE','DELTACHI2'\
        ,'COADD_FIBERSTATUS','FIBERASSIGN_X','FIBERASSIGN_Y','COADD_NUMEXP','COADD_EXPTIME','COADD_NUMNIGHT'\
        ,'MEAN_DELTA_X','MEAN_DELTA_Y','RMS_DELTA_X','RMS_DELTA_Y','MEAN_PSF_TO_FIBER_SPECFLUX']
        dailydir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/'
        prog = 'dark'
        if tp[:3] == 'BGS':
            prog = 'bright'

        specdat = fitsio.read(dailydir+'datcomb_'+prog+'_spec_zdone.fits',columns=spec_cols)
        dz = join(dz,specdat,keys=['TARGETID','TILEID','LOCATION'],join_type='left')
    
    if len(imbits) > 0:
        dz = common.cutphotmask(dz,imbits)
        common.printlog('length after imaging mask; should not have changed '+str(len(dz)),logger)


    if tp[:3] == 'ELG' and azf != '' and azfm == 'cumul':# or tp == 'ELG_HIP':
        arz = Table(fitsio.read(azf,columns=['TARGETID','LOCATION','TILEID','OII_FLUX','OII_FLUX_IVAR']))
        arz['TILEID'] = arz['TILEID'].astype(int)
        dz = join(dz,arz,keys=['TARGETID','LOCATION','TILEID'],join_type='left')#,uniq_col_name='{col_name}{table_name}',table_names=['', '_OII'])
        o2c = np.log10(dz['OII_FLUX'] * np.sqrt(dz['OII_FLUX_IVAR']))+0.2*np.log10(dz['DELTACHI2'])
        w = (o2c*0) != 0
        w |= dz['OII_FLUX'] < 0
        o2c[w] = -20
        dz['o2c'] = o2c
        print('check length after merge with OII strength file:' +str(len(dz)))

    if tp[:3] == 'QSO' and azf != '' and azfm == 'cumul':
        arz = Table(fitsio.read(azf))
        arz.keep_columns(['TARGETID','LOCATION','TILEID','Z','Z_QN'])
        arz['TILEID'] = arz['TILEID'].astype(int)
        print(arz.dtype.names)
        dz = join(dz,arz,keys=['TARGETID','TILEID','LOCATION'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
        dz['Z'].name = 'Z_RR' #rename the original redrock redshifts
        dz['Z_QF'].name = 'Z' #the redshifts from the quasar file should be used instead

    if tp[:3] == 'ELG' and azf != '':
        print('number of masked oII row (hopefully matches number not assigned) '+ str(np.sum(dz['o2c'].mask)))
    if tp[:3] == 'QSO' and azf != '' and azfm == 'hp':
        arz = Table(fitsio.read(azf))
        sel = arz['SURVEY'] == 'main'
        sel &= arz['PROGRAM'] == 'dark'
        arz = arz[sel]
        arz.keep_columns(['TARGETID','Z','ZERR','Z_QN','TSNR2_LYA','TSNR2_QSO','QSO_MASKBITS'])
        
        print(arz.dtype.names)
        #arz['TILE'].name = 'TILEID'
        print('length of dz before QSO join '+str(len(dz)))
        dz = join(dz,arz,keys=['TARGETID'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
        print('length of dz after QSO join (shoudl be the same)'+str(len(dz)))
        dz['Z'].name = 'Z_RR' #rename the original redrock redshifts
        dz['Z_QF'].name = 'Z' #the redshifts from the quasar file should be used instead

    
    #needs to change because mocks actually need real spec info as well
    if mockz: #specver == 'mock':
        dz[mockz].name = 'Z' 
        
    if tp == 'QSO' and azf != '':
        print('number of good z according to qso file '+str(len(dz)-np.sum(dz['Z'].mask)))
    try:
        dz['Z'] = dz['Z'].filled(999999)
    except:
        common.printlog('filling masked Z rows did not succeed',logger)
    selm = dz['Z'] == 999999
    common.printlog('999999s for Z '+str(len(dz[selm])),logger)
    selo = dz['LOCATION_ASSIGNED'] == True
    #print('unique Z for unassigned:')
    #print(np.unique(dz[~selo]['Z']))

    common.printlog('length after cutting to unique targetid '+str(len(dz)),logger)
    #print('LOCATION_ASSIGNED numbers')
    #print(np.unique(dz['LOCATION_ASSIGNED'],return_counts=True))

    #print('TILELOCID_ASSIGNED numbers')
    #print(np.unique(dz['TILELOCID_ASSIGNED'],return_counts=True))

    probl = np.zeros(len(dz))

    #get completeness based on unique sets of tiles
    compa = []
    tll = []
    ti = 0
    common.printlog('getting completeness',logger)
    dz['TILES'] = dz['TILES'].filled('0')
    dz.sort('TILES')
    tlsl = dz['TILES']
    #tlsl.sort()
    nts = len(tlsl)
    
    tlslu = np.unique(tlsl)
    laa = dz['LOCATION_ASSIGNED']

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
            common.printlog('at tiles '+str(ti)+' of '+str(nts),logger)

        if nli == 0:
            common.printlog('no data for '+str(tlslu[ti]),logger)
            cp = 0
        else:
            cp = nai/nli#no/nt
        
        compa.append(cp)
        tll.append(tlslu[ti])
        ti += 1
    comp_dicta = dict(zip(tll, compa))
    fcompa = []
    for tl in dz['TILES']:
        fcompa.append(comp_dicta[tl])
    dz['COMP_TILE'] = np.array(fcompa)
    wc0 = dz['COMP_TILE'] == 0
    common.printlog('number of targets in 0 completeness regions '+str(len(dz[wc0])),logger)

    locl,nlocl = np.unique(dz['TILELOCID'],return_counts=True)
    wz = dz['LOCATION_ASSIGNED'] == 1
    dzz = dz[wz]

    loclz,nloclz = np.unique(dzz['TILELOCID'],return_counts=True)
    #print(np.max(nloclz),np.min(loclz))
    
    #print(len(locl),len(nloclz),sum(nlocl),sum(nloclz))
    natloc = ~np.isin(dz['TILELOCID'],loclz)
    common.printlog('number of unique targets around unassigned locations is '+str(np.sum(natloc)),logger)

    common.printlog('getting fraction assigned for each tilelocid',logger)
    nm = 0
    nmt =0
    pd = []
    nloclt = len(locl)
    lzs = np.isin(locl,loclz)
    for i in range(0,len(locl)):
        if i%100000 == 0:
            common.printlog('at row '+str(i)+' of '+str(nloclt),logger)
        nt = nlocl[i]
        nz = lzs[i]
        loc = locl[i]
        pd.append((loc,nz/nt))
    pd = dict(pd)
    for i in range(0,len(dz)):
        probl[i] = pd[dz['TILELOCID'][i]]
    
    #print('number of fibers with no observation, number targets on those fibers')
    #print(nm,nmt)

    #dz['FRACZ_TILELOCID'] = probl
    #print('sum of 1/FRACZ_TILELOCID, 1/COMP_TILE, and length of input; no longer rejecting unobserved loc, so wont match')
    #print(np.sum(1./dz[wz]['FRACZ_TILELOCID']),np.sum(1./dz[wz]['COMP_TILE']),len(dz))

    #print(np.unique(dz['NTILE']))

    common.printlog('number of fibers with no observation, number targets on those fibers: '+str(nm)+','+str(nmt),logger)
    #print(nm,nmt)

    dz['FRACZ_TILELOCID'] = probl
    common.printlog('sum of 1/FRACZ_TILELOCID, 1/COMP_TILE, and length of input; no longer rejecting unobserved loc, so wont match',logger)
    common.printlog(str(np.sum(1./dz[wz]['FRACZ_TILELOCID']))+','+str(np.sum(1./dz[wz]['COMP_TILE']))+','+str(len(dz)),logger)

    common.printlog('number of unique tileid: '+str(np.unique(dz['NTILE'])),logger)

    
    #needs to change, because specver should still point to real data
    if mockz:
        dz['PHOTSYS'] = 'N'
        sel = dz['DEC'] < 32.375
        wra = (dz['RA'] > 100-dz['DEC'])
        wra &= (dz['RA'] < 280 +dz['DEC'])
        sel |= ~wra
        dz['PHOTSYS'][sel] = 'S'
               
    
    common.write_LSS_scratchcp(dz,outf,logger=logger)
    #common.write_LSS(dz,outf)

def mkfulldat(zf,imbits,ftar,tp,bit,outf,ftiles,maxp=3400,azf='',azfm='cumul',emlin_fn=None,desitarg='DESI_TARGET',survey='Y1',specver='daily',notqso='',qsobit=4,min_tsnr2=0,badfib=None,badfib_status=None,gtl_all=None,mockz=None, mask_coll=False,logger=None, mocknum=None, mockassigndir=None,return_array='n',calc_ctile='y'):
    import LSS.common_tools as common
    """Make 'full' data catalog, contains all targets that were reachable, with columns denoted various vetos to apply
    ----------
    zf : :class:`str` path to the file containing merged potential targets and redshift 
        info
    imbits : :class:`list`, the list of imaging bits to mask against; ignored in None
        is passed
    ftar : :class:`~numpy.array` or`~astropy.table.Table`, contains extra target info
        to merge to. Ignored if None is passed.
    tp : :class:`str`, the target class
    bit : :class:`int`, the targeting bit corresponding to the targeting class and desitarg
        argument.
    outf : :class:`str`, path to write output to
    ftiles : :class:`str`, path to file containing information on how and where each target
    azf : :class:`str`, path to where to find extra redshift info for ELG/QSO catalogs
    azfm : :class:`str`, whether to use per tile ('cumul') or healpix redshifts ('hp')
    desitarg : :class:`str`, column to use when selecting on targeting bit
    specver : :class:`str`, version of spectroscopic reductions
    notqso : :class:`str`, if 'notqso', quasar targets are rejected
    qsobit : :class:`int`, targeting bit to select quasar targets
    min_tsnr2 : :class: `float`, minimum TSNR2_ value to cut on
    badfib : :class: `str`, path to list of bad fibers to cut agains
    Returns
    -------
    nothing
    Notes
    -----
    """


    if tp[:3] == 'BGS' or tp[:3] == 'MWS':
        pd = 'bright'
        tscol = 'TSNR2_BGS'
        #CHANGE TO HANDLE MOCK PATHS PROPERLY
        collf = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/collisions-BRIGHT.fits'
    else:
        pd = 'dark'
        tscol = 'TSNR2_ELG'
        collf = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/collisions-DARK.fits'

    
    if mockz and mask_coll:
        collf = mask_coll

    dz = Table(fitsio.read(zf))
    wtype = ((dz[desitarg] & bit) > 0)
    if notqso == 'notqso':
        if logger is not None:
            logger.info('removing QSO targets')
        else:
            print('removing QSO targets')
        wtype &= ((dz[desitarg] & qsobit) == 0)

    if logger is not None:
        logger.info('length before cut is '+str(len(dz))+', length of input cut after to type is '+str(len(dz[wtype])))
    else:
        print(len(dz[wtype]))
    dz = dz[wtype]

    if mask_coll:
        coll = Table(fitsio.read(collf))
        if logger is not None:
            logger.info('length before masking collisions '+str(len(dz)))
        else:
            print('length before masking collisions '+str(len(dz)))
        dz = setdiff(dz,coll,keys=['TARGETID','LOCATION','TILEID'])
        if logger is not None:
            logger.info('length after masking collisions '+str(len(dz)))
        else:        
            print('length after masking collisions '+str(len(dz)))

    #instead of full spec data, we are going to get type specific data and cut to unique entries
    #in the end, we can only use the data associated with an observation
    #NOTE, this is not what we want to do for randoms, where instead we want to keep all of the
    #locations where it was possible a target could have been assigned
    #changing behavior back, load file that is spec info including zmtl
    specdir = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/'+specver+'/'
    prog = 'dark'
    if tp[:3] == 'BGS':
        prog = 'bright'
    if 'LGE' in tp:
        prog = 'dark1b'


    if mockz:
        common.printlog('getting mock assignment info',logger)
        assignf = os.path.join(mockassigndir, 'datcomb_{PROG}assignwdup.fits').format(PROG=prog)
        fs = fitsio.read(assignf.replace('global', 'dvs_ro'))
        fs = Table(fs)
        fs['TILELOCID'] = 10000*fs['TILEID'] +fs['LOCATION']
    else:
        specf = specdir+'datcomb_'+prog+'_spec_zdone.fits'
        if logger is not None:
            logger.info('reading from spec file '+specf)
        else:
            print(specf)
        fs = fitsio.read(specf)
        #common.printlog('badfib type is '+str(type(badfib).__name__),logger)
        #common.printlog('badfib type row 0 '+str(type(badfib[0]).__name__),logger)
        if specver == 'daily':
            fs = common.cut_specdat(fs,badfib,tsnr_min=min_tsnr2,tsnr_col=tscol,fibstatusbits=badfib_status,remove_badfiber_spike_nz=False,mask_petal_nights=False,logger=logger)
        else:
            fs = common.cut_specdat(fs,badfib,tsnr_min=min_tsnr2,tsnr_col=tscol,fibstatusbits=badfib_status,remove_badfiber_spike_nz=True,mask_petal_nights=True,logger=logger)
        fs = Table(fs)
        fs['TILELOCID'] = 10000*fs['TILEID'] +fs['LOCATION']
        gtl = np.unique(fs['TILELOCID'])
    
    #print(len(gtl))
    fs.keep_columns(['TILELOCID','PRIORITY'])
    #''' FOR MOCKS with fiberassign, PUT IN SOMETHING TO READ FROM MOCK FIBERASSIGN INFO'''
    dz = join(dz,fs,keys=['TILELOCID'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_ASSIGNED'])
    if logger is not None:
        logger.info('columns after join to spec info '+str(dz.dtype.names))
    else:
        print(dz.dtype.names)
    del fs
    dz['PRIORITY_ASSIGNED'] = dz['PRIORITY_ASSIGNED'].filled(999999)
    dz['GOODPRI'] = np.zeros(len(dz)).astype('bool')
    selp = dz['PRIORITY_ASSIGNED'] <= maxp
    selp |=  dz['PRIORITY_ASSIGNED'] == 999999
    dz['GOODPRI'][selp] = 1
    
    if mockz:
        wg = np.ones(len(dz),dtype=bool)
        logger.info('good hardware all set to true because mock should already have been masked')
    else:
        #specf = specdir+'datcomb_'+prog+'_spec_zdone.fits'
        #if logger is not None:
        #    logger.info('reading from spec file '+specf)
        #else:
        #    print(specf)
    
        #fs = fitsio.read(specf)
        #fs = common.cut_specdat(fs,badfib,tsnr_min=min_tsnr2,tsnr_col=tscol,fibstatusbits=badfib_status,remove_badfiber_spike_nz=True,mask_petal_nights=True,logger=logger)
        #fs = Table(fs)
        #fs['TILELOCID'] = 10000*fs['TILEID'] +fs['LOCATION']
        #gtl = np.unique(fs['TILELOCID'])
    
        wg = np.isin(dz['TILELOCID'],gtl)
        #print(len(dz[wg]))
        if gtl_all is not None:
            wg &= np.isin(dz['TILELOCID'],gtl_all)
        if logger is not None:
            logger.info('number at good hardware '+str(len(dz[wg])))
        else:
            print(len(dz[wg]))
    #print(len(dz[wg]))
    dz['GOODHARDLOC'] = np.zeros(len(dz)).astype('bool')
    dz['GOODHARDLOC'][wg] = 1
    #print('length after selecting type '+str(len(dz)))

    wz = dz['ZWARN'] != 999999 #this is what the null column becomes
    wz &= dz['ZWARN']*0 == 0 #just in case of nans
    dz['LOCATION_ASSIGNED'] = np.zeros(len(dz)).astype('bool')
    dz['LOCATION_ASSIGNED'][wz] = 1
    if logger is not None:
        logger.info('number assigned '+str(np.sum(dz['LOCATION_ASSIGNED'])))
        logger.info('number assigned at good priority '+str(np.sum(dz['LOCATION_ASSIGNED']*dz['GOODPRI']*1.)))
        logger.info('number assigned at good priority and good hardware '+str(np.sum(dz['LOCATION_ASSIGNED']*dz['GOODPRI']*dz['GOODHARDLOC']*1.)))

    else:
        print('number assigned',np.sum(dz['LOCATION_ASSIGNED']))
        print('number assigned at good priority',np.sum(dz['LOCATION_ASSIGNED']*dz['GOODPRI']*1.))
        print('number assigned at good priority and good hardware',np.sum(dz['LOCATION_ASSIGNED']*dz['GOODPRI']*dz['GOODHARDLOC']*1.))
    tlids = np.unique(dz['TILELOCID'][wz])
    wtl = np.isin(dz['TILELOCID'],tlids)
    dz['TILELOCID_ASSIGNED'] = np.zeros(len(dz)).astype('bool')
    dz['TILELOCID_ASSIGNED'][wtl] = 1
    if logger is not None:
        logger.info('number of unique targets at assigned tilelocid: '+str(len(np.unique(dz[wtl]['TARGETID']))) )
    else:
        print('number of unique targets at assigned tilelocid:')
        print(len(np.unique(dz[wtl]['TARGETID'])))

    cols = list(dz.dtype.names)
    #if tscol not in cols:
    #    dz[tscol] = np.ones(len(dz))
    #    print('added '+tscol+' and set all to 1') 


    #wnts = dz[tscol]*0 != 0
    #wnts |= dz[tscol] == 999999
    #dz[tscol][wnts] = 0
    #print(np.max(dz[tscol]))
    #dz['GOODTSNR'] = np.ones(len(dz)).astype('bool')
    #if min_tsnr2 > 0:
    #    sel = dz[tscol] > min_tsnr2
    #    dz['GOODTSNR'][sel] = 1
    
    common.printlog('getting tile counts',logger)
    if ftiles is None:
        dtl = count_tiles_input(dz[wg],logger=logger)
        #dtl = count_tiles_input_alt(dz[wg],logger=logger) #alt was faster when run in notebook but is much slower on test on node...
    else:
        dtl = Table.read(ftiles)
    
    #if tp[:3] != 'QSO':
    if tp[:3] == 'QSO':
        selnp = dz['LOCATION_ASSIGNED'] == 0
        pv = dz['PRIORITY'] #we will multiply by priority in order to keep priority 3400 over lya follow-up
        pv[selnp] = 0
        #dz['sort'] = dz['LOCATION_ASSIGNED']*dz['GOODTSNR']*dz['GOODHARDLOC']*dz['GOODPRI']*pv+dz['TILELOCID_ASSIGNED']*dz['GOODHARDLOC']*dz['GOODPRI']*1  + dz['GOODHARDLOC']*1 + dz['GOODPRI']*1#*(1+np.clip(dz[tscol],0,200))*1+dz['TILELOCID_ASSIGNED']*dz['GOODHARDLOC']*1+dz['GOODHARDLOC']*1
        dz['sort'] = dz['LOCATION_ASSIGNED']*dz['GOODHARDLOC']*dz['GOODPRI']*pv+dz['TILELOCID_ASSIGNED']*dz['GOODHARDLOC']*dz['GOODPRI']*1  + dz['GOODHARDLOC']*1 + dz['GOODPRI']*1#
    else:
        #dz['sort'] = dz['LOCATION_ASSIGNED']*dz['GOODTSNR']*dz['GOODHARDLOC']*dz['GOODPRI']*1+dz['TILELOCID_ASSIGNED']*dz['GOODHARDLOC']*dz['GOODPRI']*1  + dz['GOODHARDLOC']*1 + dz['GOODPRI']*1#*(1+np.clip(dz[tscol],0,200))*1+dz['TILELOCID_ASSIGNED']*dz['GOODHARDLOC']*1+dz['GOODHARDLOC']*1
        dz['sort'] = dz['LOCATION_ASSIGNED']*dz['GOODHARDLOC']*dz['GOODPRI']*1+dz['TILELOCID_ASSIGNED']*dz['GOODHARDLOC']*dz['GOODPRI']*1  + dz['GOODHARDLOC']*1 + dz['GOODPRI']*1
    #else:
    #    selnp = dz['LOCATION_ASSIGNED'] == 0
    #    pv = dz['PRIORITY']
    #    pv[selnp] = 0
    #    dz['sort'] = dz['LOCATION_ASSIGNED']*dz['GOODTSNR']*dz['GOODHARDLOC']*1+dz['TILELOCID_ASSIGNED']*dz['GOODHARDLOC']*1+dz['GOODHARDLOC']*1/(dz['PRIORITY_ASSIGNED']+2)
    if logger is not None:
        logger.info('about to sort')
    else:
        print('about to sort')

    dz.sort('sort')
    if logger is not None:
        logger.info('sorted')
    else:
        print('sorted')
    
    dz = unique(dz,keys=['TARGETID'],keep='last')
    common.printlog('cut to unique targetid',logger)
    dz.remove_column('sort')
    
    if logger is not None:
        logger.info('cut number assigned '+str(np.sum(dz['LOCATION_ASSIGNED'])))
        logger.info('cut number assigned at good priority '+str(np.sum(dz['LOCATION_ASSIGNED']*dz['GOODPRI'])))
        logger.info('cut number assigned at good priority and good hardware '+str(np.sum(dz['LOCATION_ASSIGNED']*dz['GOODPRI']*dz['GOODHARDLOC'])))

    else:
        print('cut number assigned',np.sum(dz['LOCATION_ASSIGNED']))
        print('cut number assigned at good priority',np.sum(dz['LOCATION_ASSIGNED']*dz['GOODPRI']))
        print('cut number assigned at good priority and good hardwared',np.sum(dz['LOCATION_ASSIGNED']*dz['GOODPRI']*dz['GOODHARDLOC']))


    if logger is not None:
        logger.info('length after cutting to unique targets '+str(len(dz)))
    else:
        print('length after cutting to unique targets '+str(len(dz)))
    #dtl = Table.read(ftiles)

    common.printlog('joining to ntile info',logger=logger)
    dtl.keep_columns(['TARGETID','NTILE','TILES'])#,'TILELOCIDS'])
    dz = join(dz,dtl,keys='TARGETID',join_type='left')
    tin = np.isin(dz['TARGETID'],dtl['TARGETID'])
    dz['NTILE'][~tin] = 0
    #print(np.unique(dz['NTILE']))
    if ftar is not None:
        common.printlog('joining to full imaging',logger)
        remcol = ['RA','DEC','DESI_TARGET','BGS_TARGET']
    
        for col in remcol:
            if col in cols:
                dz.remove_columns([col]) #these come back in with merge to full target file
        dz = join(dz,ftar,keys=['TARGETID'])
    
    if specver == 'daily':
        spec_cols = ['TARGETID','TILEID','LOCATION','Z','ZERR','SPECTYPE','DELTACHI2'\
        ,'COADD_FIBERSTATUS','FIBERASSIGN_X','FIBERASSIGN_Y','COADD_NUMEXP','COADD_EXPTIME','COADD_NUMNIGHT'\
        ,'MEAN_DELTA_X','MEAN_DELTA_Y','RMS_DELTA_X','RMS_DELTA_Y','MEAN_PSF_TO_FIBER_SPECFLUX']
        dailydir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/'
        prog = 'dark'
        if tp[:3] == 'BGS':
            prog = 'bright'
        if tp == 'LGE':
            prog = 'dark1b'
        specdat = fitsio.read(dailydir+'datcomb_'+prog+'_spec_zdone.fits',columns=spec_cols)
        dz = join(dz,specdat,keys=['TARGETID','TILEID','LOCATION'],join_type='left')
    
    if len(imbits) > 0:
        dz = common.cutphotmask(dz,imbits,logger=logger)
        if logger is not None:
            logger.info('length after imaging mask; should not have changed '+str(len(dz)))
        else:
            print('length after imaging mask; should not have changed '+str(len(dz)))


    if tp[:3] == 'ELG' and azf != '' and azfm == 'cumul':# or tp == 'ELG_HIP':
        arz = Table(fitsio.read(azf,columns=['TARGETID','LOCATION','TILEID','OII_FLUX','OII_FLUX_IVAR']))
        arz['TILEID'] = arz['TILEID'].astype(int)
        dz = join(dz,arz,keys=['TARGETID','LOCATION','TILEID'],join_type='left')#,uniq_col_name='{col_name}{table_name}',table_names=['', '_OII'])
        o2c = np.log10(dz['OII_FLUX'] * np.sqrt(dz['OII_FLUX_IVAR']))+0.2*np.log10(dz['DELTACHI2'])
        w = (o2c*0) != 0
        w |= dz['OII_FLUX'] < 0
        o2c[w] = -20
        dz['o2c'] = o2c
        if logger is not None:
            logger.info('check length after merge with OII strength file:' +str(len(dz)))
        else:
            print('check length after merge with OII strength file:' +str(len(dz)))

    if tp[:3] == 'QSO' and azf != '' and azfm == 'cumul':
        arz = Table(fitsio.read(azf))
        arz.keep_columns(['TARGETID','LOCATION','TILEID','Z','Z_QN'])
        arz['TILEID'] = arz['TILEID'].astype(int)
        #print(arz.dtype.names)
        dz = join(dz,arz,keys=['TARGETID','TILEID','LOCATION'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
        dz['Z'].name = 'Z_RR' #rename the original redrock redshifts
        dz['Z_QF'].name = 'Z' #the redshifts from the quasar file should be used instead
        if emlin_fn is not None:
            emcat =  Table(fitsio.read(emlin_fn,columns=['TARGETID','LOCATION','TILEID','OII_FLUX','OII_FLUX_IVAR','OIII_FLUX','OIII_FLUX_IVAR']))
            emcat['TILEID'] = emcat['TILEID'].astype(int)
            dz = join(dz,emcat,keys=['TARGETID','LOCATION','TILEID'],join_type='left')

    if tp[:3] == 'ELG' and azf != '':
        if logger is not None:
            logger.info('number of masked oII row (hopefully matches number not assigned) '+ str(np.sum(dz['o2c'].mask)))
        else:
            print('number of masked oII row (hopefully matches number not assigned) '+ str(np.sum(dz['o2c'].mask)))
    if tp[:3] == 'QSO' and azf != '' and azfm == 'hp':
        message = 'adding healpix based QSO info'
        common.printlog(message,logger)
        arz = Table(fitsio.read(azf))
        sel = arz['SURVEY'] == 'main'
        sel &= arz['PROGRAM'] == 'dark'
        arz = arz[sel]
        arz.keep_columns(['TARGETID','Z','ZERR','Z_QN','TSNR2_LYA','TSNR2_QSO','QSO_MASKBITS'])
        
        #print(arz.dtype.names)
        #arz['TILE'].name = 'TILEID'
        common.printlog('length of dz before QSO join '+str(len(dz)),logger)
        dz = join(dz,arz,keys=['TARGETID'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
        common.printlog('length of dz after QSO join (shoudl be the same)'+str(len(dz)),logger)
        dz['Z'].name = 'Z_RR' #rename the original redrock redshifts
        dz['Z_QF'].name = 'Z' #the redshifts from the quasar file should be used instead

    
    #needs to change because mocks actually need real spec info as well
    if mockz: #specver == 'mock':
        dz[mockz].name = 'Z' 

    #if specver == 'mock':
    #    dz[mockz].name = 'Z' 
        
    if tp == 'QSO' and azf != '':
        common.printlog('number of good z according to qso file '+str(len(dz)-np.sum(dz['Z'].mask)),logger)
    try:
    #if dz.masked:
        common.printlog('filling masked Z rows with 999999',logger)
        dz['Z'] = dz['Z'].filled(999999)
    #else:
        common.printlog('table is not masked, no masked rows to fill',logger)
    except:
        common.printlog('filling masked Z rows did not succeed, perhaps it is not masked',logger)
    selm = dz['Z'] == 999999
    common.printlog('999999s for Z '+str(len(dz[selm])),logger)
    selo = dz['LOCATION_ASSIGNED'] == True
    common.printlog('unique Z for unassigned: '+str(np.unique(dz[~selo]['Z'])),logger)

    common.printlog('length after cutting to unique targetid '+str(len(dz)),logger)
    common.printlog('LOCATION_ASSIGNED numbers',logger)
    common.printlog(str(np.unique(dz['LOCATION_ASSIGNED'],return_counts=True)),logger)

    common.printlog('TILELOCID_ASSIGNED numbers',logger)
    common.printlog(str(np.unique(dz['TILELOCID_ASSIGNED'],return_counts=True)),logger)

    probl = np.zeros(len(dz))

    #get completeness based on unique sets of tiles
    compa = []
    tll = []
    ti = 0
    common.printlog('getting completeness',logger)
    #if dz['TILES'].masked:
    try:
        dz['TILES'] = dz['TILES'].filled('0')
        common.printlog('filling masked TILES values',logger)
    except:
        common.printlog('filling masked TILES values did not succeed, perhaps it is not masked',logger)
    #dz.sort('TILES')
    tlsl = np.array(dz['TILES'])
    #common.printlog(str(tlsl.dtype),logger)
    #tlsl.sort()
    nts = len(tlsl)
    

    if calc_ctile == 'y':
        tlslu,indices,cnts= np.unique(tlsl,return_inverse=True,return_counts=True)
        n_of_tiles = len(tlslu)
        laa = dz['LOCATION_ASSIGNED']
        acnts = np.bincount(indices,laa)
        compa = acnts/cnts
        #i = 0
        #while i < len(dz):
        #    tls  = []
#             tlis = []
#             nli = 0
#             nai = 0
#     
#             while tlsl[i] == tlslu[ti]:
#                 nli += 1
#                 nai += laa[i]
#                 i += 1
#                 if i == len(dz):
#                     break
#     
#             if ti%1000 == 0:
#                 common.printlog('at tiles '+str(ti)+' of '+str(n_of_tiles),logger)
#     
#             if nli == 0:
#                 common.printlog('no data for '+str(tlslu[ti])+' and '+str(tlsl[i]),logger)
#                 cp = 0
#                 i += 1
#             else:
#                 cp = nai/nli#no/nt
            
#            compa.append(cp)
#            tll.append(tlslu[ti])
#            ti += 1
        comp_dicta = dict(zip(tlslu, compa))
        fcompa = []
        for tl in dz['TILES']:
            fcompa.append(comp_dicta[tl])
        dz['COMP_TILE'] = np.array(fcompa)
        wc0 = dz['COMP_TILE'] == 0
        common.printlog('number of targets in 0 completeness regions '+str(len(dz[wc0])),logger)
    else:
        dz['COMP_TILE'] = 1

    locl,nlocl = np.unique(dz['TILELOCID'],return_counts=True)
    wz = dz['LOCATION_ASSIGNED'] == 1
    dzz = dz[wz]

    loclz,nloclz = np.unique(dzz['TILELOCID'],return_counts=True)
    #print(np.max(nloclz),np.min(loclz))
    
    #print(len(locl),len(nloclz),sum(nlocl),sum(nloclz))
    natloc = ~np.isin(dz['TILELOCID'],loclz)
    common.printlog('number of unique targets around unassigned locations is '+str(np.sum(natloc)),logger)

    common.printlog('getting fraction assigned for each tilelocid',logger)
    nm = 0
    nmt =0
    pd = []
    nloclt = len(locl)
    lzs = np.isin(locl,loclz)
    for i in range(0,len(locl)):
        if i%100000 == 0:
            common.printlog('at row '+str(i)+' of '+str(nloclt),logger)
        nt = nlocl[i]
        nz = lzs[i]
        loc = locl[i]
        pd.append((loc,nz/nt))
    pd = dict(pd)
    for i in range(0,len(dz)):
        probl[i] = pd[dz['TILELOCID'][i]]
    common.printlog('number of fibers with no observation, number targets on those fibers: '+str(nm)+','+str(nmt),logger)
    #print(nm,nmt)

    dz['FRACZ_TILELOCID'] = probl
    common.printlog('sum of 1/FRACZ_TILELOCID, 1/COMP_TILE, and length of input; no longer rejecting unobserved loc, so wont match',logger)
    common.printlog(str(np.sum(1./dz[wz]['FRACZ_TILELOCID']))+','+str(np.sum(1./dz[wz]['COMP_TILE']))+','+str(len(dz)),logger)

    common.printlog('number of unique tileid: '+str(np.unique(dz['NTILE'])),logger)
    
    #needs to change, because specver should still point to real data
    if mockz:
        dz['PHOTSYS'] = 'N'
        sel = dz['DEC'] < 32.375
        wra = (dz['RA'] > 100-dz['DEC'])
        wra &= (dz['RA'] < 280 +dz['DEC'])
        sel |= ~wra
        dz['PHOTSYS'][sel] = 'S'
               
    
    common.write_LSS_scratchcp(dz,outf,logger=logger)
    if return_array == 'y':
        return dz


def get_ELG_SSR_tile(ff,o2c_thresh,zmin=.6,zmax=1.5,tsnrcut=80):
    ff['relSSR_tile'] = np.zeros(len(ff))
    wo = ff['ZWARN']*0 == 0
    wo &= ff['ZWARN'] != 999999
    wo &= ff['LOCATION_ASSIGNED'] == 1
    wo &= ff['TSNR2_ELG'] > tsnrcut
    wno = ff['SPECTYPE'] == 'QSO'
    wno &= ff['SPECTYPE'] == 'STAR'
    wno &= ((ff['ZWARN'] == 0) & (ff['Z']<0.6))
    fall = ff[wo&~wno]
    wg = ff['o2c'] > o2c_thresh
    ssr_all = len(ff[wo&~wno&wg])/len(fall)
    print('overall success rate: '+str(ssr_all))

    tids = np.unique(ff['TILEID'])
    for tid in tids:
        selt = ff['TILEID'] == tid
        ssr_t = len(ff[wo&~wno&wg&selt])/len(ff[wo&~wno&selt])
        ff['relSSR_tile'][selt] = ssr_t/ssr_all
        print(tid,ssr_t,ssr_t/ssr_all)
    return ff

def add_zfail_weight2fullQSO(indir,version,qsocat,tsnrcut=80,readpars=False,logger=None):
    import LSS.common_tools as common
    from LSS import ssr_tools_new

    ff = fitsio.read(indir+'datcomb_QSO_tarspecwdup_zdone.fits')
    selobs = ff['ZWARN'] != 999999
    selobs &= ff['TSNR2_ELG'] > tsnrcut
    ff = ff[selobs]
    azf = qsocat
    arz = Table(fitsio.read(azf))
    arz.keep_columns(['TARGETID','LOCATION','TILEID','Z','Z_QN'])
    arz['TILEID'] = arz['TILEID'].astype(int)
    print(arz.dtype.names)
    ff = join(ff,arz,keys=['TARGETID','TILEID','LOCATION'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
    ff['Z'].name = 'Z_RR' #rename the original redrock redshifts
    ff['Z_QF'].name = 'Z_not4clus' #the redshifts from the quasar file should be used instead
    ff = common.addNS(ff)
    ff = common.cut_specdat(ff)
    outdir = indir+'LSScats/'+version+'/'
    tp = 'QSO'
    ffv = Table.read(outdir+tp+'_full_noveto.dat.fits')
    ffv.keep_columns(['TARGETID','FIBERFLUX_R','MW_TRANSMISSION_R','EBV'])
    ff = join(ff,ffv,keys=['TARGETID'],join_type='left')
    print(ff.dtype.names)
    mintsnr=450/(8.60/0.255)
    maxtsnr=1800/(8.60/0.255)
    band = 'R'
    
    
    selp = ff['PRIORITY'] == 3400

    s = 0
    modl =[]
    regl = ['S','N']
    for reg in regl:
        # success rate dips at TSNR > 40 in NGC, probably due to chance alignment with low galactic latitude
        # and high contamination.
        # removing TSNR > 40 from the fits makes very little difference, however (probably because
        # the model cannot model a decreasing success rate with increasing TSNR, so the dip just
        # doesn't change anything)
        #if reg == 'N':
        #    maxtsnr = 40.
        #else:
        #    maxtsnr = 1800/(8.60/0.255)
        mod = ssr_tools_new.model_ssr(ff[selp],tsnr_min=mintsnr,tsnr_max=maxtsnr,tracer='QSO',reg=reg,outdir=outdir,band=band,outfn_root='QSO',readpars=readpars)
        modl.append(mod)    

    ff = Table.read(outdir+tp+'_full_noveto.dat.fits')
    wzf = np.ones(len(ff))
    msr = np.ones(len(ff))
    selobs = ff['ZWARN']*0 == 0
    selobs &= ff['ZWARN'] != 999999
    #selobs &= ff['GOODHARDLOC'] == 1
    selobs &= ff['TSNR2_'+tp[:3]]*0 == 0
    selp = ff['PRIORITY'] == 3400
    selgz = common.goodz_infull(tp[:3],ff,zcol='Z')
    
    print('check that ~98% fulfill priority cut:')
    print(np.sum(selp&selobs)/np.sum(selobs))
    
    
    for reg,mod in zip(regl,modl):
        selreg = ff['PHOTSYS'] == reg
        wts,md = mod.add_modpre(ff[selobs&selreg&selp])
        
        wzf[selobs&selreg&selp] = wts
        msr[selobs&selreg&selp] = md
        
        print('compare good z frac to sum of model')
        print(len(ff[selgz&selobs&selreg])/len(ff[selobs&selreg]),np.sum(msr[selobs&selreg])/len(ff[selobs&selreg]))

    ff['WEIGHT_ZFAIL'] = wzf
    ff['mod_success_rate'] = msr
    
    #print(len(ff[wz]),len(ff))

    print('min/max of zfail weights:')
    print(np.min(ff[selobs]['WEIGHT_ZFAIL']),np.max(ff[selobs]['WEIGHT_ZFAIL']))
 
        


    #plt.plot(ff[selgz&selobs]['TSNR2_'+tp[:3]],ff[selgz&selobs]['WEIGHT_ZFAIL'],'k,')
    #plt.xlim(np.percentile(ff[selgz]['TSNR2_'+tp[:3]],0.5),np.percentile(ff[selgz]['TSNR2_'+tp[:3]],99))
    #plt.show()
    
    common.write_LSS_scratchcp(ff,outdir+tp+'_full_noveto.dat.fits',logger=logger)
    ff.keep_columns(['TARGETID','WEIGHT_ZFAIL','mod_success_rate'])
    ffc = Table.read(outdir+tp+'_full.dat.fits')
    cols = list(ffc.dtype.names)
    if 'WEIGHT_ZFAIL' in cols:
        ffc.remove_columns(['WEIGHT_ZFAIL'])
    if 'mod_success_rate' in cols:
        ffc.remove_columns(['mod_success_rate'])
    ffc = join(ffc,ff,keys=['TARGETID'],join_type='left')
    common.write_LSS_scratchcp(ffc,outdir+tp+'_full.dat.fits',comments='added ZFAIL weight',logger=logger)

    fname_mapveto = outdir+tp+'_full_HPmapcut.dat.fits'
    if os.path.isfile(fname_mapveto):
        ff.keep_columns(['TARGETID','WEIGHT_ZFAIL','mod_success_rate'])
        ffc = Table.read(fname_mapveto)
        cols = list(ffc.dtype.names)
        if 'WEIGHT_ZFAIL' in cols:
            ffc.remove_columns(['WEIGHT_ZFAIL'])
        if 'mod_success_rate' in cols:
            ffc.remove_columns(['mod_success_rate'])
        ffc = join(ffc,ff,keys=['TARGETID'],join_type='left')
        common.write_LSS_scratchcp(ffc,fname_mapveto,logger=logger)#,comments='added ZFAIL weight')
 


def add_zfail_weight2full(indir,tp='',tsnrcut=80,readpars=False,hpmapcut='_HPmapcut',logger=None):
    import LSS.common_tools as common
    from LSS import ssr_tools_new
    '''
    fl is the root of the input/output file
    weighttileloc determines whether to include 1/FRACZ_TILELOCID as a completeness weight
    zmask determines whether to apply a mask at some given redshift
    tp is the target type
    dchi2 is the threshold for keeping as a good redshift
    tnsrcut determines where to mask based on the tsnr2 value (defined below per tracer)

    '''
    
    ff = Table.read(indir+tp+'_full'+hpmapcut+'.dat.fits')
    cols = list(ff.dtype.names)
    if 'Z' in cols:
        #print('Z column already in full file')
    #else:
        #ff['Z_not4clus'].name = 'Z'
        ff['Z'].name = 'Z_not4clus'
        common.write_LSS(ff,indir+tp+'_full'+hpmapcut+'.dat.fits',comments='changed Z column back to Z_not4clus')

    #selobs = ff['ZWARN'] == 0
    selobs = ff['ZWARN']*0 == 0
    selobs &= ff['ZWARN'] != 999999
    
    regl = ['S','N']

    
    if tp == 'QSO':
        selobs &= ff['TSNR2_ELG'] > tsnrcut
        #might want to add something to remove stars/galaxies from selobs
        mintsnr=450/(8.60/0.255)
        maxtsnr=1800/(8.60/0.255)
        band = 'R'
        
    if tp[:3] == 'ELG':
        selobs &= ff['TSNR2_ELG'] > tsnrcut
        mintsnr = 80
        maxtsnr = 200
        band = 'G'

    if tp == 'LRG':
        selobs &= ff['TSNR2_ELG'] > tsnrcut
        mintsnr=500/12.15
        #maxtsnr =2000/12.15
        maxtsnr =1700/12.15
        band = 'Z'
        
    if tp[:3] == 'BGS':
        selobs &= ff['TSNR2_BGS'] > tsnrcut
        mintsnr=120/(12.15/89.8)
        maxtsnr =300/(12.15/89.8)

        band = 'R'


    #ffz = ff[wz]
    #print('length after cutting to good z '+str(len(ff[wz])))
    #ff['GOODZ'] = wz
    #if tp != 'LRG':
    #    ff['WEIGHT_ZFAIL'] = np.ones(len(ff))
    #    ff['mod_success_rate'] = np.ones(len(ff))
    #selobs = ff['ZWARN'] != 999999
    s = 0
    modl =[]
    for reg in regl:
        #selreg = np.ones(len(ff),dtype='bool')
        #if reg is not None:
        #    print('working with data from region '+reg)
        #    gal = func(surveys=[survey],specrels=[specrel],versions=[version],efftime_min=minefftime,efftime_max=maxefftime,reg=reg)
        #    selreg = ff['PHOTSYS'] == reg            
        #else:    
        #    print('working with the full data, no region split')
        #    gal = func(surveys=[survey],specrels=[specrel],versions=[version],efftime_min=minefftime,efftime_max=maxefftime)
            
        if tp[:3] == 'ELG':
            mod = ssr_tools_new.model_ssr_zfac(ff[selobs],reg=reg,outdir=indir)
        else:
            mod = ssr_tools_new.model_ssr(ff[selobs],tsnr_min=mintsnr,tsnr_max=maxtsnr,tracer=tp[:3],reg=reg,outdir=indir,band=band,outfn_root=tp,readpars=readpars)
        modl.append(mod)
        #ffwz = gal.add_modpre(ff[selobs&selreg])
        #print(min(ffwz['mod_success_rate']),max(ffwz['mod_success_rate']))
        #print(min(ffwz['WEIGHT_ZFAIL']),max(ffwz['WEIGHT_ZFAIL']))
        #ffwz['WEIGHT_ZFAIL'] = 1./ffwz['mod_success_rate']
        #ffwz.keep_columns(['TARGETID','WEIGHT_ZFAIL','mod_success_rate'])
        #if s == 0:
        #    rem_cols = ['WEIGHT_ZFAIL','mod_success_rate']
        #    for col in rem_cols:
        #        try:
        #            ff.remove_columns([col])
        #            print(col +' was in full file and will be replaced')
        #        except:
        #            print(col +' was not yet in full file')    
        #ff = join(ff,ffwz,keys=['TARGETID'],join_type='left')
        #wzf[selobs&selreg] = np.array(ffwz['WEIGHT_ZFAIL'])
        #msr[selobs&selreg] =  np.array(ffwz['mod_success_rate'])
        #print(min(wzf),max(wzf))
        #s = 1
    common.printlog('model obtained',logger)
    if tp == 'BGS_BRIGHT-21.5':
        fullname = indir+tp+'_full'+hpmapcut+'.dat.fits'
    else:
        fullname = indir+tp+'_full_noveto.dat.fits'
    ff = Table.read(fullname)
    wzf = np.ones(len(ff))
    msr = np.ones(len(ff))
    selobs = ff['ZWARN']*0 == 0
    selobs &= ff['ZWARN'] != 999999
    #selobs &= ff['GOODHARDLOC'] == 1
    selobs &= ff['TSNR2_'+tp[:3]]*0 == 0
    selgz = common.goodz_infull(tp[:3],ff,zcol='Z')
    for reg,mod in zip(regl,modl):
        selreg = ff['PHOTSYS'] == reg
        wts,md = mod.add_modpre(ff[selobs&selreg])
        wzf[selobs&selreg] = wts
        msr[selobs&selreg] = md
        common.printlog('compare good z frac to sum of model',logger)
        common.printlog(str(len(ff[selgz&selobs&selreg])/len(ff[selobs&selreg]))+','+str(np.sum(msr[selobs&selreg])/len(ff[selobs&selreg])),logger)

    ff['WEIGHT_ZFAIL'] = wzf
    ff['mod_success_rate'] = msr
    
    #print(len(ff[wz]),len(ff))

    common.printlog('min/max of zfail weights:',logger)
    common.printlog(str(np.min(ff[selobs]['WEIGHT_ZFAIL']))+','+str(np.max(ff[selobs]['WEIGHT_ZFAIL'])),logger)
 
        


    #plt.plot(ff[selgz&selobs]['TSNR2_'+tp[:3]],ff[selgz&selobs]['WEIGHT_ZFAIL'],'k,')
    #plt.xlim(np.percentile(ff[selgz]['TSNR2_'+tp[:3]],0.5),np.percentile(ff[selgz]['TSNR2_'+tp[:3]],99))
    #plt.show()
    
    common.write_LSS_scratchcp(ff,fullname,logger=logger)#,comments='added ZFAIL weight')
    if tp != 'BGS_BRIGHT-21.5':
        ff.keep_columns(['TARGETID','WEIGHT_ZFAIL','mod_success_rate'])
        ffc = Table.read(indir+tp+'_full.dat.fits')
        cols = list(ffc.dtype.names)
        if 'WEIGHT_ZFAIL' in cols:
            ffc.remove_columns(['WEIGHT_ZFAIL'])
        if 'mod_success_rate' in cols:
            ffc.remove_columns(['mod_success_rate'])
        ffc = join(ffc,ff,keys=['TARGETID'],join_type='left')
        common.write_LSS_scratchcp(ffc,indir+tp+'_full.dat.fits',logger=logger)#,comments='added ZFAIL weight')
        fname_mapveto = indir+tp+'_full_HPmapcut.dat.fits'
        if os.path.isfile(fname_mapveto):
            ff.keep_columns(['TARGETID','WEIGHT_ZFAIL','mod_success_rate'])
            ffc = Table.read(fname_mapveto)
            cols = list(ffc.dtype.names)
            if 'WEIGHT_ZFAIL' in cols:
                ffc.remove_columns(['WEIGHT_ZFAIL'])
            if 'mod_success_rate' in cols:
                ffc.remove_columns(['mod_success_rate'])
            ffc = join(ffc,ff,keys=['TARGETID'],join_type='left')
            common.write_LSS_scratchcp(ffc,fname_mapveto,logger=logger)#,comments='added ZFAIL weight')
    
    
#     if dchi2 is not None:
#         if tp[:3] == 'LRG':
#             lrg = ssr_tools.LRG_ssr(surveys=[survey],specrels=[specrel],versions=[version])
#             ffwz = lrg.add_modpre(ff[wz])
#             print(min(ffwz['mod_success_rate']),max(ffwz['mod_success_rate']))
#             ffwz['WEIGHT_ZFAIL'] = 1./ffwz['mod_success_rate']
#             ffwz.keep_columns(['TARGETID','WEIGHT_ZFAIL','mod_success_rate'])
#             ff = join(ff,ffwz,keys=['TARGETID'],join_type='left')
#             #print(min(zf),max(zf))
#             wz = ff['GOODZ']
#             print(len(ff[wz]),len(ff))
#             #ff[wz]['WEIGHT_ZFAIL'] = zf
# 
#         if tp == 'BGS_BRIGHT':
#             bgs = ssr_tools.BGS_ssr(surveys=[survey],specrels=[specrel],versions=[version])
#             ff[wz] = bgs.add_modpre(ff[wz],fl)
#             ff[wz]['WEIGHT_ZFAIL'] = np.clip(1./ff[wz]['mod_success_rate'],1,1.2)
# 
# 
#         if tp == 'ELG_LOP':
#             elg = ssr_tools.ELG_ssr(surveys=[survey],specrels=[specrel],versions=[version])
#             ff[wz] = elg.add_modpre(ff[wz])
# 
#         if tp == 'QSO':
#             qso = ssr_tools.QSO_ssr(surveys=[survey],specrels=[specrel],versions=[version])
#             ff[wz] = qso.add_modpre(ff[wz],fl)
#             print(np.min(ff['WEIGHT_ZFAIL']),np.max(ff['WEIGHT_ZFAIL']))
#             ff['WEIGHT_ZFAIL'] = np.clip(ff['WEIGHT_ZFAIL'],1,2)
#         



def mkclusdat(fl,weighttileloc=True,zmask=False,correct_zcmb='n',tp='',dchi2=9,rcut=None,ntilecut=0,ccut=None,ebits=None,zmin=0,zmax=6,write_cat='y',splitNS='n',return_cat='n',compmd='ran',kemd='',wsyscol=None,use_map_veto='',subfrac=1,zsplit=None, ismock=False,logger=None,extradir='', extracols=None):
    import LSS.common_tools as common
    from LSS import ssr_tools
    '''
    fl is the root of the input/output file
    weighttileloc determines whether to include 1/FRACZ_TILELOCID as a completeness weight
    zmask determines whether to apply a mask at some given redshift
    tp is the target type
    dchi2 is the threshold for keeping as a good redshift
    

    '''
    if write_cat == 'n' and return_cat == 'n':
        print('writing and returning both set to n, this will do nothing so exiting!')
        return True
    
    wzm = '_'
    if ccut is not None:
        wzm = ccut+'_' #you could change this to however you want the file names to turn out

    if zmask:
        wzm += 'zmask_'
    if rcut is not None:
        wzm += 'rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
    if ntilecut > 0:
        wzm += 'ntileg'+str(ntilecut)+'_'
    outf = (fl+wzm+'clustering.dat.fits').replace(tp,extradir+tp)
    ff = Table.read(fl+'_full'+use_map_veto+'.dat.fits'.replace('global','dvs_ro'))
    if wsyscol is not None:
        ff['WEIGHT_SYS'] = np.copy(ff[wsyscol])
    cols = list(ff.dtype.names)
    if 'Z' in cols:
        print('Z column already in full file')
    else:
        ff['Z_not4clus'].name = 'Z'
    if tp[:3] == 'QSO':
        #good redshifts are currently just the ones that should have been defined in the QSO file when merged in full
        wz = ff['Z']*0 == 0
        wz &= ff['Z'] != 999999
        wz &= ff['Z'] != 1.e20
        wz &= ff['ZWARN'] != 999999
        #if not ismock:
        #    wz &= ff['TSNR2_ELG'] > tsnrcut

    if tp[:3] == 'ELG':
        #ff = get_ELG_SSR_tile(ff,dchi2,tsnrcut=tsnrcut)
        wz = ff['ZWARN']*0 == 0
        
        wz &= ff['ZWARN'] != 999999
        if dchi2 is not None:
            wz &= ff['o2c'] > dchi2
            common.printlog('length after oII cut '+str(len(ff[wz])),logger)
        wz &= ff['LOCATION_ASSIGNED'] == 1
        common.printlog('length after also making sure location assigned '+str(len(ff[wz])),logger)
        #if not ismock:
        #    wz &= ff['TSNR2_ELG'] > tsnrcut
        #    common.printlog('length after tsnrcut '+str(len(ff[wz])),logger)

    if tp[:3] == 'LRG':
        print('applying extra cut for LRGs')
        # Custom DELTACHI2 vs z cut from Rongpu
        wz = ff['ZWARN'] == 0
        wz &= ff['ZWARN']*0 == 0
        wz &= ff['ZWARN'] != 999999
        if ismock:
            wz &= ff['Z']<1.5
        else:
            if dchi2 is not None:
                selg = ssr_tools.LRG_goodz(ff)
                wz &= selg

        #wz &= ff['DELTACHI2'] > dchi2
        common.printlog('length after Rongpu cut '+str(len(ff[wz])),logger)
        
        #if not ismock:
        #    wz &= ff['TSNR2_ELG'] > tsnrcut
        #    print('length after tsnrcut '+str(len(ff[wz])))

    if tp[:3] == 'BGS':
        wz = ff['ZWARN'] == 0
        wz &= ff['ZWARN']*0 == 0
        wz &= ff['ZWARN'] != 999999

        if dchi2 is not None:
            print('applying extra cut for BGS')
            wz &= ff['DELTACHI2'] > dchi2
            common.printlog('length after dchi2 cut '+str(len(ff[wz])),logger)
        #wz &= ff['TSNR2_BGS'] > tsnrcut
        #print('length after tsnrcut '+str(len(ff[wz])))

    if correct_zcmb == 'y':
        zcmb = common.get_zcmbdipole(ff['RA'],ff['DEC'])
        newz = (1+ff['Z'])*(1+zcmb)-1
        ff['Z'] = newz
        wzm += 'zcmb_'
        common.printlog('corrected redshifts to cmb frame',logger)


    if subfrac != 1:
        subfracl = np.ones(len(ff))
        sub_array = np.random.random(len(ff))
        if zsplit is not None:
            #subfrac = np.ones(len(ff))
            selzsub = ff['Z'] < zsplit
            
            subfracl[selzsub] = subfrac[0]
            subfracl[~selzsub] = subfrac[1]
        else:
            subfracl *= subfrac
        keep = sub_array < subfracl
        wz &= keep
        
    ff = ff[wz]
    common.printlog('length after cutting to good z '+str(len(ff)),logger)
    ff['WEIGHT'] = np.ones(len(ff))#ff['WEIGHT_ZFAIL']
    if 'WEIGHT_ZFAIL' not in cols:
        ff['WEIGHT_ZFAIL'] = np.ones(len(ff))
    ff['WEIGHT'] *= ff['WEIGHT_ZFAIL']    
#     if dchi2 is not None and 'WEIGHT_ZFAIL' not in cols:
#         if tp[:3] == 'LRG':
#             lrg = ssr_tools.LRG_ssr()
#             ff = lrg.add_modpre(ff)
#             ff['WEIGHT_ZFAIL'] = 1./ff['mod_success_rate']
#             print('min/max of zfail weights:')
#             print(np.min(ff['WEIGHT_ZFAIL']),np.max(ff['WEIGHT_ZFAIL']))
# 
#             print('checking sum of zfail weights compared to length of good z')
#             print(len(ff),np.sum(ff['WEIGHT_ZFAIL']))
#             ff['WEIGHT'] *= ff['WEIGHT_ZFAIL']
# 
#         if tp == 'BGS_BRIGHT':
#             bgs = ssr_tools.BGS_ssr()
#             ff = bgs.add_modpre(ff,fl)
#             ff['WEIGHT_ZFAIL'] = np.clip(1./ff['mod_success_rate'],1,1.2)
#             print('min/max of zfail weights:')
#             print(np.min(ff['WEIGHT_ZFAIL']),np.max(ff['WEIGHT_ZFAIL']))
#             print('checking sum of zfail weights compared to length of good z')
#             print(len(ff),np.sum(ff['WEIGHT_ZFAIL']))
#             ff['WEIGHT'] *= ff['WEIGHT_ZFAIL']
# 
# 
#         if tp == 'ELG_LOP':
#             elg = ssr_tools.ELG_ssr()
#             ff = elg.add_modpre(ff)
#             print('min/max of zfail weights:')
#             print(np.min(ff['WEIGHT_ZFAIL']),np.max(ff['WEIGHT_ZFAIL']))
# 
#             print('checking sum of zfail weights compared to length of good z')
#             print(len(ff),np.sum(ff['WEIGHT_ZFAIL']))
#             ff['WEIGHT'] *= ff['WEIGHT_ZFAIL']
# 
#         if tp == 'QSO':
#             qso = ssr_tools.QSO_ssr()
#             ff = qso.add_modpre(ff,fl)
#             print(np.min(ff['WEIGHT_ZFAIL']),np.max(ff['WEIGHT_ZFAIL']))
#             ff['WEIGHT_ZFAIL'] = np.clip(ff['WEIGHT_ZFAIL'],1,2)
#             print('min/max of zfail weights:')
#             print(np.min(ff['WEIGHT_ZFAIL']),np.max(ff['WEIGHT_ZFAIL']))
#             print('checking sum of zfail weights compared to length of good z')
#             print(len(ff),np.sum(ff['WEIGHT_ZFAIL']))
#             ff['WEIGHT'] *= ff['WEIGHT_ZFAIL']

        
    
    #if tp[:3] == 'ELG':
    #    ff['WEIGHT_ZFAIL'] = 1./ff['relSSR_tile']
    
    if weighttileloc == True:
        ff['WEIGHT_COMP'] = 1./ff['FRACZ_TILELOCID']
        if 'FRAC_TLOBS_TILES' in cols and compmd == 'dat':
            ff['WEIGHT_COMP'] *= 1/ff['FRAC_TLOBS_TILES']

        ff['WEIGHT'] *= ff['WEIGHT_COMP']
    else:
        print('using PROB_OBS for WEIGHT_COMP')
        ff['WEIGHT_COMP'] = 129/(1+128*ff['PROB_OBS'])
        ff['WEIGHT'] *= ff['WEIGHT_COMP']
#    if 'WEIGHT_SYS' not in cols:
#        ff['WEIGHT_SYS'] =  np.ones(len(ff)) #need to initialize these at 1
#    ff['WEIGHT'] *= ff['WEIGHT_SYS']
    if 'WEIGHT_SYS' not in cols:
        ff['WEIGHT_SYS'] =  np.ones(len(ff)) #need to initialize these at 1
    sel = ff['WEIGHT_SYS']*0 != 0
    print(str(len(ff[sel]))+ ' with nan weight_sys being give a value of 1')
    ff['WEIGHT_SYS'][sel] = 1
    print('weightsys bounds',min(ff['WEIGHT_SYS']),max(ff['WEIGHT_SYS']))
    ff['WEIGHT'] *= ff['WEIGHT_SYS']

    #weights for imaging systematic go here
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

#    ff['WEIGHT'] *= ff['WEIGHT_SYS']


    #apply cut on ntile
    if ntilecut > 0:
        print('length before ntile cut '+str(len(ff)))
        wt = ff['NTILE'] > ntilecut
        ff = ff[wt]
        print('length after ntile cut '+str(len(ff)))
    if ccut == 'zQSO':
        wc = ff['SPECTYPE'] ==  'QSO'
        print('length before cutting to spectype QSO '+str(len(ff)))
        ff = ff[wc]
        print('length after cutting to spectype QSO '+str(len(ff)))

    #select down to specific columns below and then also split N/S
    
    common.printlog('cutting to z range '+str(zmin)+' '+str(zmax),logger)
    selz = ff['Z'] > zmin
    selz &= ff['Z'] < zmax
    ff = ff[selz]


    kl = ['RA','DEC','Z','WEIGHT','TARGETID','NTILE','WEIGHT_SYS','WEIGHT_COMP','WEIGHT_ZFAIL','PHOTSYS','FRAC_TLOBS_TILES','TILEID']#,'WEIGHT_FKP']
    if 'WEIGHT_FKP' in cols:
        kl.append('WEIGHT_FKP')
    if 'WEIGHT_SN' in cols:
        kl.append('WEIGHT_SN')
    if 'WEIGHT_RF' in cols:
        kl.append('WEIGHT_RF')
    if 'WEIGHT_IMLIN' in cols:
        kl.append('WEIGHT_IMLIN')
    if 'BITWEIGHTS' in cols:
        kl.append('BITWEIGHTS')
        kl.append('PROB_OBS')
    if 'WEIGHT_NT_MISSPW' in cols:
        kl.append('WEIGHT_NT_MISSPW')
    if extracols is not None:
        if isinstance(extracols, list):
            for ex in extracols:
                kl.append(ex)
        else:
            kl.append(extracols)
    if tp[:3] == 'BGS':
        #ff['flux_r_dered'] = ff['FLUX_R']/ff['MW_TRANSMISSION_R']
        #kl.append('flux_r_dered')
        #print(kl)
        if not ismock:
            fcols = ['G','R','Z','W1','W2']
            ff = common.add_dered_flux(ff,fcols)
            for col in fcols:
                kl.append('flux_'+col.lower()+'_dered')
            print(kl)
            if kemd == 'phot':
                restcols = ['REST_GMR_0P1','REST_GMR_0P0','ABSMAG_RP0','ABSMAG_RP1']
                for col in restcols:
                    kl.append(col)

        if ccut == '-21.5':
            if ismock:
                abr = ff['R_MAG_ABS']
            else:
                from LSS.tabulated_cosmo import TabulatedDESI
                cosmo = TabulatedDESI()
                dis_dc = cosmo.comoving_radial_distance
                dm = 5.*np.log10(dis_dc(ff['Z'])*(1.+ff['Z'])) + 25.
                r_dered = 22.5 - 2.5*np.log10(ff['flux_r_dered'])
                abr = r_dered -dm
            sel = abr < float(ccut)
            print('comparison before/after abs mag cut')
            print(len(ff),len(ff[sel]))
            ff = ff[sel]
        if 'G_R_OBS' in cols:
            kl.append('G_R_OBS')
        if 'G_R_REST' in cols:
            kl.append('G_R_REST')
        if 'R_MAG_ABS' in cols:
            kl.append('R_MAG_ABS')
        if 'R_MAG_APP' in cols:
            kl.append('R_MAG_APP')

    wn = ff['PHOTSYS'] == 'N'
    kll = []
    data_cols = list(ff.dtype.names)
    for name in kl:
        if name in data_cols:
            kll.append(name)
        else:
            print(name+' not found in input and will not be in clustering catalog')
    ff.keep_columns(kll)
    print('minimum,maximum weight')
    print(np.min(ff['WEIGHT']),np.max(ff['WEIGHT']))

    #comments = ["DA02 'clustering' LSS catalog for data, all regions","entries are only for data with good redshifts"]
    #common.write_LSS(ff,outf,comments)

    if write_cat == 'y':
        if splitNS == 'y':
            outfn = fl+wzm+'N_clustering.dat.fits'
            comments = ["DA02 'clustering' LSS catalog for data, BASS/MzLS region","entries are only for data with good redshifts"]
            common.write_LSS(ff[wn],outfn.replace(tp,extradir+tp),comments)

            outfn = fl+wzm+'S_clustering.dat.fits'
            comments = ["DA02 'clustering' LSS catalog for data, DECaLS region","entries are only for data with good redshifts"]
            ffs = ff[~wn]
            common.write_LSS(ffs,outfn.replace(tp,extradir+tp),comments)
        else:
            outfn = fl+wzm+'clustering.dat.fits'
            common.write_LSS_scratchcp(ff,outfn.replace(tp,extradir+tp))
    if return_cat == 'y':
        if splitNS == 'y':
            return ff[wn],ff[~wn]
        else:
            return ff
#     for reg,com in zip(['DS','DN'],[' SGC ',' NGC ']): #split DECaLS NGC/SGC
#         outfn = fl+wzm+reg+'_clustering.dat.fits'
#         sel = densvar.sel_reg(ffs['RA'],ffs['DEC'],reg)
#         comments = ["DA02 'clustering' LSS catalog for data, DECaLS"+com+"region","entries are only for data with good redshifts"]
#         common.write_LSS(ffs[sel],outfn,comments)

def add_tlobs_ran(fl,rann,hpmapcut='',wo=True,logger=None):
    import LSS.common_tools as common
    rf_name = fl+str(rann)+'_full'+hpmapcut+'.ran.fits'
    ranf = Table(fitsio.read(rf_name.replace('global','dvs_ro')))
    tlf = fitsio.read(fl+'frac_tlobs.fits')
    tldic = dict(zip(tlf['TILES'],tlf['FRAC_TLOBS_TILES']))
    #tlarray = np.zeros(len(ranf))
    tlarray = []
    nt = 0
    utls = np.unique(ranf['TILES'])
    gtls = np.isin(utls,tlf['TILES'])
    #for gd,tls in zip(gtls,utls):
    #    sel = ranf['TILES'] == tls
    #    if gd:
    #        fr = tldic[tls]
    #        tlarray[sel] = fr
    #    nt += 1
    #    if nt%1000 == 0:
    #        print(nt,len(utls))
    
    #for i in range(0,len(ranf)):
    for tls in ranf['TILES']:
    #    tls = ranf['TILES'][i]
    #    if tls in tlf['TILES']:
        try:    
            fr = tldic[tls]
        except:
            #print(tls +' not found')
            fr = 1
        tlarray.append(fr)
        if nt%100000 == 0:
           common.printlog(str(nt)+','+str(len(ranf)),logger)  
        nt += 1  
    #        tlarray[i] = fr
    tlarray = np.array(tlarray)
    sel = tlarray == 0
    common.printlog(str(len(tlarray[sel]))+', number with 0 frac; fraction '+str(len(tlarray[sel])/len(tlarray)),logger)
    ranf['FRAC_TLOBS_TILES'] = tlarray
    outf = rf_name#fl+str(rann)+'_full.ran.fits'
    
    common.write_LSS_scratchcp(ranf,outf,logger=logger)
    del ranf
    return True
  
def add_tlobs_ran_array(ranf,tlf,logger=None):
    import LSS.common_tools as common
    tldic = dict(zip(tlf['TILES'],tlf['FRAC_TLOBS_TILES']))
    common.printlog('adding FRAC_TLOBS_TILES',logger)
    tlarray = []
    nt = 0
    utls = np.unique(ranf['TILES'])
    gtls = np.isin(utls,tlf['TILES'])
    nnf = 0
    for tls in ranf['TILES']:
        try:    
            fr = tldic[tls]
        except:
            fr = 1
            nnf += 1
        tlarray.append(fr)
        if nt%100000 == 0:
           common.printlog(str(nt)+','+str(len(ranf)),logger)  
        nt += 1  
    tlarray = np.array(tlarray)
    sel = tlarray == 0
    common.printlog('number of tiles not found in the data '+str(nnf),logger)
    common.printlog(str(len(tlarray[sel]))+' number with 0 frac',logger)
    ranf['FRAC_TLOBS_TILES'] = tlarray
    return ranf
  
    
def mkclusran(flin,fl,rann,rcols=['Z','WEIGHT'],zmask=False,utlid=False,ebits=None,write_cat='y',nosplit='y',return_cat='n',compmd='ran',clus_arrays=None,use_map_veto='',add_tlobs='n',logger=None,extradir='',tp=''):#,tsnrcut=80,tsnrcol='TSNR2_ELG'
    import LSS.common_tools as common
    rng = np.random.default_rng(seed=rann)
    #first find tilelocids where fiber was wanted, but none was assigned; should take care of all priority issues
    wzm = ''
    if zmask:
        wzm += 'zmask_'
    ws = ''
    if utlid:
        ws = 'utlid_'
    if isinstance(flin, str):
        in_fname = flin+str(rann)+'_full'+use_map_veto+'.ran.fits'
        
        ran_cols = ['RA','DEC','TARGETID','TILEID','NTILE','PHOTSYS']#,tsnrcol]
        #if add_tlobs == 'n':
        #    try:
        #        with fitsio.read(in_fname,columns=['FRAC_TLOBS_TILES'],rows=1) as test_f:
        #            ran_cols.append('FRAC_TLOBS_TILES')
        #    except:
        #        common.printlog('failed to find FRAC_TLOBS_TILES, will need to add it',logger)
        #        add_tlobs = 'y'
        #        ran_cols.append('TILES')    
        #else:
        if add_tlobs == 'y':
            ran_cols.append('TILES')
        else:
            ran_cols.append('FRAC_TLOBS_TILES')
        #ffc = Table(fitsio.read(in_fname.replace('global','dvs_ro'),columns=ran_cols))
        ffc = Table(fitsio.read(in_fname,columns=ran_cols))
        common.printlog('loaded '+in_fname,logger)
        #wz = ffr[tsnrcol] > tsnrcut
        #ffc = ffr#[wz]
        #common.printlog(str(rann)+' length after,before tsnr cut:'+' '+str(len(ffc))+','+str(len(ffr)),logger)
        #print(len(ffc),len(ffr))
        #del ffr
        if add_tlobs == 'y':
            
            tlf = fitsio.read(flin+'frac_tlobs.fits')
            ffc = add_tlobs_ran_array(ffc,tlf,logger)
    else:
        ffc = flin
        del flin
        ran_cols = ['RA','DEC','TARGETID','TILEID','NTILE','PHOTSYS','FRAC_TLOBS_TILES']
        ffc.keep_columns(ran_cols)
        
    
    if return_cat == 'y' and nosplit=='y':
        tempcols = ['RA','DEC','TARGETID','NTILE','FRAC_TLOBS_TILES','PHOTSYS']
        if 'WEIGHT_NT_MISSPW' in ffc.columns:
            tempcols.append('WEIGHT_NT_MISSPW')
        #this option will then pass the arrays to the clusran_resamp_arrays function
        ffc.keep_columns(tempcols)
        return ffc
    if utlid:
        ffc = unique(ffc,keys=['TILELOCID'])
        print('length after cutting to unique tilelocid '+str(len(ffc)))
    #def _resamp(selregr,selregd,ffr,fcdn):
    def _resamp(rand_sel,dat_sel,ffr,fcdn):
        for col in rcols:
            ffr[col] =  np.zeros_like(fcdn[col],shape=len(ffr))
        #rand_sel = [selregr,~selregr]
        #dat_sel = [ selregd,~selregd]
        for dsel,rsel in zip(dat_sel,rand_sel):
            #print('aqui len fcdn dsel and ffr rsel', len(fcdn[dsel]),len(ffr[rsel]))
            #inds = np.random.choice(len(fcdn[dsel]),len(ffr[rsel]))
            inds = rng.choice(len(fcdn[dsel]),len(ffr[rsel]))
            dshuf = fcdn[dsel][inds]
            for col in rcols:
                ffr[col][rsel] = dshuf[col]
        if compmd == 'ran':
            ffr['WEIGHT'] *= ffr['FRAC_TLOBS_TILES'] 
        rdl = []
        for dsel,rsel in zip(dat_sel,rand_sel):
            rd = np.sum(ffr[rsel]['WEIGHT'])/np.sum(fcdn[dsel]['WEIGHT'])
            rdl.append(rd)
        for i in range(0,len(rand_sel)-1):
            rdr = rdl[0]/rdl[i+1]
            #print('norm factor is '+str(rdr))
            ffr['WEIGHT'][rand_sel[i+1]] *= rdr
        for dsel,rsel in zip(dat_sel,rand_sel):
            rd = np.sum(ffr[rsel]['WEIGHT'])/np.sum(fcdn[dsel]['WEIGHT'])
            common.printlog(str(rann)+' data/random weighted ratio after resampling:'+str(rd),logger)

#         tabsr = []
#         ffrn = ffr[selregr]
#         ffrs = ffr[~selregr]
#         fcdnn = fcdn[selregd]
#         fcdns = fcdn[~selregd]
#         tabsr = [ffrn,ffrs]
#         tabsd = [fcdnn,fcdns]
#         rdl =[]
#         for i in range(0,len(tabsr)):
#             inds = np.random.choice(len(tabsd[i]),len(tabsr[i]))
#             dshuf = tabsd[i][inds]
#             for col in rcols:
#                 tabsr[i][col] =  dshuf[col]
#             if compmd == 'ran':
#                 tabsr[i]['WEIGHT'] *= tabsr[i]['FRAC_TLOBS_TILES']
#             rd = np.sum(tabsr[i]['WEIGHT'])/np.sum(tabsd[i]['WEIGHT'])
#             rdl.append(rd)
#         rdr = rdl[0]/rdl[1]
#         print('norm factor is '+str(rdr))
#         tabsr[1]['WEIGHT'] *= rdr
#         #print(np.sum(tabsr[0]['WEIGHT'])/np.sum(tabsd[0]['WEIGHT']),np.sum(tabsr[1]['WEIGHT'])/np.sum(tabsd[1]['WEIGHT']))
#         ffr = vstack(tabsr)   
        #print(len(ffr),len_o)     
        return ffr


    if nosplit == 'n':
        regl = ['N_','S_']
    else:
        regl = ['']
    tabl = []
    for ind in range(0,len(regl)):
        reg = regl[ind]
        if clus_arrays is None:
            fcdn = Table.read(fl+wzm+reg+'clustering.dat.fits')#.replace(tp,extradir+tp))
        else:
            fcdn = Table(np.copy(clus_arrays[ind]))
        fcdn.rename_column('TARGETID', 'TARGETID_DATA')
        kc = ['RA','DEC','Z','WEIGHT','TARGETID','NTILE','FRAC_TLOBS_TILES','PHOTSYS']
        if 'WEIGHT_NT_MISSPW' in ffc.columns:
            kc.append('WEIGHT_NT_MISSPW')

        rcols = np.array(rcols)
        wc = np.isin(rcols,list(fcdn.dtype.names))
        rcols = rcols[wc]
        common.printlog(str(rann)+' columns sampled from data are:' +str(rcols),logger)
        #print(rcols)

        if reg != '':
            wn = ffc['PHOTSYS'] == reg.strip('_')
            ffcn = ffc[wn]
        else:
            ffcn = ffc
        outfn =  fl+ws+wzm+reg+str(rann)+'_clustering.ran.fits'#).replace(tp,extradir+tp)  
        
        des_resamp = False
        if 'QSO' in tp:
            if 'S' in reg or reg == '':
                des_resamp = True
        if reg == '' and des_resamp == False: #N/S resampling
            selregr = ffcn['PHOTSYS'] ==  'N'
            selregd = fcdn['PHOTSYS'] ==  'N'
            rand_sel = [selregr,~selregr]
            dat_sel = [selregd,~selregd]
            #print('rand_sel', set(rand_sel[0]), set(rand_sel[1]))
            #print('dat_sel', set(dat_sel[0]), set(dat_sel[1]))
            ffcn = _resamp(rand_sel,dat_sel,ffcn,fcdn)

        if des_resamp and reg == '':
            common.printlog(str(rann) +' resampling in DES region',logger)
            from regressis import footprint
            foot = footprint.DR9Footprint(256, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
            north, south, des = foot.get_imaging_surveys()
            th_ran,phi_ran = (-ffcn['DEC']+90.)*np.pi/180.,ffcn['RA']*np.pi/180.
            th_dat,phi_dat = (-fcdn['DEC']+90.)*np.pi/180.,fcdn['RA']*np.pi/180.
            pixr = hp.ang2pix(256,th_ran,phi_ran,nest=True)
            selregr_des = des[pixr]
            pixd = hp.ang2pix(256,th_dat,phi_dat,nest=True)
            selregd_des = des[pixd]
            selregr = ffcn['PHOTSYS'] ==  'N'
            selregd = fcdn['PHOTSYS'] ==  'N'
            rand_sel = [selregr,selregr_des,~selregr&~selregr_des]
            dat_sel = [ selregd,selregd_des,~selregd&~selregd_des]

            ffcn = _resamp(rand_sel,dat_sel,ffcn,fcdn)
        no_resamp = False
        if reg == 'N' or (reg == 'S' and des_resamp == False):
            no_resamp = True
        if no_resamp:
            common.printlog('Not doing any re-sampling',logger)
            inds = np.random.choice(len(fcdn),len(ffcn))
            dshuf = fcdn[inds]
            for col in rcols:
                ffcn[col] = dshuf[col]
            if compmd == 'ran':
                ffcn['WEIGHT'] *= ffcn['FRAC_TLOBS_TILES']

        for col in rcols:
            kc.append(col)

        ffcn.keep_columns(kc)
    
        if write_cat == 'y':
            #comments seem to cause I/O issues
            #comments = ["'clustering' LSS catalog for random number "+str(rann)+", "+reg+" region","entries are only for data with good redshifts"]
            #common.write_LSS(ffcn,outfn)#,comments)
            common.write_LSS_scratchcp(ffcn,outfn,logger=logger)
        tabl.append(ffcn)
    #outfs =  fl+ws+wzm+'S_'+str(rann)+'_clustering.ran.fits'
    #if clus_arrays is None:
    #    fcds = Table.read(fl+wzm+'S_clustering.dat.fits')
    #else:
    #    fcds = Table(np.copy(clus_arrays[1]))
    #fcds.rename_column('TARGETID', 'TARGETID_DATA')
    #ffcs = ffc[~wn]
    #inds = np.random.choice(len(fcds),len(ffcs))
    #dshuf = fcds[inds]
    #for col in rcols:
    #    ffcs[col] = dshuf[col]
    #ffcs.keep_columns(kc)
    #if write_cat == 'y':
    #    comments = ["'clustering' LSS catalog for random number "+str(rann)+", DECaLS region","entries are only for data with good redshifts"]
    #    common.write_LSS(ffcs,outfs,comments)
    
    if return_cat == 'y':
        return tabl#ffcn,ffcs
#     for reg,com in zip(['DS','DN'],[' SGC ',' NGC ']): #split DECaLS NGC/SGC
#         outfn = fl+wzm+reg+'_'+str(rann)+'_clustering.ran.fits'
#         sel = densvar.sel_reg(ffcs['RA'],ffcs['DEC'],reg)
#         fcd = Table.read(fl+wzm+reg+'_clustering.dat.fits')
#         ffss = ffcs[sel]
#         inds = np.random.choice(len(fcd),len(ffss))
#         dshuf = fcd[inds]
#         for col in rcols:
#             ffss[col] = dshuf[col]
# 
#         comments = ["DA02 'clustering' LSS catalog for random number "+str(rann)+", DECaLS"+com+"region","entries are only for data with good redshifts"]
#         common.write_LSS(ffss,outfn,comments)

def clusran_resamp(flin,rann,rcols=['Z','WEIGHT'],write_cat='y',compmd='ran',logger=None):
    #take existing data/random clustering catalogs and re-sample redshift dependent quantities to assign to randoms
    import LSS.common_tools as common
    rng = np.random.default_rng(seed=rann)
    ffr = Table.read(flin+'_'+str(rann)+'_clustering.ran.fits')
    for col in rcols:
        try:
            ffr.remove_columns([col])
        except:
            print(col+' not in original randoms')
    fcdn = Table.read(flin+'_clustering.dat.fits')
    fcdn.rename_column('TARGETID', 'TARGETID_DATA')
    kc = ['RA','DEC','Z','WEIGHT','TARGETID','NTILE','FRAC_TLOBS_TILES','PHOTSYS']
    rcols = np.array(rcols)
    wc = np.isin(rcols,list(fcdn.dtype.names))
    rcols = rcols[wc]
    print('columns sampled from data are:')
    print(rcols)
    for col in rcols:
        kc.append(col)


    outfn =  flin+'_'+str(rann)+'_clustering.ran.fits'
    
    len_o = len(ffr)
    def _resamp(selregr,selregd,ffr,fcdn):
        for col in rcols:
            ffr[col] =  np.zeros_like(fcdn[col],shape=len(ffr))
        rand_sel = [selregr,~selregr]
        dat_sel = [ selregd,~selregd]
        for dsel,rsel in zip(dat_sel,rand_sel):
            #inds = np.random.choice(len(fcdn[dsel]),len(ffr[rsel]))
            inds = rng.choice(len(fcdn[dsel]),len(ffr[rsel]))
            print(len(fcdn[dsel]),len(inds),np.max(inds))
            dshuf = fcdn[dsel][inds]
            for col in rcols:
                ffr[col][rsel] = dshuf[col]
        ffr['WEIGHT'] *= ffr['FRAC_TLOBS_TILES'] 
        rdl = []
        for dsel,rsel in zip(dat_sel,rand_sel):
            rd = np.sum(ffr[rsel]['WEIGHT'])/np.sum(fcdn[dsel]['WEIGHT'])
            rdl.append(rd)
        rdr = rdl[0]/rdl[1]
        print('norm factor is '+str(rdr))
        ffr['WEIGHT'][rand_sel[1]] *= rdr
        #check that everything worked
        for dsel,rsel in zip(dat_sel,rand_sel):
            rd = np.sum(ffr[rsel]['WEIGHT'])/np.sum(fcdn[dsel]['WEIGHT'])
            print('data/random weighted ratio after resampling:'+str(rd))
        

#         tabsr = []
#         ffrn = ffr[selregr]
#         ffrs = ffr[~selregr]
#         fcdnn = fcdn[selregd]
#         fcdns = fcdn[~selregd]
#         tabsr = [ffrn,ffrs]
#         tabsd = [fcdnn,fcdns]
#         rdl =[]
#         for i in range(0,len(tabsr)):
#             inds = np.random.choice(len(tabsd[i]),len(tabsr[i]))
#             dshuf = tabsd[i][inds]
#             for col in rcols:
#                 tabsr[i][col] =  dshuf[col]
#             if compmd == 'ran':
#                 tabsr[i]['WEIGHT'] *= tabsr[i]['FRAC_TLOBS_TILES']
#             rd = np.sum(tabsr[i]['WEIGHT'])/np.sum(tabsd[i]['WEIGHT'])
#             rdl.append(rd)
#         rdr = rdl[0]/rdl[1]
#         print('norm factor is '+str(rdr))
#         tabsr[1]['WEIGHT'] *= rdr
#         #print(np.sum(tabsr[0]['WEIGHT'])/np.sum(tabsd[0]['WEIGHT']),np.sum(tabsr[1]['WEIGHT'])/np.sum(tabsd[1]['WEIGHT']))
#         ffr = vstack(tabsr)   
        #print(len(ffr),len_o)     
        return ffr

    if 'NGC' in flin:
        #need to split N/S when sampling
        print('doing N/S re-sampling')
        selregr = ffr['DEC'] > 32.375
        selregd = fcdn['DEC'] > 32.375
        ffr = _resamp(selregr,selregd,ffr,fcdn)
    elif 'SGC' in flin and 'QSO' in flin:
        print('resampling in DES region')
        from regressis import footprint
        foot = footprint.DR9Footprint(256, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
        north, south, des = foot.get_imaging_surveys()
        th_ran,phi_ran = (-ffr['DEC']+90.)*np.pi/180.,ffr['RA']*np.pi/180.
        th_dat,phi_dat = (-fcdn['DEC']+90.)*np.pi/180.,fcdn['RA']*np.pi/180.
        pixr = hp.ang2pix(256,th_ran,phi_ran,nest=True)
        selregr = des[pixr]
        pixd = hp.ang2pix(256,th_dat,phi_dat,nest=True)
        selregd = des[pixd]
        ffr = _resamp(selregr,selregd,ffr,fcdn)
    else:
        inds = np.random.choice(len(fcdn),len(ffr))
        dshuf = fcdn[inds]
        for col in rcols:
            ffr[col] = dshuf[col]
        if compmd == 'ran':
            ffr['WEIGHT'] *= ffr['FRAC_TLOBS_TILES']
        #kc.append(col)
    ffr.keep_columns(kc)
    
    if write_cat == 'y':
        #comments = ["'clustering' LSS catalog for random number "+str(rann)+", BASS/MzLS region","entries are only for data with good redshifts"]
        common.write_LSS_scratchcp(ffr,outfn,logger=logger)

def clusran_resamp_arrays(ffr,fcdn,reg,tracer,rcols=['Z','WEIGHT'],compmd='ran'):
    #take existing data/random clustering catalogs and re-sample redshift dependent quantities to assign to randoms
    import LSS.common_tools as common
    
    ffr = Table(ffr)
    for col in rcols:
        try:
            ffr.remove_columns([col])
        except:
            print(col+' not in original randoms')
    fcdn.rename_column('TARGETID', 'TARGETID_DATA')
    kc = ['RA','DEC','Z','WEIGHT','TARGETID','NTILE','FRAC_TLOBS_TILES']
    for col in rcols:
        kc.append(col)
    rcols = np.array(rcols)
    wc = np.isin(rcols,list(fcdn.dtype.names))
    rcols = rcols[wc]
    print('columns sampled from data are:')
    print(rcols)


    
    len_o = len(ffr)
    def _resamp(selregr,selregd,ffr,fcdn):
        tabsr = []
        ffrn = ffr[selregr]
        ffrs = ffr[~selregr]
        fcdnn = fcdn[selregd]
        fcdns = fcdn[~selregd]
        tabsr = [ffrn,ffrs]
        tabsd = [fcdnn,fcdns]
        rdl =[]
        for i in range(0,len(tabsr)):
            inds = np.random.choice(len(tabsd[i]),len(tabsr[i]))
            dshuf = tabsd[i][inds]
            #print(dshuf.dtype.names)
            for col in rcols:
                tabsr[i][col] =  dshuf[col]
            if compmd == 'ran':
                tabsr[i]['WEIGHT'] *= tabsr[i]['FRAC_TLOBS_TILES']
            rd = np.sum(tabsr[i]['WEIGHT'])/np.sum(tabsd[i]['WEIGHT'])
            rdl.append(rd)
        rdr = rdl[0]/rdl[1]
        print('norm factor is '+str(rdr))
        tabsr[1]['WEIGHT'] *= rdr
        #print(np.sum(tabsr[0]['WEIGHT'])/np.sum(tabsd[0]['WEIGHT']),np.sum(tabsr[1]['WEIGHT'])/np.sum(tabsd[1]['WEIGHT']))
        ffr = vstack(tabsr)   
        #print(len(ffr),len_o)     
        return ffr

    if reg == 'NGC':
        #need to split N/S when sampling
        selregr = ffr['DEC'] > 32.375
        selregd = fcdn['DEC'] > 32.375
        ffr = _resamp(selregr,selregd,ffr,fcdn)
        
    elif reg == 'SGC' and tracer == 'QSO':
        print('resampling in DES region')
        from regressis import footprint
        foot = footprint.DR9Footprint(256, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
        north, south, des = foot.get_imaging_surveys()
        th_ran,phi_ran = (-ffr['DEC']+90.)*np.pi/180.,ffr['RA']*np.pi/180.
        th_dat,phi_dat = (-fcdn['DEC']+90.)*np.pi/180.,fcdn['RA']*np.pi/180.
        pixr = hp.ang2pix(256,th_ran,phi_ran,nest=True)
        selregr = des[pixr]
        pixd = hp.ang2pix(256,th_dat,phi_dat,nest=True)
        selregd = des[pixd]
        ffr = _resamp(selregr,selregd,ffr,fcdn)
        
    else:
        inds = np.random.choice(len(fcdn),len(ffr))
        dshuf = fcdn[inds]
        for col in rcols:
            ffr[col] = dshuf[col]
        if compmd == 'ran':
            ffr['WEIGHT'] *= ffr['FRAC_TLOBS_TILES']
        #kc.append(col)
    ffr.keep_columns(kc)
    return ffr




def clusNStoGC(flroot,nran=1):
    import LSS.common_tools as common
    '''
    combine N and S then split NGC/SGC
    randoms get re-normalized so that Ndata/Nrandom is the same N and S; they will no longer be able to be used to get areas
    flroot is the common root of the file
    nran is the number of random files to process
    '''
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    fn = Table(fitsio.read(flroot+'N_clustering.dat.fits'))
    nn = np.sum(fn['WEIGHT'])
    fs = Table(fitsio.read(flroot+'S_clustering.dat.fits'))
    ns = np.sum(fs['WEIGHT'])
    if 'QSO' in flroot:
        print('getting data numbers in DES region')
        from regressis import footprint
        foot = footprint.DR9Footprint(256, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
        north, south, des = foot.get_imaging_surveys()
        th_dat,phi_dat = (-fs['DEC']+90.)*np.pi/180.,fs['RA']*np.pi/180.
        pixd = hp.ang2pix(256,th_dat,phi_dat,nest=True)
        selregd = des[pixd]
        ns_des =  np.sum(fs[selregd]['WEIGHT'])
        ns_notdes =  np.sum(fs[~selregd]['WEIGHT'])

    fc = vstack((fn,fs))
    print(np.sum(fc['WEIGHT']),nn,ns)
    c = SkyCoord(fc['RA']* u.deg,fc['DEC']* u.deg,frame='icrs')
    gc = c.transform_to('galactic')
    sel_ngc = gc.b > 0
    outf_ngc = flroot+'NGC_clustering.dat.fits'
    common.write_LSS(fc[sel_ngc],outf_ngc)
    outf_sgc = flroot+'SGC_clustering.dat.fits'
    common.write_LSS(fc[~sel_ngc],outf_sgc)
    
    for rann in range(0,nran):
        fn = Table(fitsio.read(flroot+'N_'+str(rann)+'_clustering.ran.fits'))
        nnr = np.sum(fn['WEIGHT'])
        fs = Table(fitsio.read(flroot+'S_'+str(rann)+'_clustering.ran.fits'))
        if 'QSO' in flroot:
            print('resampling in DES region')
            th_ran,phi_ran = (-fs['DEC']+90.)*np.pi/180.,fs['RA']*np.pi/180.
            pixr = hp.ang2pix(256,th_ran,phi_ran,nest=True)
            selregr = des[pixr]
            nsr_des = np.sum(fs[selregr]['WEIGHT'])
            nsr_notdes = np.sum(fs[~selregr]['WEIGHT'])
        

        nsr = np.sum(fs['WEIGHT'])
        rn = nn/nnr
        rs = ns/nsr
        #we want rn and rs to be the same, this can be effectively accomplished by 
        fac = rs/rn
        fs['WEIGHT'] *= fac
        #double check
        nsr = np.sum(fs['WEIGHT'])
        rs = ns/nsr
        print('checking that random ratios are now the same size',rn,rs)            
        fc = vstack((fn,fs))
        print(np.sum(fc['WEIGHT']),nnr,nsr)

        c = SkyCoord(fc['RA']* u.deg,fc['DEC']* u.deg,frame='icrs')
        gc = c.transform_to('galactic')
        sel_ngc = gc.b > 0
        outf_ngc = flroot+'NGC_'+str(rann)+'_clustering.ran.fits'
        common.write_LSS(fc[sel_ngc],outf_ngc)
        outf_sgc = flroot+'SGC_'+str(rann)+'_clustering.ran.fits'
        common.write_LSS(fc[~sel_ngc],outf_sgc)
   
def splitclusGC(flroot,nran=1,par='n'):
    import LSS.common_tools as common
    '''
    split full clustering catalog by Galactic cap; should already have been re-sampled N/S (and DES for QSO)
    '''
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    fc = Table(fitsio.read(flroot+'clustering.dat.fits'))
    c = SkyCoord(fc['RA']* u.deg,fc['DEC']* u.deg,frame='icrs')
    gc = c.transform_to('galactic')
    sel_ngc = gc.b > 0
    outf_ngc = flroot+'NGC_clustering.dat.fits'
    common.write_LSS(fc[sel_ngc],outf_ngc)
    outf_sgc = flroot+'SGC_clustering.dat.fits'
    common.write_LSS(fc[~sel_ngc],outf_sgc)
    def _ranparfun(rann):
    
        fc = Table(fitsio.read(flroot+str(rann)+'_clustering.ran.fits'))
        c = SkyCoord(fc['RA']* u.deg,fc['DEC']* u.deg,frame='icrs')
        gc = c.transform_to('galactic')
        sel_ngc = gc.b > 0
        outf_ngc = flroot+'NGC_'+str(rann)+'_clustering.ran.fits'
        common.write_LSS(fc[sel_ngc],outf_ngc)
        outf_sgc = flroot+'SGC_'+str(rann)+'_clustering.ran.fits'
        common.write_LSS(fc[~sel_ngc],outf_sgc)
    inds = np.arange(nran)
    if par == 'n':
        for rann in inds:
            _ranparfun(rann)
    if par == 'y':
        from multiprocessing import Pool
        with Pool() as pool:
            res = pool.map(_ranparfun, inds)
    


def random_mtl(rd,outf ):
    '''
    rd is the table containing the randoms
    outf is the file name to save to
    '''
    rmtl = Table(rd)
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


def randomtiles_allmain(tiles,dirout='/global/cfs/cdirs/desi/survey/catalogs/main/LSS/random',imin=0,imax=18,rann=1,dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/' ):
    import desimodel.focalplane
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

def randomtiles_allmain_pix_2step(tiles,dirout='/global/cfs/cdirs/desi/survey/catalogs/main/LSS/random',ii=0,dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/',logger=None ):
    '''
    tiles should be a table containing the relevant info
    '''
    from desitarget.io import read_targets_in_tiles
    import desimodel.focalplane
    import desimodel.footprint
    import LSS.common_tools as common
    common.printlog('making random target files for tiles',logger)
    trad = desimodel.focalplane.get_tile_radius_deg()*1.1 #make 10% greater just in case
    #print(trad)

    nd = 0
    sel_tile = np.zeros(len(tiles),dtype=bool)
    for i in range(0,len(tiles)):

        #print('length of tile file is (expected to be 1):'+str(len(tiles)))
        #tile = tiles[tiles['TILEID']==tiles['TILEID'][i]]
        fname = dirout+str(ii)+'/tilenofa-'+str(tiles['TILEID'][i])+'.fits'
        if os.path.isfile(fname):
            #print(fname +' already exists')
            pass
        else:
            sel_tile[i] = True
    tiles = tiles[sel_tile]
    if len(tiles) == 0:
        common.printlog('no tiles to process for '+str(ii),logger)
        return True
    rtall = read_targets_in_tiles(dirrt,tiles)
    common.printlog('read targets on all tiles',logger)

    common.printlog('creating files for '+str(len(tiles))+' tiles',logger)
    #for i in range(0,len(tiles)):
    def _create_rantile(ind):
        fname = dirout+str(ii)+'/tilenofa-'+str(tiles['TILEID'][ind])+'.fits'
        #print('creating '+fname)
        tdec = tiles['DEC'][ind]
        decmin = tdec - trad
        decmax = tdec + trad
        wdec = (rtall['DEC'] > decmin) & (rtall['DEC'] < decmax)
        #print(len(rt[wdec]))
        inds = desimodel.footprint.find_points_radec(tiles['RA'][ind], tdec,rtall[wdec]['RA'], rtall[wdec]['DEC'])
        #print('got indexes')
        rtw = rtall[wdec][inds]
        rmtl = Table(rtw)
        #print('made table for '+fname)
        del rtw
        #rmtl['TARGETID'] = np.arange(len(rmtl))
        #print(len(rmtl['TARGETID'])) #checking this column is there
        rmtl['DESI_TARGET'] = np.ones(len(rmtl),dtype=int)*2
        rmtl['NUMOBS_INIT'] = np.zeros(len(rmtl),dtype=int)
        rmtl['NUMOBS_MORE'] = np.ones(len(rmtl),dtype=int)
        rmtl['PRIORITY'] = np.ones(len(rmtl),dtype=int)*3400
        rmtl['OBSCONDITIONS'] = np.ones(len(rmtl),dtype=int)*516#tiles['OBSCONDITIONS'][i]
        rmtl['SUBPRIORITY'] = np.random.random(len(rmtl))
        #print('added columns for '+fname)
        rmtl.write(fname,format='fits', overwrite=True)
        del rmtl
        common.printlog('added columns, wrote to '+fname,logger)
        #nd += 1
        #print(str(nd),len(tiles))
    inds = np.arange(len(tiles))
    for ind in inds:
        _create_rantile(ind)
        common.printlog('done with '+str(ii)+' tile '+str(ind)+' out of '+str(len(inds)),logger)
    #from multiprocessing import Pool
    #with Pool() as pool:
    #    res = pool.map(_create_rantile, inds)


def randomtiles_allmain_pix(tiles,dirout='/global/cfs/cdirs/desi/survey/catalogs/main/LSS/random',imin=0,imax=18,dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/' ):
    '''
    tiles should be a table containing the relevant info
    '''
    from desitarget.io import read_targets_in_tiles
    for ii in range(imin,imax):
        nd = 0
        for i in range(0,len(tiles)):

            #print('length of tile file is (expected to be 1):'+str(len(tiles)))
            tile = tiles[tiles['TILEID']==tiles['TILEID'][i]]
            fname = dirout+str(ii)+'/tilenofa-'+str(tiles['TILEID'][i])+'.fits'
            if os.path.isfile(fname):
                #print(fname +' already exists')
                pass
            else:
                print('creating '+fname)
                rtw = read_targets_in_tiles(dirrt,tile)
                print('read targets for '+fname)
                rmtl = Table(rtw)
                del rtw
                #rmtl['TARGETID'] = np.arange(len(rmtl))
                #print(len(rmtl['TARGETID'])) #checking this column is there
                rmtl['DESI_TARGET'] = np.ones(len(rmtl),dtype=int)*2
                rmtl['NUMOBS_INIT'] = np.zeros(len(rmtl),dtype=int)
                rmtl['NUMOBS_MORE'] = np.ones(len(rmtl),dtype=int)
                rmtl['PRIORITY'] = np.ones(len(rmtl),dtype=int)*3400
                rmtl['OBSCONDITIONS'] = np.ones(len(rmtl),dtype=int)*516#tiles['OBSCONDITIONS'][i]
                rmtl['SUBPRIORITY'] = np.random.random(len(rmtl))
                print('added columns for '+fname)
                rmtl.write(fname,format='fits', overwrite=True)
                del rmtl
                print('added columns, wrote to '+fname)
                nd += 1
                print(str(nd))


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
