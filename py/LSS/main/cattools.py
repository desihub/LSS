'''
python functions to do various useful date processing/manipulation
'''
import numpy as np
import glob
import os
from random import random

import astropy.io.fits as fits
from astropy.table import Table,join,unique,vstack

import fitsio

from matplotlib import pyplot as plt

import desimodel.footprint as foot
import desimodel.focalplane

from desitarget.io import read_targets_in_tiles
from desitarget.targetmask import obsmask, obsconditions, zwarn_mask

from desispec.io.emlinefit import read_emlines_inputs
from desispec.emlinefit import get_emlines

import healpy as hp

#from LSS.Cosmo import distance
from LSS.imaging import densvar
from LSS.common_tools import find_znotposs
import LSS.common_tools as common
from LSS import ssr_tools

import logging
logging.getLogger("QSO_CAT_UTILS").setLevel(logging.ERROR)



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



def combtile_spec(tiles,outf='',md='',specver='daily',redo='n',specrel='guadalupe',prog='dark'):
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

            #if s == 0:
            #    specd = tspec
            #    s = 1
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
        tspec = combEMdata_guad(tile,tdate,coaddir=coaddir)
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



def combEMdata_guad(tile,tdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/guadalupe/tiles/cumulative/'):
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
        dt['TILEID'] = tile
        dt.remove_columns(remcol)
        return dt
    else:
        return None

def combEMdata_daily(tile,zdate,tdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/archive/',outf='temp.fits'):
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

# def cutphotmask(aa,bits):
#     print(str(len(aa)) +' before imaging veto' )
#     keep = (aa['NOBS_G']>0) & (aa['NOBS_R']>0) & (aa['NOBS_Z']>0)
#     for biti in bits:
#         keep &= ((aa['MASKBITS'] & 2**biti)==0)
#     aa = aa[keep]
#     print(str(len(aa)) +' after imaging veto' )
#     return aa

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

    tarsn = None
    tls = foot.pix2tiles(8,[hpx],tiles)
    if os.path.isfile(fout):
        tarsn = Table.read(fout)
        s = 1
        tdone = np.unique(tarsn['TILEID'])
        tmask = ~np.isin(tls['TILEID'],tdone)
    else:
        tmask = np.ones(len(tls)).astype('bool')
    print('there are potentially '+str(len(tls[tmask]))+' to get updates from, out of a possible '+str(len(tls))+' overlapping this pixel')
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

            print(tile,n,len(tls[tmask]),len(tarsn))

        else:
            print('no overlapping targetid')
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
    #indir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+specrel
    if ver == 'everest' or ver == 'guadalupe':
        zf = indir+'/datcomb_'+pd+'_tarspecwdup_zdone.fits'
    if ver == 'daily':
        zf = indir+'/datcomb_'+pd+'_spec_zdone.fits'
    print(zf)
    dz = Table.read(zf)
    #dz = fitsio.read(zf)
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
    if badfib is not None:
        bad = np.isin(fs['FIBER'],badfib)
        print('number at bad fibers '+str(sum(bad)))
        wfqa &= ~bad

    return fs[wfqa]

def cut_specdat(dz):
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


def count_tiles_better(dr,pd,rann=0,specrel='daily',fibcol='COADD_FIBERSTATUS',px=False,survey='main'):
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

    indir = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/'+specrel
    ps = pd
    if pd[:3] == 'LRG' or pd[:3] == 'ELG' or pd[:3] =='QSO':
        ps = 'dark'
    if pd[:3] == 'BGS' or pd[:3] == 'MWS_ANY':
        ps = 'bright'
    fs = get_specdat(indir,ps)

    stlid = 10000*fs['TILEID'] +fs['LOCATION']
    gtl = np.unique(stlid)

    if dr == 'dat':
        fj = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/'+specrel+'/datcomb_'+pd+'_tarspecwdup_zdone.fits')
        #outf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/datcomb_'+pd+'ntileinfo.fits'
    if dr == 'ran':
        if px:
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

def count_tiles_better_px(dr,pd,gtl,rann=0,specrel='daily',fibcol='COADD_FIBERSTATUS',px=None):
    '''
    from files with duplicates that have already been sorted by targetid, quickly go
    through and get the multi-tile information
    dr is either 'dat' or 'ran'
    returns file with TARGETID,NTILE,TILES,TILELOCIDS
    '''

    if dr == 'dat':
        fj = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+specrel+'/datcomb_'+pd+'_tarspecwdup_zdone.fits')
        #outf = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/datcomb_'+pd+'ntileinfo.fits'
    if dr == 'ran':
        if px is not None:
            fj = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+specrel+'/healpix/rancomb_'+str(rann)+pd+'_'+str(px)+'_wdupspec_zdone.fits')
        else:
            fj = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+specrel+'/rancomb_'+str(rann)+pd+'wdupspec_zdone.fits')

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
        rv = True
    else:
        rv = False
    #specf = Table.read(lspecdir+'datcomb_'+tp+'_spec_zdone.fits')
    specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    specf.keep_columns(keepcols)
    #specf.keep_columns(['ZWARN','LOCATION','TILEID','TILELOCID','FIBERSTATUS','FIBERASSIGN_X','FIBERASSIGN_Y','PRIORITY','DELTA_X','DELTA_Y','EXPTIME','PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B','TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z','TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
    fgu = join(fgu,specf,keys=['LOCATION','TILEID','FIBER'],join_type='left')
    fgu.sort('TARGETID')
    outf = lspecdir+'/rancomb_'+str(rann)+tp+'wdupspec_zdone.fits'
    print(outf)
    fgu.write(outf,format='fits', overwrite=True)
    return rv

def combran_wdup_hp(hpx,tiles,rann,randir,tp,lspecdir,specf,keepcols=[],outf='',redos=False):

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
        fgu = Table.read(outf)
        s = 1
        tdone = np.unique(fgu['TILEID'])
        tmask = ~np.isin(tls['TILEID'],tdone)
    else:
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

def mkfullran(gtl,lznp,indir,rann,imbits,outf,tp,pd,notqso='',maxp=3400,min_tsnr2=0,tlid_full=None,badfib=None):

    if pd == 'bright':
        tscol = 'TSNR2_BGS'
    else:
        tscol = 'TSNR2_ELG'

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
    #zf = indir+'/datcomb_'+pd+'_tarspecwdup_zdone.fits'
    #dz = Table.read(zf)

    #fs = get_specdat(indir,pd,verspec)
#     stlid = 10000*fs['TILEID'] +fs['LOCATION']
#     gtl = np.unique(stlid)
#
#     wtype = ((dz[desitarg] & bit) > 0)
#     if notqso == 'notqso':
#         wtype &= ((dz[desitarg] & qsobit) == 0)
#
#     wg = np.isin(dz['TILELOCID'],gtl)
#     dz = dz[wtype&wg]
#     print('length after selecting type and good hardware '+str(len(dz)))
#     lznp = find_znotposs(dz)


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
    print('length including duplicates '+str(len(dz)))
    wk = ~np.isin(dz['TILELOCID'],lznp)
    dz['ZPOSSLOC'] = np.zeros(len(dz)).astype('bool')
    dz['ZPOSSLOC'][wk] = 1

    wg = np.isin(dz['TILELOCID'],gtl)
    if badfib is not None:
        bad = np.isin(dz['FIBER'],badfib)
        print('number at bad fibers '+str(sum(bad)))
        wg &= ~bad


    dz['GOODHARDLOC'] = np.zeros(len(dz)).astype('bool')
    dz['GOODHARDLOC'][wg] = 1

    dz['LOCFULL'] = np.zeros(len(dz)).astype('bool')
    if tlid_full is not None:
        wf = np.isin(dz['TILELOCID'],tlid_full)
        dz['LOCFULL'][wf] = 1

    #dz = dz[wk]
    #print('length after cutting to good positions '+str(len(dz)))
    print('joining with original randoms to get mask properties')
    dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'
    tcol = ['TARGETID','MASKBITS','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z'] #only including what are necessary for mask cuts for now
    #tcol = ['TARGETID','EBV','WISEMASK_W1','WISEMASK_W2','BRICKID','PSFDEPTH_G','PSFDEPTH_R','PSFDEPTH_Z','GALDEPTH_G',\
    #'GALDEPTH_R','GALDEPTH_Z','PSFDEPTH_W1','PSFDEPTH_W2','PSFSIZE_G','PSFSIZE_R','PSFSIZE_Z','MASKBITS','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z']
    tarf = fitsio.read(dirrt+'/randoms-1-'+str(rann)+'.fits',columns=tcol)
    dz = join(dz,tarf,keys=['TARGETID'])
    del tarf
    dz = common.cutphotmask(dz,imbits)
    print('length after cutting to based on imaging veto mask '+str(len(dz)))
#     pl = np.copy(dz['PRIORITY']).astype(float)#dz['PRIORITY']
#     sp = pl <= 0
#     pl[sp] = .1
# 
#     dz['sort'] = dz[tsnr]*dz['GOODHARDLOC']*dz['ZPOSSLOC']+dz['GOODHARDLOC']*dz['ZPOSSLOC']+dz['GOODHARDLOC']*dz['ZPOSSLOC']/pl

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
    dz['sort'] =  dz['GOODPRI']*dz['GOODHARDLOC']*dz['ZPOSSLOC']*dz['GOODTSNR']*1-0.5*dz['LOCFULL']#*(1+dz[tsnr])

    #dz['sort'] =  dz['GOODPRI']*dz['GOODHARDLOC']*dz['ZPOSSLOC']#*(1+dz[tsnr])


    dz.sort('sort') #should allow to later cut on tsnr for match to data
    dz = unique(dz,keys=['TARGETID'],keep='last')
    print('length after cutting to unique TARGETID '+str(len(dz)))
    print(np.unique(dz['NTILE']))

    dz.write(outf,format='fits', overwrite=True)
    print('wrote to '+outf)
    del dz

def mkfullran_px(indir,rann,imbits,outf,tp,pd,gtl,lznp,px,dirrt,maxp=3400,min_tsnr2=0,tlid_full=None):
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
            tarf = fitsio.read(dirrt+'/randoms-1-hp-'+str(px)+'.fits',columns=tcol)
            dz = join(dz,tarf,keys=['TARGETID'])
            del tarf

            dz = common.cutphotmask(dz,imbits)
            #print('length after cutting to based on imaging veto mask '+str(len(dz)))
            if len(dz) > 0:
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
                dz['sort'] =  dz['GOODPRI']*dz['GOODHARDLOC']*dz['ZPOSSLOC']*dz['GOODTSNR']*1-0.5*dz['LOCFULL']#*(1+dz[tsnr])

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



def mkfulldat(zf,imbits,ftar,tp,bit,outf,ftiles,azf='',azfm='cumul',desitarg='DESI_TARGET',specver='daily',notqso='',qsobit=4,min_tsnr2=0,badfib=None):
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


    #indir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+specver
    #fs = get_specdat(indir,pd)
    #stlid = 10000*fs['TILEID'] +fs['LOCATION']
    #gtl = np.unique(stlid)

    dz = Table(fitsio.read(zf))
    wtype = ((dz[desitarg] & bit) > 0)
    if notqso == 'notqso':
        print('removing QSO targets')
        wtype &= ((dz[desitarg] & qsobit) == 0)

    #wg = np.isin(dz['TILELOCID'],gtl)
    print(len(dz[wtype]))
    #dz = dz[wtype&wg]
    dz = dz[wtype]

    #instead of full spec data, we are going to get type specific data and cut to unique entries
    #in the end, we can only use the data associated with an observation
    #NOTE, this is not what we want to do for randoms, where instead we want to keep all of the
    #locations where it was possible a target could have been assigned

    fs = common.cut_specdat(dz,badfib)
    #fs['sort'] = fs['TSNR2_LRG']
    #fs.sort('sort')
    #fsu = unique(fs,keys=['TARGETID'],keep='last')
    #gtl = np.unique(fsu['TILELOCID'])
    gtl = np.unique(fs['TILELOCID'])

    wg = np.isin(dz['TILELOCID'],gtl)
    print(len(dz[wg]))
    dz['GOODHARDLOC'] = np.zeros(len(dz)).astype('bool')
    dz['GOODHARDLOC'][wg] = 1
    print('length after selecting type '+str(len(dz)))
    #print('length after selecting to locations where target type was observed '+str(len(dz)))
    #These steps are not needed if we cut already to only locations where the target type was observed
    #lznp = find_znotposs(dz)
    #wk = ~np.isin(dz['TILELOCID'],lznp)#dz['ZPOSS'] == 1
    #dz = dz[wk]
    #print('length after priority veto '+str(len(dz)))
    dtl = Table.read(ftiles)
    dtl.keep_columns(['TARGETID','NTILE','TILES','TILELOCIDS'])
    dz = join(dz,dtl,keys='TARGETID')
    wz = dz['ZWARN'] != 999999 #this is what the null column becomes
    wz &= dz['ZWARN']*0 == 0 #just in case of nans
    dz['LOCATION_ASSIGNED'] = np.zeros(len(dz)).astype('bool')
    dz['LOCATION_ASSIGNED'][wz] = 1
    tlids = np.unique(dz['TILELOCID'][wz])
    wtl = np.isin(dz['TILELOCID'],tlids)
    dz['TILELOCID_ASSIGNED'] = np.zeros(len(dz)).astype('bool')
    dz['TILELOCID_ASSIGNED'][wtl] = 1
    print('number of unique targets at assigned tilelocid:')
    print(len(np.unique(dz[wtl]['TARGETID'])))

    wnts = dz[tscol]*0 != 0
    wnts |= dz[tscol] == 999999
    dz[tscol][wnts] = 0
    print(np.max(dz[tscol]))
    dz['GOODTSNR'] = np.zeros(len(dz)).astype('bool')
    sel = dz[tscol] > min_tsnr2
    dz['GOODTSNR'][sel] = 1

    dz['sort'] = dz['LOCATION_ASSIGNED']*dz['GOODTSNR']*dz['GOODHARDLOC']*(1+np.clip(dz[tscol],0,200))*1+dz['TILELOCID_ASSIGNED']*dz['GOODHARDLOC']*1+dz['GOODHARDLOC']*1

    dz.sort('sort')
    dz = unique(dz,keys=['TARGETID'],keep='last')

    print('length after cutting to unique targets '+str(len(dz)))
    print('joining to full imaging')
    dz.remove_columns(['RA','DEC','DESI_TARGET','BGS_TARGET']) #these come back in with merge to full target file
    dz = join(dz,ftar,keys=['TARGETID'])
    #print('length after join to full targets (should be same) '+str(len(dz)))
    dz = common.cutphotmask(dz,imbits)
    print('length after imaging mask; should not have changed '+str(len(dz)))


    if tp[:3] == 'ELG' and azf != '' and azfm == 'cumul':# or tp == 'ELG_HIP':
        #arz = fitsio.read(azf,columns=['TARGETID','LOCATION','TILEID','OII_FLUX','OII_FLUX_IVAR','SUBSET','DELTACHI2'])
        arz = Table(fitsio.read(azf,columns=['TARGETID','LOCATION','TILEID','OII_FLUX','OII_FLUX_IVAR']))
        dz = join(dz,arz,keys=['TARGETID','LOCATION','TILEID'],join_type='left')#,uniq_col_name='{col_name}{table_name}',table_names=['', '_OII'])
        #st = []
        #for i in range(0,len(arz)):
        #    st.append(arz['SUBSET'][i][:4])
        #st = np.array(st)
        #wg = arz[fbcol] == 0
        #wg = st == "thru"
        #arz = arz[wg]
        #o2c = np.log10(arz['OII_FLUX'] * np.sqrt(arz['OII_FLUX_IVAR']))+0.2*np.log10(arz['DELTACHI2'])
        o2c = np.log10(dz['OII_FLUX'] * np.sqrt(dz['OII_FLUX_IVAR']))+0.2*np.log10(dz['DELTACHI2'])
        w = (o2c*0) != 0
        w |= dz['OII_FLUX'] < 0
        o2c[w] = -20
        #arz.keep_columns(['TARGETID','LOCATION','TILEID','o2c','OII_FLUX','OII_SIGMA'])#,'Z','ZWARN','TSNR2_ELG'])
        #arz = Table(arz)
        #arz['o2c'] = o2c
        dz['o2c'] = o2c
        #dz = join(dz,arz,keys=['TARGETID','LOCATION','TILEID'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['', '_OII'])

        #dz.remove_columns(['SUBSET','DELTACHI2_OII'])#,fbcol+'_OII'])
        print('check length after merge with OII strength file:' +str(len(dz)))

    if tp[:3] == 'QSO' and azf != '' and azfm == 'cumul':
        arz = Table(fitsio.read(azf))
        arz.keep_columns(['TARGETID','LOCATION','TILEID','Z','Z_QN'])
        arz['TILEID'] = arz['TILEID'].astype(int)
        print(arz.dtype.names)
        #arz['TILE'].name = 'TILEID'
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

    if tp == 'QSO':
        #if azfm == 'hp':
        #    print('number of good z according to qso file '+str(len(dz)-np.sum(dz['Z_HP'].mask)))
        #    dz['Z_HP'] = dz['Z_HP'].filled(999999)
        #else:
        print('number of good z according to qso file '+str(len(dz)-np.sum(dz['Z'].mask)))
    try:
        dz['Z'] = dz['Z'].filled(999999)
    except:
        print('filling masked Z rows did not succeed')
    selm = dz['Z'] == 999999
    print('999999s for Z',len(dz[selm]))
    selo = dz['LOCATION_ASSIGNED'] == True
    print('unique Z for unassigned:')
    print(np.unique(dz[~selo]['Z']))

    print('length after cutting to unique targetid '+str(len(dz)))
    print('LOCATION_ASSIGNED numbers')
    print(np.unique(dz['LOCATION_ASSIGNED'],return_counts=True))

    print('TILELOCID_ASSIGNED numbers')
    print(np.unique(dz['TILELOCID_ASSIGNED'],return_counts=True))
    #print('length after join to file with tiles info is '+str(len(dz)))
    #NT = np.zeros(len(dz))
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
    print('number of unique targets around unassigned locations is '+str(np.sum(natloc)))

    #should not be necessary any more
#     locs = np.copy(dz['TILELOCID'])
#     print('reassigning TILELOCID for duplicates ')
#     nch = 0
#     nbl = 0
#     tlids = dz['TILELOCIDS']
#     for ii in range(0,len(dz['TILEID'])): #not sure why, but this only works when using loop for Table.read but array option works for fitsio.read
#         ti = dz[ii]['TILEID']
#
#         if natloc[ii]:# == False:
#             nbl += 1
#             s = 0
#             tids = tlids[ii].split('-')
#             if s == 0:
#                 for tl in tids:
#                     ttlocid  = int(tl)
#                     if np.isin(ttlocid,loclz):
#                         locs[ii] = ttlocid #use below instead and assign at end, maybe faster
#                         nch += 1
#                         s = 1
#                         break
#         if ii%100000 == 0:
#             print(ii,len(dz['TILEID']),ti,nch,nbl)
#
#     dz['TILELOCID'] = locs
#     locl,nlocl = np.unique(dz['TILELOCID'],return_counts=True)
#     loclz,nloclz = np.unique(dzz['TILELOCID'],return_counts=True)

    #print(np.unique(dz['NTILE']))

    print('getting fraction assigned for each tilelocid')
    nm = 0
    nmt =0
    pd = []
    nloclt = len(locl)
    lzs = np.isin(locl,loclz)
    for i in range(0,len(locl)):
        if i%100000 == 0:
            print('at row '+str(i)+' of '+str(nloclt))
        nt = nlocl[i]
        nz = lzs[i]
        loc = locl[i]
        pd.append((loc,nz/nt))
    pd = dict(pd)
    for i in range(0,len(dz)):
        probl[i] = pd[dz['TILELOCID'][i]]
    print('number of fibers with no observation, number targets on those fibers')
    print(nm,nmt)

    dz['FRACZ_TILELOCID'] = probl
    print('sum of 1/FRACZ_TILELOCID, 1/COMP_TILE, and length of input; no longer rejecting unobserved loc, so wont match')
    print(np.sum(1./dz[wz]['FRACZ_TILELOCID']),np.sum(1./dz[wz]['COMP_TILE']),len(dz))

    print(np.unique(dz['NTILE']))
    common.write_LSS(dz,outf)
    #dz.write(outf,format='fits', overwrite=True)
    #print('wrote '+outf)


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


def mkclusdat(fl,weighttileloc=True,zmask=False,tp='',dchi2=9,tsnrcut=80,rcut=None,ntilecut=0,ccut=None,ebits=None,zmin=0,zmax=6):
    '''
    fl is the root of the input/output file
    weighttileloc determines whether to include 1/FRACZ_TILELOCID as a completeness weight
    zmask determines whether to apply a mask at some given redshift
    tp is the target type
    dchi2 is the threshold for keeping as a good redshift
    tnsrcut determines where to mask based on the tsnr2 value (defined below per tracer)

    '''
    wzm = '_'
    if ccut is not None:
        wzm = ccut+'_' #you could change this to however you want the file names to turn out

    if zmask:
        wzm += 'zmask_'
    if rcut is not None:
        wzm += 'rmin'+str(rcut[0])+'rmax'+str(rcut[1])+'_'
    if ntilecut > 0:
        wzm += 'ntileg'+str(ntilecut)+'_'
    outf = fl+wzm+'clustering.dat.fits'
    ff = Table.read(fl+'_full.dat.fits')
    cols = list(ff.dtype.names)
    if 'Z' in cols:
        print('Z column already in full file')
    else:
        ff['Z_not4clus'].name = 'Z'
    if tp == 'QSO':
        #good redshifts are currently just the ones that should have been defined in the QSO file when merged in full
        wz = ff['Z']*0 == 0
        wz &= ff['Z'] != 999999
        wz &= ff['Z'] != 1.e20
        wz &= ff['ZWARN'] != 999999
        wz &= ff['TSNR2_QSO'] > tsnrcut

    if tp[:3] == 'ELG':
        ff = get_ELG_SSR_tile(ff,dchi2,tsnrcut=tsnrcut)
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
        wz = ff['ZWARN'] == 0
        wz &= ff['ZWARN']*0 == 0
        wz &= ff['ZWARN'] != 999999

        selg = ssr_tools.LRG_goodz(ff)
        #drz = (10**(3 - 3.5*ff['Z']))
        #mask_bad = (drz>30) & (ff['DELTACHI2']<30)
        #mask_bad |= (drz<30) & (ff['DELTACHI2']<drz)
        #mask_bad |= (ff['DELTACHI2']<10)
        #wz &= ff['Z']<1.4
        #wz &= (~mask_bad)
        wz &= selg

        #wz &= ff['DELTACHI2'] > dchi2
        print('length after Rongpu cut '+str(len(ff[wz])))
        wz &= ff['TSNR2_ELG'] > tsnrcut
        print('length after tsnrcut '+str(len(ff[wz])))

    if tp[:3] == 'BGS':
        wz = ff['ZWARN'] == 0
        wz &= ff['ZWARN']*0 == 0
        wz &= ff['ZWARN'] != 999999

        print('applying extra cut for BGS')
        wz &= ff['DELTACHI2'] > dchi2
        print('length after dchi2 cut '+str(len(ff[wz])))
        wz &= ff['TSNR2_BGS'] > tsnrcut
        print('length after tsnrcut '+str(len(ff[wz])))


    ff = ff[wz]
    print('length after cutting to good z '+str(len(ff)))
    ff['WEIGHT'] = np.ones(len(ff))#ff['WEIGHT_ZFAIL']
    ff['WEIGHT_ZFAIL'] = np.ones(len(ff))
    if tp[:3] == 'LRG':
        lrg = ssr_tools.LRG_ssr()
        ff = lrg.add_modpre(ff)
        ff['WEIGHT_ZFAIL'] = 1./ff['mod_success_rate']
        print('min/max of zfail weights:')
        print(np.min(ff['WEIGHT_ZFAIL']),np.max(ff['WEIGHT_ZFAIL']))

        print('checking sum of zfail weights compared to length of good z')
        print(len(ff),np.sum(ff['WEIGHT_ZFAIL']))
        ff['WEIGHT'] *= ff['WEIGHT_ZFAIL']

    if tp == 'BGS_BRIGHT':
        bgs = ssr_tools.BGS_ssr()
        ff = bgs.add_modpre(ff,fl)
        ff['WEIGHT_ZFAIL'] = np.clip(1./ff['mod_success_rate'],1,1.2)
        print('min/max of zfail weights:')
        print(np.min(ff['WEIGHT_ZFAIL']),np.max(ff['WEIGHT_ZFAIL']))
        print('checking sum of zfail weights compared to length of good z')
        print(len(ff),np.sum(ff['WEIGHT_ZFAIL']))
        ff['WEIGHT'] *= ff['WEIGHT_ZFAIL']


    if tp == 'ELG_LOP':
        elg = ssr_tools.ELG_ssr()
        ff = elg.add_modpre(ff)
        print('min/max of zfail weights:')
        print(np.min(ff['WEIGHT_ZFAIL']),np.max(ff['WEIGHT_ZFAIL']))

        print('checking sum of zfail weights compared to length of good z')
        print(len(ff),np.sum(ff['WEIGHT_ZFAIL']))
        ff['WEIGHT'] *= ff['WEIGHT_ZFAIL']

    if tp == 'QSO':
        qso = ssr_tools.QSO_ssr()
        ff = qso.add_modpre(ff,fl)
        print(np.min(ff['WEIGHT_ZFAIL']),np.max(ff['WEIGHT_ZFAIL']))
        ff['WEIGHT_ZFAIL'] = np.clip(ff['WEIGHT_ZFAIL'],1,2)
        print('min/max of zfail weights:')
        print(np.min(ff['WEIGHT_ZFAIL']),np.max(ff['WEIGHT_ZFAIL']))
        print('checking sum of zfail weights compared to length of good z')
        print(len(ff),np.sum(ff['WEIGHT_ZFAIL']))
        ff['WEIGHT'] *= ff['WEIGHT_ZFAIL']

        
    #if tp[:3] == 'ELG':
    #    ff['WEIGHT_ZFAIL'] = 1./ff['relSSR_tile']
    #    ff['WEIGHT'] *= ff['WEIGHT_ZFAIL']
    if weighttileloc == True:
        ff['WEIGHT_COMP'] = 1./ff['FRACZ_TILELOCID']
        ff['WEIGHT'] *= ff['WEIGHT_COMP']

    #weights for imaging systematic go here
    ff['WEIGHT_SYS'] =  np.ones(len(ff)) #need to initialize these at 1
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
    
    selz = ff['Z'] > zmin
    selz &= ff['Z'] < zmax
    ff = ff[selz]


    kl = ['RA','DEC','Z','WEIGHT','TARGETID','NTILE','TILES','WEIGHT_SYS','WEIGHT_COMP','WEIGHT_ZFAIL']
    if tp[:3] == 'BGS':
        #ff['flux_r_dered'] = ff['FLUX_R']/ff['MW_TRANSMISSION_R']
        #kl.append('flux_r_dered')
        #print(kl)
        fcols = ['G','R','Z','W1','W2']
        ff = common.add_dered_flux(ff,fcols)
        for col in fcols:
            kl.append('flux_'+col.lower()+'_dered')
        print(kl)
        if ccut == '-21.5':
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
        
    wn = ff['PHOTSYS'] == 'N'

    ff.keep_columns(kl)
    print('minimum,maximum weight')
    print(np.min(ff['WEIGHT']),np.max(ff['WEIGHT']))

    #comments = ["DA02 'clustering' LSS catalog for data, all regions","entries are only for data with good redshifts"]
    #common.write_LSS(ff,outf,comments)

    outfn = fl+wzm+'N_clustering.dat.fits'
    comments = ["DA02 'clustering' LSS catalog for data, BASS/MzLS region","entries are only for data with good redshifts"]
    common.write_LSS(ff[wn],outfn,comments)

    outfn = fl+wzm+'S_clustering.dat.fits'
    comments = ["DA02 'clustering' LSS catalog for data, DECaLS region","entries are only for data with good redshifts"]
    ffs = ff[~wn]
    common.write_LSS(ffs,outfn,comments)

#     for reg,com in zip(['DS','DN'],[' SGC ',' NGC ']): #split DECaLS NGC/SGC
#         outfn = fl+wzm+reg+'_clustering.dat.fits'
#         sel = densvar.sel_reg(ffs['RA'],ffs['DEC'],reg)
#         comments = ["DA02 'clustering' LSS catalog for data, DECaLS"+com+"region","entries are only for data with good redshifts"]
#         common.write_LSS(ffs[sel],outfn,comments)

def mkclusran(flin,fl,rann,rcols=['Z','WEIGHT'],zmask=False,tsnrcut=80,tsnrcol='TSNR2_ELG',ebits=None):
    #first find tilelocids where fiber was wanted, but none was assigned; should take care of all priority issues
    wzm = ''
    if zmask:
        wzm = 'zmask_'

    #ffd = Table.read(fl+'full.dat.fits')
    #fcd = Table.read(fl+wzm+'clustering.dat.fits')
    ffr = Table.read(flin+str(rann)+'_full.ran.fits')

    #if type[:3] == 'ELG' or type == 'LRG':
    wz = ffr[tsnrcol] > tsnrcut
    #wif = np.isin(ffr['TILELOCID'],ffd['TILELOCID'])
    #wic = np.isin(ffr['TILELOCID'],fcd['TILELOCID'])
    #wb = wif & ~wic #these are the tilelocid in the full but not in clustering, should be masked
    #ffc = ffr[~wb]
    ffc = ffr[wz]
    print('length after,before tsnr cut:')
    print(len(ffc),len(ffr))
    #inds = np.random.choice(len(fcd),len(ffc))
    #dshuf = fcd[inds]
    fcdn = Table.read(fl+wzm+'N_clustering.dat.fits')
    kc = ['RA','DEC','Z','WEIGHT','TARGETID','NTILE','TILES']
    rcols = np.array(rcols)
    wc = np.isin(rcols,list(fcdn.dtype.names))
    rcols = rcols[wc]
    print('columns sampled from data are:')
    print(rcols)

    #for col in rcols:
    #    ffc[col] = dshuf[col]
    #    kc.append(col)
    wn = ffc['PHOTSYS'] == 'N'

    #ffc.keep_columns(kc)
    #outf =  fl+wzm+str(rann)+'_clustering.ran.fits'
    #comments = ["DA02 'clustering' LSS catalog for random number "+str(rann)+", all regions","entries are only for data with good redshifts"]
    #common.write_LSS(ffc,outf,comments)

    outfn =  fl+wzm+'N_'+str(rann)+'_clustering.ran.fits'
    
    ffcn = ffc[wn]
    inds = np.random.choice(len(fcdn),len(ffcn))
    dshuf = fcdn[inds]
    for col in rcols:
        ffcn[col] = dshuf[col]
        kc.append(col)
    ffcn.keep_columns(kc)
    
    comments = ["DA02 'clustering' LSS catalog for random number "+str(rann)+", BASS/MzLS region","entries are only for data with good redshifts"]
    common.write_LSS(ffcn,outfn,comments)

    outfs =  fl+wzm+'S_'+str(rann)+'_clustering.ran.fits'
    fcds = Table.read(fl+wzm+'S_clustering.dat.fits')
    ffcs = ffc[~wn]
    inds = np.random.choice(len(fcds),len(ffcs))
    dshuf = fcds[inds]
    for col in rcols:
        ffcs[col] = dshuf[col]
    ffcs.keep_columns(kc)
    comments = ["DA02 'clustering' LSS catalog for random number "+str(rann)+", DECaLS region","entries are only for data with good redshifts"]
    common.write_LSS(ffcs,outfs,comments)

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

def randomtiles_allmain_pix(tiles,dirout='/global/cfs/cdirs/desi/survey/catalogs/main/LSS/random',imin=0,imax=18,dirrt='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/' ):
    '''
    tiles should be a table containing the relevant info
    '''

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
