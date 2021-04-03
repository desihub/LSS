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
import desimodel.footprint
import desimodel.focalplane
from random import random
from desitarget.io import read_targets_in_tiles

from LSS.Cosmo import distance

def combspecdata(tile,zdate,coaddir='/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/' ):
    #put data from different spectrographs together, one table for fibermap, other for z
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
    for i in range(1,len(specs)):
        tn = Table.read(coaddir+str(tile)+'/'+zdate+'/zbest-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='ZBEST')
        tnf = Table.read(coaddir+str(tile)+'/'+zdate+'/zbest-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='FIBERMAP')
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


def gettarinfo_type(faf,tars,goodloc,tarbit,pdict,tp='SV2_DESI_TARGET'):
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

    tt = unique(tt,keys=['TARGETID']) #cut to unique target ids
    
    wgt = (np.isin(tt['LOCATION'],goodloc)) 
    print(str(len(np.unique(tt[wgt]['LOCATION']))) + ' good locations')
    print('comparison of number targets, number of targets with good locations')
    print(len(tt),len(tt[wgt]))
    
    tt = tt[wgt]
    #print(tarf)
    #tars = Table.read(tarf)
    #tars.remove_columns(['Z','ZWARN'])#,'PRIORITY','SUBPRIORITY','OBSCONDITIONS'])
    #we want to get these from the zbest file that is specific to the tile and thus when it was observed
    tars = tars[[b for b in list(tars.dtype.names) if b != 'Z']]
    tars = tars[[b for b in list(tars.dtype.names) if b != 'ZWARN']]
    tars = tars[[b for b in list(tars.dtype.names) if b != 'PRIORITY']]
    
    tt = join(tt,tars,keys=['TARGETID'])
    
    #tfa = unique(tfa[wgt],keys=['TARGETID'])
    wtype = ((tt[tp] & 2**tarbit) > 0)
    tt = tt[wtype]
    
    #tfa = join(tfa,tt,keys=['TARGETID'])
    #tft = join(tft,tt,keys=['TARGETID'])
    #print(str(len(tfa)) +' unique targets with good locations and  at '+str(len(np.unique(tfa['LOCATION'])))+' unique locations and '+str(len(tft))+ ' total unique targets at '+str(len(np.unique(tft['LOCATION']))) +' unique locations ')

    #Mark targets that actually got assigned fibers
    tfall = Table.read(faf,hdu='FIBERASSIGN')
    
    tfall.keep_columns(['TARGETID','LOCATION'])
    
    tt = join(tt,tfall,keys=['TARGETID'],join_type='left',table_names = ['', '_ASSIGNED'], uniq_col_name='{col_name}{table_name}')
    
    #wgl = np.isin(tfa['LOCATION_ASSIGNED'],goodloc)
    #wtype = ((tfa[tp] & 2**tarbit) > 0)
    #wtfa = wgl & wtype
    #print('number of assigned fibers at good locations '+str(len(tfa[wtfa])))

    wal = tt['LOCATION_ASSIGNED']*0 == 0
    print('number of assigned fibers '+str(len(tt[wal])))
    tt['LOCATION_ASSIGNED'] = np.zeros(len(tt),dtype=int)
    tt['LOCATION_ASSIGNED'][wal] = 1
    wal = tt['LOCATION_ASSIGNED'] == 1
    print('number of assigned fibers '+str(len(tt[wal]))+' (check to match agrees with above)')
    tt['PRIORITY_ASSIGNED'] = np.vectorize(pdict.__getitem__)(tt['LOCATION'])

    return tt

def combtiles(tiles,catdir,tp):
    '''
    For list of tileids, combine, taking care of overlaps
    
    '''

    s = 0
 
    for tile in tiles:
        fl = catdir+tp+str(tile)+'_full.dat.fits'
        fgun = Table.read(fl)
        aa = np.chararray(len(fgun),unicode=True,itemsize=100)
        aa[:] = str(tile)
        fgun['TILE'] = aa
        fgun['TILELOCID'] = 10000*tile +fgun['LOCATION']
        if s == 0:
            fgu = fgun
            s =1
        #wm = np.ma.getmaskarray(fgun['LOCATION_ASSIGNED'])
        #fgun['TILELOCID_ASSIGNED'] = 0
        #fgun['TILELOCID_ASSIGNED'][~wm] = tile*10000+fgun['LOCATION_ASSIGNED'][~wm]
        else:
            fgu = vstack([fgu,fgun])

    print(len(np.unique(fgu['TARGETID'])),np.sum(fgu['LOCATION_ASSIGNED']))
    #fgu.write(e2eout+outf,format='fits', overwrite=True)    


def randomtiles_allSV2(tiles,dirout='/global/cfs/cdirs/desi/survey/catalogs/SV2/LSS/random',imin=0,imax=18,dirr='/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/'):
    '''
    tiles should be a table containing the relevant info
    '''
    trad = desimodel.focalplane.get_tile_radius_deg()*1.1 #make 10% greater just in case
    print(trad)
    for ii in range(imin,imax):
        rt = fitsio.read(dirr+'/randoms-1-'+str(ii)+'.fits',columns=['RA','DEC','TARGETID','PHOTSYS','NOBS_G','NOBS_R','NOBS_Z','MASKBITS'])
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
                rmtl['SV1_DESI_TARGET'] = np.ones(len(rmtl),dtype=int)*2
                rmtl['NUMOBS_INIT'] = np.zeros(len(rmtl),dtype=int)
                rmtl['NUMOBS_MORE'] = np.ones(len(rmtl),dtype=int)
                rmtl['PRIORITY'] = np.ones(len(rmtl),dtype=int)*3400
                rmtl['OBSCONDITIONS'] = np.ones(len(rmtl),dtype=int)*516#tiles['OBSCONDITIONS'][i]
                rmtl['SUBPRIORITY'] = np.random.random(len(rmtl))
                rmtl.write(fname,format='fits', overwrite=True)
                print('added columns, wrote to '+fname)
