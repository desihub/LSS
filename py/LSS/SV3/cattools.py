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
    tf.keep_columns(['TARGETID','LOCATION','FIBERSTATUS','PRIORITY','DELTA_X','DELTA_Y','PSF_TO_FIBER_SPECFLUX','EXPTIME','OBJTYPE'])
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
    tf = Table.read(coaddir+str(tile)+'/'+zdate+'/zbest-'+str(specs[0])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='FIBERMAP')
    for i in range(1,len(specs)):
        tnf = Table.read(coaddir+str(tile)+'/'+zdate+'/zbest-'+str(specs[i])+'-'+str(tile)+'-thru'+zdate+'.fits',hdu='FIBERMAP')
        tf = vstack([tf,tnf],metadata_conflicts='silent')
        
    
    tf = unique(tf,keys=['TARGETID'])
    tf.keep_columns(['TARGETID','LOCATION','FIBERSTATUS','PRIORITY','DELTA_X','DELTA_Y','PSF_TO_FIBER_SPECFLUX','EXPTIME','OBJTYPE'])
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
    
    tt = join(tt,tars,keys=['TARGETID'],table_names = ['_AVAIL', ''], uniq_col_name='{col_name}{table_name}')
    
    #tfa = unique(tfa[wgt],keys=['TARGETID'])
    #wtype = ((tt[tp] & 2**tarbit) > 0) #don't cut by type here any more
    #tt = tt[wtype]
    
    #tfa = join(tfa,tt,keys=['TARGETID'])
    #tft = join(tft,tt,keys=['TARGETID'])
    #print(str(len(tfa)) +' unique targets with good locations and  at '+str(len(np.unique(tfa['LOCATION'])))+' unique locations and '+str(len(tft))+ ' total unique targets at '+str(len(np.unique(tft['LOCATION']))) +' unique locations ')

    #Mark targets that actually got assigned fibers
    tfall = Table.read(faf,hdu='FIBERASSIGN')
    
    tfall.keep_columns(['TARGETID','LOCATION','PRIORITY'])
    
    tt = join(tt,tfall,keys=['TARGETID'],join_type='left',table_names = ['', '_ASSIGNED'], uniq_col_name='{col_name}{table_name}')
    print(tt.columns)
    #wgl = np.isin(tfa['LOCATION_ASSIGNED'],goodloc)
    #wtype = ((tfa[tp] & 2**tarbit) > 0)
    #wtfa = wgl & wtype
    #print('number of assigned fibers at good locations '+str(len(tfa[wtfa])))

    wal = tt['LOCATION_ASSIGNED']*0 == 0
    print('number of assigned fibers '+str(len(tt[wal])))
    print('number of unique target id '+str(len(np.unique(tt[wal]['TARGETID']))))
    print('max priority of assigned '+str(np.max(tt[wal]['PRIORITY_ASSIGNED'])))
    tt[wal]['LOCATION'] = tt[wal]['LOCATION_ASSIGNED']
    tt[wal]['LOCATION_AVAIL'] = tt[wal]['LOCATION_ASSIGNED']
    tt['LOCATION_ASSIGNED'] = np.zeros(len(tt),dtype=int)
    tt['LOCATION_ASSIGNED'][wal] = 1
    wal = tt['LOCATION_ASSIGNED'] == 1
    print('number of assigned fibers '+str(len(tt[wal]))+' (check to match agrees with above)')
    wal = tt['LOCATION']*0 == 0
    print('number of locations from z file '+str(len(tt[wal]))+' (check to match agrees with above)')
    #tt['PRIORITY_ASSIGNED'] = np.vectorize(pdict.__getitem__)(tt['LOCATION'])

    return tt

def combtiles(tiles,catdir,pd,tp='ALL'):
    '''
    For list of tileids, combine data generated per tile , taking care of overlaps
    
    '''

    s = 0
 
    for tile in tiles:
        fl = catdir+tp+str(tile)+'_full.dat.fits'
        fgun = Table.read(fl)
        aa = np.chararray(len(fgun),unicode=True,itemsize=100)
        aa[:] = str(tile)
        fgun['TILE'] = aa
        fgun['TILELOCID'] = 10000*tile +fgun['LOCATION_AVAIL']
        if s == 0:
            fgu = fgun
            s =1
        #wm = np.ma.getmaskarray(fgun['LOCATION_ASSIGNED'])
        #fgun['TILELOCID_ASSIGNED'] = 0
        #fgun['TILELOCID_ASSIGNED'][~wm] = tile*10000+fgun['LOCATION_ASSIGNED'][~wm]
        else:
            fgu = vstack([fgu,fgun])

    print(len(np.unique(fgu['TARGETID'])),np.sum(fgu['LOCATION_ASSIGNED']))
    
    wn = fgu['PRIORITY_ASSIGNED']*0 != 0
    wn |= fgu['PRIORITY_ASSIGNED'] == 999999
    print(len(fgu[~wn]),np.max(fgu[~wn]['PRIORITY_ASSIGNED']),'max priority assigned')
    fgu[wn]['PRIORITY_ASSIGNED'] = 0
    fgu['sort'] = -1.*fgu['LOCATION_ASSIGNED']*fgu['PRIORITY_ASSIGNED'] #create this column so assigned always show up in order of highest priority
    wa = fgu['LOCATION_ASSIGNED'] == 1
    #wa &= fgu['PRIORITY_ASSIGNED'] >= 2000 #this was put SV2 to ignore BGS repeats
    fa = fgu[wa]
    print(len(fa),len(np.unique(fa['TARGETID'])))
    fgu.sort('sort')
    fu = unique(fgu,keys=['TARGETID'])
    print(np.sum(fu['LOCATION_ASSIGNED']))
    tidsu = fu['TARGETID']
    tids = fgu['TARGETID']
    tiles = fgu['TILE']
    tilesu = fu['TILE']
    for ii in range(0,len(tidsu)): #this takes a long time and something more efficient will be necessary
        tid = tidsu[ii]#fu[ii]['TARGETID']
        wt = tids == tid
        ot = tilesu[ii]
        
        tt = tiles[wt]
        for tl in tt:
            if tl != ot:
                tilesu[ii] += '-'+str(tl)
        if ii%1000 == 0:
            print(ii)        
    fu['TILE'] = tilesu
    print(np.unique(fu['TILE']))
    #wa = fu['LOCATION_ASSIGNED'] == 1
    #wa &= fu['PRIORITY_ASSIGNED'] >= 2000
    print(np.sum(fu['LOCATION_ASSIGNED']))
    fu.write(catdir+tp+'Alltiles_'+pd+'_full.dat.fits',format='fits', overwrite=True)    

def countloc(aa):
    locs = aa['LOCATION']
    locsa = aa['LOCATION_ASSIGNED']
    la = np.max(locs)+1
    nl = np.zeros(la)
    nla = np.zeros(la)
    for i in range(0,len(aa)):
        nl[locs[i]] += 1
        nla[locs[i]] += locsa[i]
    return nl,nla


def combran(tiles,rann,randir,ddir,tp,tmask,tc='SV3_DESI_TARGET',maskzfail=True):

    s = 0
    #tiles.sort('ZDATE')
    for tile,zdate in zip(tiles['TILEID'],tiles['ZDATE']):
        #tspec = combfibmap(tile,zdate)
        #pdict,gloc = goodlocdict(tspec)
        dt = ddir+'ALL'+str(tile)+'_full.dat.fits'
        ffa = randir+str(rann)+'/fba-'+str(tile).zfill(6)+'.fits'
        ffna = randir+str(rann)+'/tilenofa-'+str(tile)+'.fits'
        if os.path.isfile(ffa):
            fd = Table.read(dt)
            gloc = np.unique(fd['LOCATION']) #bad locations already removed from this files
            wt = (dt[tc] & tmask[tp]) > 0
            fd = fd[wt]
            wzf = fd['ZWARN'] != 0 
            wzf &= fd['ZWARN'] != 999999
            wzf &= fd['ZWARN']*0 == 0
            loc_fail = np.unique(fd[wzf]['LOCATION'])
            nl,nla = countloc(fd)
        # 
            #find the locations that were requested by LRGs but not assigned
            locsna = []
            for i in range(0,len(nla)):
                if nla[i] == 0 and nl[i] > 0:
                    locsna.append(i)

            fa = Table.read(ffa,hdu='FAVAIL')
            wg = np.isin(fa['LOCATION'],gloc)
            wg &= ~np.isin(fa['LOCATION'],locsna)
            if maskzfail:
                wg &= np.isin(fa['LOCATION'],loc_fail)
            
            #wzt = wpr & ~wzf & ~wna

            fg = fa[wg]
            print('before,after vetoing locations:')
            print(len(fa),len(fg))
            fgun = unique(fg,keys=['TARGETID'])
            ffna = Table.read(ffna)
            fgun = join(fgun,ffna,keys=['TARGETID'])
            print(str(len(fgun))+' unique new randoms')
            aa = np.chararray(len(fgun),unicode=True,itemsize=100)
            aa[:] = str(tile)
            fgun['TILE'] = aa
            fgun['TILELOCID'] = 10000*tile +fgun['LOCATION']
            if s == 0:
                fgu = fgun
                s = 1
            else:   
                fv = vstack([fgu,fgun])
                fgo = fgu.copy()
                fgu = unique(fv,keys='TARGETID')#,keep='last') 
                
                dids = np.isin(fgun['TARGETID'],fgo['TARGETID']) #get the rows with target IDs that were duplicates in the new file
                didsc = np.isin(fgu['TARGETID'],fgun['TARGETID'][dids]) #get the row in the concatenated table that had dup IDs
                #print(len(fgu),len(fgo),len(fgun),len(fgu[didsc]),len(fgun[dids]))
                fgu['TILELOCID'][didsc] = fgun['TILELOCID'][dids] #give the repeats the new tilelocids, since those are the most likely to be available to low priority targets

                aa = np.chararray(len(fgu['TILE']),unicode=True,itemsize=20)
                aa[:] = '-'+str(tile)
                #rint(aa)
                ms = np.core.defchararray.add(fgu['TILE'][didsc],aa[didsc])
                #print(ms)
                fgu['TILE'][didsc] = ms #add the tile info
                print(str(len(fgu))+' unique total randoms')
    fgu.write(randir+str(rann)+'/rancomb_'+tp+'_Alltiles.fits',format='fits', overwrite=True)

def mkfullran(randir,rann,imbits,outf,tp):
    zf = randir+str(rann)+'/rancomb_'+tp+'_Alltiles.fits'
    dz = Table.read(zf)
    
    dz = cutphotmask(dz,imbits)
    NT = np.zeros(len(dz))
    for ii in range(0,len(dz['TILE'])): #not sure why, but this only works when using loop for Table.read but array option works for fitsio.read
        NT[ii] = np.char.count(dz['TILE'][ii],'-')+1
    
    #NT = np.char.count(dz['TILE'],'-')
    #NT += 1
    print(np.unique(NT))
    dz['NTILE'] = NT
    dz.write(outf,format='fits', overwrite=True)
    


def mkfulldat(zf,imbits,tdir,tp,bit,outf):
    #from desitarget.mtl import inflate_ledger
    dz = Table.read(zf) 
    wtype = ((dz[tp] & bit) > 0)
    dz = dz[wtype]
    dz = cutphotmask(dz,imbits)
    
    NT = np.zeros(len(dz))
    for ii in range(0,len(dz['TILE'])): #not sure why, but this only works when using loop for Table.read but array option works for fitsio.read
        NT[ii] = np.char.count(dz['TILE'][ii],'-')+1
    #NT = np.char.count(dz['TILE'],'-')
    #NT += 1
    print(np.unique(NT))
    wz = dz['ZWARN'] == 0
    dzz = dz[wz]
    probl = np.zeros(len(dz))
    #dr = fitsio.read(e2eout+ program+'/'+type+'_oneper_full.ran.fits')
    locl,nlocl = np.unique(dz['TILELOCID'],return_counts=True)
    wa = dzz['LOCATION_ASSIGNED'] == 1
    if len(dzz[wa]) != len(dzz):
        print('!found some zwarn = 0 without location_assigned = 1!')
    loclz,nloclz = np.unique(dzz['TILELOCID'],return_counts=True)
    print(np.max(nloclz),np.min(loclz))
    print(len(locl),len(nloclz))
    nm = 0
    nmt =0
    pd = []
    for i in range(0,len(locl)):
        if i%10000 == 0:
            print('at row '+str(i))
        nt = nlocl[i]
        loc = locl[i]
        w = loclz == loc
        nz = 0
        if len(loclz[w]) == 1:
            nz = 1.#nloclz[w]
            
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
    print('number of fibers with no good z, number targets on those fibers')
    print(nm,nmt)
    #print(np.min(probl),np.max(probl))
    #dz = Table.read(zf) #table is slow, so using fitsio above, Table here
    dz['FRACZ_TILELOCID'] = probl
    #print(np.unique(dz['TILE']))
    dz['NTILE']  = NT
    print(np.unique(dz['NTILE']))
    dz.write(outf,format='fits', overwrite=True)

def mkclusdat(fl,weighttileloc=True):
    '''
    take full catalog, cut to ra,dec,z add any weight
    program is dark,gray, or bright
    type is 'LRG', 'QSO', 'ELG', or 'BGS'

    '''    
    ff = Table.read(fl+'full.dat.fits')
    outf = fl+'clustering.dat.fits'
    wz = ff['ZWARN'] == 0
    ff = ff[wz]
    ff['WEIGHT'] = np.ones(len(ff))
    if weighttileloc == True:
        ff['WEIGHT'] = 1./ff['FRACZ_TILELOCID']

    ff.keep_columns(['RA','DEC','Z','WEIGHT','TARGETID','NTILE','TILELOCID'])
    print('minimum,maximum weight')
    print(np.min(ff['WEIGHT']),np.max(ff['WEIGHT']))


    ff.write(outf,format='fits', overwrite=True)

def mkclusran(fl,rann,rcols=['Z','WEIGHT']):
    #first find tilelocids where fiber was wanted, but none was assigned; should take care of all priority issues
    ffd = Table.read(fl+'full.dat.fits')
    fcd = Table.read(fl+'clustering.dat.fits')
    ffr = Table.read(fl+str(rann)+'_full.ran.fits')
    #wif = np.isin(ffr['TILELOCID'],ffd['TILELOCID'])
    #wic = np.isin(ffr['TILELOCID'],fcd['TILELOCID'])
    #wb = wif & ~wic #these are the tilelocid in the full but not in clustering, should be masked
    #ffc = ffr[~wb]
    ffc = ffr
    print(len(ffc),len(ffr))
    inds = np.random.choice(len(fcd),len(ffc))
    dshuf = fcd[inds]

    for col in rcols: 
        ffc[col] = dshuf[col] 
    ffc.keep_columns(['RA','DEC','Z','WEIGHT','TARGETID','NTILE','TILELOCID'])  
    outf =  fl+str(rann)+'_clustering.ran.fits' 
    ffc.write(outf,format='fits', overwrite=True)

    

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
