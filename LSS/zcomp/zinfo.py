from astropy.table import Table,join,vstack
import numpy as np
import fitsio
# AR adding some imports
from desitarget.sv1 import sv1_targetmask
import os
from glob import glob
from astropy.io import fits
import fitsio

#'test'

def get_zfits(fl,tid):#tile,specn,subset,tid,zfitdir='blanc'):
    #fl = zfitdir+'/redrock-'+str(specn)+'-'+str(tile)+'-'+subset+'.h5'
    pt = 'zfit/'+str(tid)+'/zfit'
    zfits = Table.read(fl,path=pt)
    return zfits

def comb_subset_vert(tarbit,tp,subsets,tile,coaddir,exposures,outf,tt,mfn='temp.txt',md='rel'):
    '''
    performs a vertical concatenation of the data for a tile, so each targetid shows up N_subset times
    subsets is a list of the subsets (strings)
    tile is the particular tile (string)
    coaddir is where the data comes from (string; e.g., the directory pointing to the Blanc release)
    exposures is the file containing the information per exposure, used to get depth information (data array read by, e.g., fitsio)
    outf is where the fits file gets written out (string)
    '''
    ss = 0 #use to switch from creating to concatenating
    for night in subsets:
        #print(tile,night)
        #if tile == '80607' and night == 'subset-1':
        #    pass

        if len(night) > 0:
            coaddiru = coaddir
            if md == 'rel':
                coaddiru = coaddir+'/'+night
                tsmd='fm'
            if md == 'RZ':
                tsmd = 'cf'
            tspec = get_subset(tarbit,tp,night,tile,coaddiru,exposures,mfn=mfn,tsmd=tsmd)
            if tspec is not None:
                if ss == 0:
                    tspect = tspec
                    ss = 1
                else:
                    tspect = vstack([tspect,tspec], metadata_conflicts='silent')
                print('there are now '+str(len(tspect)) +' entries with '+str(len(np.unique(tspect['TARGETID'])))+' unique target IDs')    
                    
    if ss == 1:
        tspect.sort('TARGETID')
        tspect['TARGETS'] = tt
        tspect.write(outf,format='fits', overwrite=True) 
        print('wrote to '+outf)
        return True
    else:
        print('no data for tile '+tile)
        return False

def comb_subset_vert_denali(tarbit,tp,tile,exposures,outf,tt,mfn='temp.txt',md='rel'):
    '''
    performs a vertical concatenation of the data for a tile, so each targetid shows up N_subset times
    subsets is a list of the subsets (strings)
    tile is the particular tile (string)
    coaddir is where the data comes from (string; e.g., the directory pointing to the Blanc release)
    exposures is the file containing the information per exposure, used to get depth information (data array read by, e.g., fitsio)
    outf is where the fits file gets written out (string)
    '''
    tsmd='fm'
    ss = 0 #use to switch from creating to concatenating
    ctypes = ['cumulative','perexp','pernight']
    for ct in ctypes:
        if ct == 'cumulative':
            dirs = [ f.path for f in os.scandir('/global/cfs/cdirs/desi/spectro/redux/denali/tiles/cumulative/'+tile) if f.is_dir() ][0]
            
            night = dirs[-8:]
            nights = ['thru'+night]
            dirs = [dirs]
        if ct == 'perexp':
            dirs = [ f.path for f in os.scandir('/global/cfs/cdirs/desi/spectro/redux/denali/tiles/perexp/'+tile) if f.is_dir() ]
            nights = []
            for dr in dirs:
                nights.append('exp'+dr[-8:])
        if ct == 'pernight':
            dirs = [ f.path for f in os.scandir('/global/cfs/cdirs/desi/spectro/redux/denali/tiles/pernight/'+tile) if f.is_dir() ]
            nights = []
            for dr in dirs:
                nights.append(dr[-8:])

        for dr,night in zip(dirs,nights):
            tspec = get_subset_denali(tarbit,tp,night,tile,dr,exposures,ct,mfn=mfn,tsmd=tsmd)
            if tspec is not None:
                if ss == 0:
                    tspect = tspec
                    ss = 1
                else:
                    tspect = vstack([tspect,tspec], metadata_conflicts='silent')
                print('there are now '+str(len(tspect)) +' entries with '+str(len(np.unique(tspect['TARGETID'])))+' unique target IDs')    
                    
    if ss == 1:
        tspect.sort('TARGETID')
        tspect['TARGETS'] = tt
        tspect.write(outf,format='fits', overwrite=True) 
        print('wrote to '+outf)
        return True
    else:
        print('no data for tile '+tile)
        return False


def comb_exps_vert(tarbit,tp,tile,coaddir,exposures,outf,dirout):
    '''
    performs a vertical concatenation of the exposure data for a tile, so each targetid shows up N_subset times
    tile is the particular tile (string)
    coaddir is where the data comes from (string; e.g., the directory pointing to the Blanc release)
    exposures is the file containing the information per exposure, used to get depth information (data array read by, e.g., fitsio)
    outf is where the fits file gets written out (string)
    '''
    ss = 0 #use to switch from creating to concatenating
    tid = int(tile)
    wt = exposures['TILEID'] == tid
    exps = np.unique(exposures[wt]['EXPID'])
    for exp in exps:
        tspec = get_exp(tarbit,tp,exp,tile,coaddir,exposures)
        if tspec is not None:
            info = exposures[exposures['EXPID'] == exp]
            for name in exposures.dtype.names:
                tspec[name] = info[name][0]

            if ss == 0:
                tspect = tspec
                ss = 1
            else:
                tspect = vstack([tspect,tspec], metadata_conflicts='silent')
            print('there are now '+str(len(tspect)) +' entries with '+str(len(np.unique(tspect['TARGETID'])))+' unique target IDs')    
                    
    if ss == 1:
        tspect.sort('TARGETID')
        #tspect['TARGETS'] = tt
        deepf = Table.read(outf+'.fits')
        wd = deepf['subset'] == 'deep'
        deepf = deepf[wd]
        deepf.keep_columns(['TARGETID','Z','ZWARN','DELTACHI2'])
        for name in ['Z','ZWARN','DELTACHI2']:
            deepf.rename_column(name,name+'_deep')
        tspect = join(tspect,deepf,keys=['TARGETID'],join_type='left')
        outf += '_1exp.fits'
        tspect.write(outf,format='fits', overwrite=True) 
        print('wrote to '+outf)
        return True
    else:
        print('no data for tile '+tile)
        return False


def get_tsnrinfo(exps,spec,tsnrdir='/global/cscratch1/sd/mjwilson/desi/tsnr/blanc/exptables/v0'):
    
    es = []
    bs = []
    qs = []
    ls = []
    bands = ['b','r','z']
            
    for exp in exps:
        esv = 0
        bsv = 0
        lsv = 0
        qsv = 0
        for band in bands:
            cinfo = fitsio.read(tsnrdir+'/summary_'+band+str(spec)+'.fits')
            info = cinfo[cinfo['EXPID'] == '000'+str(exp)]    
            if len(info) == 0:
                print('did not find tsnr info for band '+band+' for expid '+str(exp))
                print(tsnrdir+'/summary_'+band+str(spec)+'.fits')
                return None
            esv += info['ELGTSNR'][0] #just get total across bands per exposure
            bsv += info['BGSTSNR'][0]
            lsv += info['LRGTSNR'][0]
            qsv += info['QSOTSNR'][0]
        es.append(esv)    
        bs.append(bsv)
        ls.append(lsv) 
        qs.append(qsv)
    return es,bs,ls,qs        

def get_subset(tarbit,tp,night,tile,coaddir,exposures,mfn='temp.txt',rel='cascades',tsmd='fm'):

    print('going through subset '+night)
    if tsmd != 'fm':
        cfdir = '/global/cfs/cdirs/desi/spectro/redux/'+rel+'/tiles/'+str(tile) #defined here to allow flexibility dealing with Rongpu's files
        print('using release '+str(rel)+' cframes, is that what you want?')
        cams = ['b','r','z']
        tsnrcols = ['TSNR2_ELG','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG']
    specs = []
    #find out which spectrograph have data
    for si in range(0,10):
        try:
            fl = coaddir+'/zbest-'+str(si)+'-'+str(tile)+'-'+night+'.fits'
            fitsio.read(fl)
            fl = coaddir+'/coadd-'+str(si)+'-'+str(tile)+'-'+night+'.fits'
            fitsio.read(fl)
            specs.append(si)
        except:
            #print(fl,specs,si)
            print('no spectrograph and/or coadd '+str(si)+ ' on subset '+night)
    if len(specs) > 2: #basically required just to reject the one night with data from only 2 specs that was in exposures
        tspec = Table.read(coaddir+'/zbest-'+str(specs[0])+'-'+str(tile)+'-'+night+'.fits',hdu='ZBEST')
        tf = Table.read(coaddir+'/coadd-'+str(specs[0])+'-'+str(tile)+'-'+night+'.fits',hdu='FIBERMAP')
        if tsmd == 'fm':
            ts = Table.read(coaddir+'/coadd-'+str(specs[0])+'-'+str(tile)+'-'+night+'.fits',hdu='SCORES')
            ts.keep_columns(['TARGETID','TSNR2_ELG','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
        #this is all to get the effective coadded exposure depth; should eventually just be in the fibermap hdu
        zfm = Table.read(coaddir+'/zbest-'+str(specs[0])+'-'+str(tile)+'-'+night+'.fits',hdu='FIBERMAP')
        exps = np.unique(zfm['EXPID'])
        bd = []
        rd = []
        zd = []
        etdark = []
        etbright = []
        etback = []
        
        bda = []
        rda = []
        zda = []
        ce = 0
        for exp in exps:
            info = exposures[exposures['EXPID'] == exp]
            if len(info) == 0:
                print('did not find info for expid '+str(exp))
                fo = open(mfn,'a')
                fo.write(str(exp)+'\n')
                return None
            else:    
                
                nt = str(info['NIGHT'][0])
                print(info['B_DEPTH'],nt)
                bd.append(info['B_DEPTH'][0])
                rd.append(info['R_DEPTH'][0])
                zd.append(info['Z_DEPTH'][0]) 
                bda.append(info['B_DEPTH_EBVAIR'][0])
                rda.append(info['R_DEPTH_EBVAIR'][0])
                zda.append(info['Z_DEPTH_EBVAIR'][0]) 
                etdark.append(info['EFFTIME_DARK'][0])
                etbright.append(info['EFFTIME_BRIGHT'][0])
                etback.append(info['EFFTIME_BACKUP'][0]) 
                
                if tsmd != 'fm':
                    for cam in cams:
                        tcols  =[]
                        for col in tsnrcols:
                            tcols.append(col+'_'+cam.upper())
                        cf = Table.read(cfdir+'/'+nt+'/cframe-'+cam+str(specs[0])+'-'+str(exp).zfill(8)+'.fits',hdu='SCORES')
                        cf.keep_columns(tcols)
                        for col in tcols:
                            cf.rename_column(col, col[:-2])
                        if ce ==0 and cam == 'b':
                            tids = Table.read(cfdir+'/'+nt+'/cframe-'+cam+str(specs[0])+'-'+str(exp).zfill(8)+'.fits',hdu='FIBERMAP')
                            tnsrt = cf.copy()
                            tnsrt['TARGETID'] = tids['TARGETID']
                        else:
                            for col in tsnrcols:
                                tnsrt[col] += cf[col]         
                ce += 1                 
        #tvs = get_tsnrinfo(exps,specs[0])    
        #if tvs is not None:
        #    es,bs,ls,qs = tvs
        #else:
        #    return tvs

       
        bdt = np.zeros(500)
        rdt = np.zeros(500)
        zdt = np.zeros(500)
        etdarkt = np.zeros(500)
        etbrightt = np.zeros(500)
        etbackt = np.zeros(500)

        bdta = np.zeros(500)
        rdta = np.zeros(500)
        zdta = np.zeros(500)
        tid = zfm[0:500]['TARGETID']
        #est = np.zeros(500)
        #bst = np.zeros(500)
        #lst = np.zeros(500)
        #qst = np.zeros(500)
        for i in range(0,len(exps)):
            sel = zfm[i*500:(i+1)*500]
            w = sel['FIBERSTATUS'] == 0
            bdt[w] += bd[i]
            rdt[w] += rd[i]
            zdt[w] += zd[i]
            bdta[w] += bda[i]
            rdta[w] += rda[i]
            zdta[w] += zda[i]
            etdarkt[w] += etdark[i]
            etbrightt[w] += etbright[i]
            etbackt[w] += etback[i]
            #est[w] += es[i]
            #bst[w] += bs[i]
            #lst[w] += ls[i]
            #qst[w] += qs[i]
    
        
        tf['EXPS'] = ",".join(exps.astype(str))
        for i in range(1,len(specs)):
            zfm = Table.read(coaddir+'/zbest-'+str(specs[i])+'-'+str(tile)+'-'+night+'.fits',hdu='FIBERMAP')
            exps = np.unique(zfm['EXPID'])
            if tsmd == 'fm':
                tsn = Table.read(coaddir+'/coadd-'+str(specs[i])+'-'+str(tile)+'-'+night+'.fits',hdu='SCORES')
                tsn.keep_columns(['TARGETID','TSNR2_ELG','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
                ts = vstack([ts,tsn], metadata_conflicts='silent')
           

            tn = Table.read(coaddir+'/zbest-'+str(specs[i])+'-'+str(tile)+'-'+night+'.fits',hdu='ZBEST')
            tnf = Table.read(coaddir+'/coadd-'+str(specs[i])+'-'+str(tile)+'-'+night+'.fits',hdu='FIBERMAP')  
            tnf['EXPS'] = ",".join(exps.astype(str))
            tspec = vstack([tspec,tn], metadata_conflicts='silent')                      
            tf = vstack([tf,tnf], metadata_conflicts='silent')
            bd = []
            rd = []
            zd = []
            etdark = []
            etbright = []
            etback = []
            bda = []
            rda = []
            zda = []
            ce = 0
            for exp in exps:
                info = exposures[exposures['EXPID'] == exp]
                if len(info) == 0:
                    print('did not find info for expid '+str(exp))
                    return None
                else:
                    nt = str(info['NIGHT'][0])
                    bd.append(info['B_DEPTH'][0])
                    rd.append(info['R_DEPTH'][0])
                    zd.append(info['Z_DEPTH'][0])        
                    bda.append(info['B_DEPTH_EBVAIR'][0])
                    rda.append(info['R_DEPTH_EBVAIR'][0])
                    zda.append(info['Z_DEPTH_EBVAIR'][0])    
                    etdark.append(info['EFFTIME_DARK'][0])
                    etbright.append(info['EFFTIME_BRIGHT'][0])
                    etback.append(info['EFFTIME_BACKUP'][0]) 
                    
                        
                    if tsmd != 'fm':
                        for cam in cams:
                            tcols = []
                            for col in tsnrcols:
                                tcols.append(col+'_'+cam.upper())

                            try:
                                cf = Table.read(cfdir+'/'+nt+'/cframe-'+cam+str(specs[i])+'-'+str(exp).zfill(8)+'.fits',hdu='SCORES')
                            except:
                                print('couldnt open '+cfdir+'/'+nt+'/cframe-'+cam+str(specs[i])+'-'+str(exp).zfill(8)+'.fits')
                                return None
                            cf.keep_columns(tcols)
                            for col in tcols:
                                cf.rename_column(col, col[:-2])

                            #print(ce,cam)
                            if ce ==0 and cam == 'b':
                                tids = Table.read(cfdir+'/'+nt+'/cframe-'+cam+str(specs[i])+'-'+str(exp).zfill(8)+'.fits',hdu='FIBERMAP')
                                tsnrtn = cf.copy()
                                tsnrtn['TARGETID'] = tids['TARGETID']
                            else:
                                for col in tsnrcols:
                                    tsnrtn[col] += cf[col]         
                    ce += 1                 

            
            if tsmd != 'fm':
                tnsrt = vstack([tnsrt,tsnrtn], metadata_conflicts='silent')
            #tvs = get_tsnrinfo(exps,specs[i]) 
            #if tvs is not None:
            #    es,bs,ls,qs = tvs
            #else:
            #    return tvs
            bdtn = np.zeros(500)
            rdtn = np.zeros(500)
            zdtn = np.zeros(500)
            bdtna = np.zeros(500)
            rdtna = np.zeros(500)
            zdtna = np.zeros(500)
            etdarktn = np.zeros(500)
            etbrighttn = np.zeros(500)
            etbacktn = np.zeros(500)

            #estn = np.zeros(500)
            #bstn = np.zeros(500)
            #lstn = np.zeros(500)
            #qstn = np.zeros(500)

            tidn = zfm[0:500]['TARGETID']
            for ii in range(0,len(exps)):
                sel = zfm[ii*500:(ii+1)*500]
                w = sel['FIBERSTATUS'] == 0
                bdtn[w] += bd[ii]
                rdtn[w] += rd[ii]
                zdtn[w] += zd[ii]
                bdtna[w] += bda[ii]
                rdtna[w] += rda[ii]
                zdtna[w] += zda[ii]
                etdarktn[w] += etdark[ii]
                etbrighttn[w] += etbright[ii]
                etbacktn[w] += etback[ii]
                #estn[w] += es[ii]
                #bstn[w] += bs[ii]
                #lstn[w] += ls[ii]
                #qstn[w] += qs[ii]



            bdt = np.concatenate([bdt,bdtn])
            rdt = np.concatenate([rdt,rdtn])
            zdt = np.concatenate([zdt,zdtn])   
            bdta = np.concatenate([bdta,bdtna])
            rdta = np.concatenate([rdta,rdtna])
            zdta = np.concatenate([zdta,zdtna])   
            etdarkt = np.concatenate([etdarkt,etdarktn])
            etbrightt = np.concatenate([etbrightt,etbrighttn])
            etbackt = np.concatenate([etbackt,etbacktn])   

            #est = np.concatenate([est,estn])
            #lst = np.concatenate([lst,lstn])
            #qst = np.concatenate([qst,qstn])   
            #bst = np.concatenate([bst,bstn])

            tid = np.concatenate([tid,tidn])
            #print(np.min(rdtn),np.max(rdtn)) 
            #print(np.min(rdt),np.max(rdt)) 
        tspec = join(tspec,tf,keys=['TARGETID'], metadata_conflicts='silent')
        if tsmd != 'fm':
            tspec = join(tspec,tnsrt,keys=['TARGETID'], metadata_conflicts='silent')
        else:
            tspec = join(tspec,ts,keys=['TARGETID'], metadata_conflicts='silent')
        td = Table([bdt,rdt,zdt,bdta,rdta,zdta,etdarkt,etbrightt,etbackt,tid],names=('B_DEPTH','R_DEPTH','Z_DEPTH','B_DEPTH_EBVAIR','R_DEPTH_EBVAIR','Z_DEPTH_EBVAIR','EFFTIME_DARK','EFFTIME_BRIGHT','EFFTIME_BACK','TARGETID'))#,'ELGTSNR','BGSTSNR','LRGTSNR','QSOTSNR'
        tspec = join(tspec,td,keys=['TARGETID'], metadata_conflicts='silent')
        if tarbit != -1:
            wtype = ((tspec[tp] & 2**tarbit) > 0)
            print(str(len(tspec))+' total entries '+str(len(tspec[wtype]))+' that are requested type entries with '+str(len(np.unique(tspec[wtype]['TARGETID'])))+' unique target IDs')
            if len(tspec[wtype]) == 0:
                return None

            tspec = tspec[wtype]
        tspec['subset'] = night
        # AR adding a weight for ELGs in the QSO+ELG and QSO+LRG tiles
        # AR to down-weight QSOs which are at a higher priority
        # AR we "rescale" to the ELGxQSO/ELG ratio of the parent input target sample per tile
        if tarbit != -1:
            tspec['elgqso_weight'] = get_elgqso_weight(tarbit,tp,tile,tspec[tp])
        return tspec
    return None    

def get_subset_denali(tarbit,tp,night,tile,coaddir,exposures,ct,mfn='temp.txt',tsmd='fm'):

    print('going through subset '+night)
    if tsmd != 'fm':
        print('tsmd needs to be fm, other options not supported, all will fail!!!')
        return None
    specs = []
    #find out which spectrograph have data
    for si in range(0,10):
        try:
            fl = coaddir+'/zbest-'+str(si)+'-'+str(tile)+'-'+night+'.fits'
            fitsio.read(fl)
            fl = coaddir+'/coadd-'+str(si)+'-'+str(tile)+'-'+night+'.fits'
            fitsio.read(fl)
            specs.append(si)
        except:
            print(fl,specs,si)
            print('no spectrograph and/or coadd '+str(si)+ ' on subset '+night)
    if len(specs) > 2: #basically required just to reject the one night with data from only 2 specs that was in exposures
        tspec = Table.read(coaddir+'/zbest-'+str(specs[0])+'-'+str(tile)+'-'+night+'.fits',hdu='ZBEST')
        tf = Table.read(coaddir+'/coadd-'+str(specs[0])+'-'+str(tile)+'-'+night+'.fits',hdu='FIBERMAP')
        if tsmd == 'fm':
            ts = Table.read(coaddir+'/coadd-'+str(specs[0])+'-'+str(tile)+'-'+night+'.fits',hdu='SCORES')
            ts.keep_columns(['TARGETID','TSNR2_ELG','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
        #this is all to get the effective coadded exposure depth; should eventually just be in the fibermap hdu
        zfm = Table.read(coaddir+'/zbest-'+str(specs[0])+'-'+str(tile)+'-'+night+'.fits',hdu='FIBERMAP')
        try:
            exps = np.unique(zfm['EXPID'])
        except:
            print('NO EXPID COLUMN in '+coaddir+'/zbest-'+str(specs[0])+'-'+str(tile)+'-'+night+'.fits')
            return None    
        etdark = []
        etbright = []
        etback = []
        
        ce = 0
        for exp in exps:
            info = exposures[exposures['EXPID'] == exp]
            if len(info) == 0:
                print('did not find info for expid '+str(exp))
                fo = open(mfn,'a')
                fo.write(str(exp)+'\n')
                return None
            else:    
                
                nt = str(info['NIGHT'][0])
                etdark.append(info['EFFTIME_DARK'][0])
                etbright.append(info['EFFTIME_BRIGHT'][0])
                etback.append(info['EFFTIME_BACKUP'][0]) 
                
                if tsmd != 'fm':
                    for cam in cams:
                        tcols  =[]
                        for col in tsnrcols:
                            tcols.append(col+'_'+cam.upper())
                        cf = Table.read(cfdir+'/'+nt+'/cframe-'+cam+str(specs[0])+'-'+str(exp).zfill(8)+'.fits',hdu='SCORES')
                        cf.keep_columns(tcols)
                        for col in tcols:
                            cf.rename_column(col, col[:-2])
                        if ce ==0 and cam == 'b':
                            tids = Table.read(cfdir+'/'+nt+'/cframe-'+cam+str(specs[0])+'-'+str(exp).zfill(8)+'.fits',hdu='FIBERMAP')
                            tnsrt = cf.copy()
                            tnsrt['TARGETID'] = tids['TARGETID']
                        else:
                            for col in tsnrcols:
                                tnsrt[col] += cf[col]         
                ce += 1                 
        etdarkt = np.zeros(500)
        etbrightt = np.zeros(500)
        etbackt = np.zeros(500)

        tid = zfm[0:500]['TARGETID']
        for i in range(0,len(exps)):
            sel = zfm[i*500:(i+1)*500]
            w = sel['FIBERSTATUS'] == 0
            etdarkt[w] += etdark[i]
            etbrightt[w] += etbright[i]
            etbackt[w] += etback[i]
    
        
        tf['EXPS'] = ",".join(exps.astype(str))
        for i in range(1,len(specs)):
            zfm = Table.read(coaddir+'/zbest-'+str(specs[i])+'-'+str(tile)+'-'+night+'.fits',hdu='FIBERMAP')
            exps = np.unique(zfm['EXPID'])
            if tsmd == 'fm':
                tsn = Table.read(coaddir+'/coadd-'+str(specs[i])+'-'+str(tile)+'-'+night+'.fits',hdu='SCORES')
                tsn.keep_columns(['TARGETID','TSNR2_ELG','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
                ts = vstack([ts,tsn], metadata_conflicts='silent')
           

            tn = Table.read(coaddir+'/zbest-'+str(specs[i])+'-'+str(tile)+'-'+night+'.fits',hdu='ZBEST')
            tnf = Table.read(coaddir+'/coadd-'+str(specs[i])+'-'+str(tile)+'-'+night+'.fits',hdu='FIBERMAP')  
            tnf['EXPS'] = ",".join(exps.astype(str))
            tspec = vstack([tspec,tn], metadata_conflicts='silent')                      
            tf = vstack([tf,tnf], metadata_conflicts='silent')
            etdark = []
            etbright = []
            etback = []
            ce = 0
            for exp in exps:
                info = exposures[exposures['EXPID'] == exp]
                if len(info) == 0:
                    print('did not find info for expid '+str(exp))
                    return None
                else:
                    nt = str(info['NIGHT'][0])
                    etdark.append(info['EFFTIME_DARK'][0])
                    etbright.append(info['EFFTIME_BRIGHT'][0])
                    etback.append(info['EFFTIME_BACKUP'][0]) 
                    
                        
                    ce += 1                 

            
            etdarktn = np.zeros(500)
            etbrighttn = np.zeros(500)
            etbacktn = np.zeros(500)

            tidn = zfm[0:500]['TARGETID']
            for ii in range(0,len(exps)):
                sel = zfm[ii*500:(ii+1)*500]
                w = sel['FIBERSTATUS'] == 0
                etdarktn[w] += etdark[ii]
                etbrighttn[w] += etbright[ii]
                etbacktn[w] += etback[ii]

            etdarkt = np.concatenate([etdarkt,etdarktn])
            etbrightt = np.concatenate([etbrightt,etbrighttn])
            etbackt = np.concatenate([etbackt,etbacktn])   


            tid = np.concatenate([tid,tidn])
        tspec = join(tspec,tf,keys=['TARGETID'], metadata_conflicts='silent')
        if tsmd != 'fm':
            tspec = join(tspec,tnsrt,keys=['TARGETID'], metadata_conflicts='silent')
        else:
            tspec = join(tspec,ts,keys=['TARGETID'], metadata_conflicts='silent')
        td = Table([etdarkt,etbrightt,etbackt,tid],names=('EFFTIME_DARK','EFFTIME_BRIGHT','EFFTIME_BACK','TARGETID'))
        tspec = join(tspec,td,keys=['TARGETID'], metadata_conflicts='silent')
        if tarbit != -1:
            wtype = ((tspec[tp] & 2**tarbit) > 0)
            print(str(len(tspec))+' total entries '+str(len(tspec[wtype]))+' that are requested type entries with '+str(len(np.unique(tspec[wtype]['TARGETID'])))+' unique target IDs')
            if len(tspec[wtype]) == 0:
                return None

            tspec = tspec[wtype]
        tspec['subset'] = night
        tspec['coadd_type'] = ct
        # AR adding a weight for ELGs in the QSO+ELG and QSO+LRG tiles
        # AR to down-weight QSOs which are at a higher priority
        # AR we "rescale" to the ELGxQSO/ELG ratio of the parent input target sample per tile
        if tarbit != -1:
            tspec['elgqso_weight'] = get_elgqso_weight(tarbit,tp,tile,tspec[tp])
        return tspec
    return None    


def get_exp(tarbit,tp,exp,tile,coaddir,exposures,mfn='temp.txt'):

    print('going through exposure '+str(exp))
    specs = []
    #find out which spectrograph have data
    for si in range(0,10):
        try:
            fl = coaddir+'zbest-'+str(si)+'-'+str(tile)+'-000'+str(exp)+'.fits'
            
            fitsio.read(fl)
            specs.append(si)
        except:
            #print(fl,specs,si)
            print(fl)
            print('no spectrograph and/or coadd '+str(si)+ ' on exposure '+str(exp))
    if len(specs) > 2: #basically required just to reject the one night with data from only 2 specs that was in exposures
        tspec = Table.read(coaddir+'zbest-'+str(specs[0])+'-'+str(tile)+'-000'+str(exp)+'.fits',hdu='ZBEST')
        tf = Table.read(coaddir+'zbest-'+str(specs[0])+'-'+str(tile)+'-000'+str(exp)+'.fits',hdu='FIBERMAP')
        tvs = get_tsnrinfo([exp],specs[0])    
        if tvs is not None:
            es,bs,ls,qs = tvs
        else:
            return tvs

        onepet = np.ones(500)
        est = onepet*es
        bst = onepet*bs
        lst = onepet*ls
        qst = onepet*qs
       
        for i in range(1,len(specs)):
            tnf = Table.read(coaddir+'zbest-'+str(specs[i])+'-'+str(tile)+'-000'+str(exp)+'.fits',hdu='FIBERMAP')
            tn = Table.read(coaddir+'zbest-'+str(specs[i])+'-'+str(tile)+'-000'+str(exp)+'.fits',hdu='ZBEST')
            tspec = vstack([tspec,tn], metadata_conflicts='silent')                      
            tf = vstack([tf,tnf], metadata_conflicts='silent')
            
            tvs = get_tsnrinfo([exp],specs[i]) 
            if tvs is not None:
                es,bs,ls,qs = tvs
            else:
                return tvs
            estn = onepet*es
            bstn = onepet*bs
            lstn = onepet*ls
            qstn = onepet*qs
            est = np.concatenate([est,estn])
            lst = np.concatenate([lst,lstn])
            qst = np.concatenate([qst,qstn])   
            bst = np.concatenate([bst,bstn])
        tf['ELGTSNR'] =est
        tf['BGSTSNR'] =bst
        tf['LRGTSNR'] =lst
        tf['QSOTSNR'] = qst  

        tspec = join(tspec,tf,keys=['TARGETID'], metadata_conflicts='silent')
        if tarbit != -1:
            wtype = ((tspec[tp] & 2**tarbit) > 0)
            print(str(len(tspec))+' total entries '+str(len(tspec[wtype]))+' that are requested type entries with '+str(len(np.unique(tspec[wtype]['TARGETID'])))+' unique target IDs')
            tspec = tspec[wtype]
        #tspec['subset'] = night
        # AR adding a weight for ELGs in the QSO+ELG and QSO+LRG tiles
        # AR to down-weight QSOs which are at a higher priority
        # AR we "rescale" to the ELGxQSO/ELG ratio of the parent input target sample per tile
        if tarbit != -1:
            tspec['elgqso_weight'] = get_elgqso_weight(tarbit,tp,tile,tspec[tp])
        return tspec
    return None    



# AR adding a weight for ELGs in the QSO+ELG and QSO+LRG tiles
# AR to down-weight QSOs which are at a higher priority
# AR we "rescale" to the ELGxQSO/ELG ratio of the parent input target sample per tile
# AR returns 1 if not relevant..
# AR weight calculation:
# AR   n_elg : nb of elg targets in the tile (including elgxqso)
# AR   n_elgqso : nb of elgxqso targets in the tile
# AR   fracp : N(elg x qso) / N(elg) on the parent sample, i.e. the input targets to fiberassign
# AR   (n_elgqso * w) / (n_elgqso * w + n_elg-n_elgqso) = fracp
# AR   n_elgqso * w = fracp * (n_elgqso * w + n_elg-n_elgqso)
# AR   (n_elgqso - fracp * n_elgqso) * w = fracp * (n_elg-n_elgqso)
# AR   w = fracp * (n_elg-n_elgqso) / (n_elgqso - fracp * n_elgqso)
def get_elgqso_weight(tarbit,tp,tile,tpval,fadir="{}/survey/fiberassign/SV1/".format(os.getenv("DESI_ROOT"))):
    # AR setting weights=1 by default
    weights = np.ones(len(tpval),dtype=float)
    # AR is the tile a QSO+ELG or QSO+LRG tile?
    fafn = glob("{}/202?????/fiberassign-{}.fits.gz".format(fadir,tile.zfill(6)))[0]
    targfn = glob("{}/202?????/{}-targ.fits".format(fadir,tile.zfill(6)))[0]
    docompute = fits.getheader(fafn)["FAFLAVOR"] in ["cmxelgqso", "cmxlrgqso", "sv1lrgqso", "sv1elgqso"]
    # AR compute only if ELG targets considered and SV1
    docompute &= tarbit == int(np.log2(sv1_targetmask.desi_mask["ELG"]))
    docompute &= tp == "SV1_DESI_TARGET"
    # AR
    if docompute:
        # AR parent sample
        tmpd = fitsio.read(targfn, columns = ["SV1_DESI_TARGET"])
        elg = (tmpd["SV1_DESI_TARGET"] & sv1_targetmask.desi_mask["ELG"]) > 0
        elgqso = (elg) & ((tmpd["SV1_DESI_TARGET"] & sv1_targetmask.desi_mask["QSO"]) > 0)
        fracp = elgqso.sum() / float(elg.sum())
        # AR assigned sample
        tmpd = fitsio.read(fafn, columns = ["SV1_DESI_TARGET"])
        elg = (tmpd["SV1_DESI_TARGET"] & sv1_targetmask.desi_mask["ELG"]) > 0
        elgqso = (elg) & ((tmpd["SV1_DESI_TARGET"] & sv1_targetmask.desi_mask["QSO"]) > 0)
        nelg, nelgqso = elg.sum(), elgqso.sum()
        # AR weight for the ELGxQSO objects in the input sample
        elgqso = ((tpval & sv1_targetmask.desi_mask["ELG"]) > 0) & ((tpval & sv1_targetmask.desi_mask["QSO"]) > 0)
        weights[elgqso] = np.round(fracp * (nelg - nelgqso) / (nelgqso - fracp*nelgqso), 3)
    # AR
    return weights

