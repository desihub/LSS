from astropy.table import Table,join,vstack
import numpy as np
import fitsio

'test'

def comb_subset_vert(tarbit,tp,subsets,tile,coaddir,exposures,outf):
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
        if len(night) > 0:
            tspec = get_subset(tarbit,tp,night,tile,coaddir,exposures)
            if tspec is not None:
                if ss == 0:
                    tspect = tspec
                    ss = 1
                else:
                    tspect = vstack([tspect,tspec])
                print('there are now '+str(len(tspect)) +' entries with '+str(len(np.unique(tspect['TARGETID'])))+' unique target IDs')    
                    
    if ss == 1:
        tspect.sort('TARGETID')
        tspect.write(outf,format='fits', overwrite=True) 
        print('wrote to '+outf)
    else:
        print('no data for tile '+tile)


def get_subset(tarbit,tp,night,tile,coaddir,exposures):

    print('going through subset '+night)
    specs = []
    #find out which spectrograph have data
    for si in range(0,10):
        try:
            fl = coaddir+'/'+night+'/zbest-'+str(si)+'-'+str(tile)+'-'+night+'.fits'
            fitsio.read(fl)
            fl = coaddir+'/'+night+'/coadd-'+str(si)+'-'+str(tile)+'-'+night+'.fits'
            fitsio.read(fl)
            specs.append(si)
        except:
            #print(fl,specs,si)
            print('no spectrograph and/or coadd '+str(si)+ ' on subset '+night)
    if len(specs) > 2: #basically required just to reject the one night with data from only 2 specs that was in exposures
        tspec = Table.read(coaddir+'/'+night+'/zbest-'+str(specs[0])+'-'+str(tile)+'-'+night+'.fits',hdu='ZBEST')
        tf = Table.read(coaddir+'/'+night+'/coadd-'+str(specs[0])+'-'+str(tile)+'-'+night+'.fits',hdu='FIBERMAP')
        #this is all to get the effective coadded exposure depth; should eventually just be in the fibermap hdu
        zfm = Table.read(coaddir+'/'+night+'/zbest-'+str(specs[0])+'-'+str(tile)+'-'+night+'.fits',hdu='FIBERMAP')
        exps = np.unique(zfm['EXPID'])
        bd = []
        rd = []
        zd = []
        for exp in exps:
            info = exposures[exposures['EXPID'] == exp]
            if len(info) == 0:
                print('did not find info for expid '+exp)
                return None
            else:    
                print(info['B_DEPTH'])
                bd.append(info['B_DEPTH'][0])
                rd.append(info['R_DEPTH'][0])
                zd.append(info['Z_DEPTH'][0])        
        bdt = np.zeros(500)
        rdt = np.zeros(500)
        zdt = np.zeros(500)
        tid = zfm[0:500]['TARGETID']
        for i in range(0,len(exps)):
            sel = zfm[i*500:(i+1)*500]
            w = sel['FIBERSTATUS'] == 0
            bdt[w] += bd[i]
            rdt[w] += rd[i]
            zdt[w] += zd[i]
    
        for i in range(1,len(specs)):
            tn = Table.read(coaddir+'/'+night+'/zbest-'+str(specs[i])+'-'+str(tile)+'-'+night+'.fits',hdu='ZBEST')
            tnf = Table.read(coaddir+'/'+night+'/coadd-'+str(specs[i])+'-'+str(tile)+'-'+night+'.fits',hdu='FIBERMAP')  
            tspec = vstack([tspec,tn])                      
            tf = vstack([tf,tnf])
            zfm = Table.read(coaddir+'/'+night+'/zbest-'+str(specs[i])+'-'+str(tile)+'-'+night+'.fits',hdu='FIBERMAP')
            exps = np.unique(zfm['EXPID'])
            bd = []
            rd = []
            zd = []
            for exp in exps:
                info = exposures[exposures['EXPID'] == exp]
                if len(info) == 0:
                    print('did not find info for expid '+str(exp))
                    return None
                else:
                    bd.append(info['B_DEPTH'][0])
                    rd.append(info['R_DEPTH'][0])
                    zd.append(info['Z_DEPTH'][0])        
            bdtn = np.zeros(500)
            rdtn = np.zeros(500)
            zdtn = np.zeros(500)
            tidn = zfm[0:500]['TARGETID']
            for ii in range(0,len(exps)):
                sel = zfm[ii*500:(ii+1)*500]
                w = sel['FIBERSTATUS'] == 0
                bdtn[w] += bd[ii]
                rdtn[w] += rd[ii]
                zdtn[w] += zd[ii]
            bdt = np.concatenate([bdt,bdtn])
            rdt = np.concatenate([rdt,rdtn])
            zdt = np.concatenate([zdt,zdtn])   
            tid = np.concatenate([tid,tidn])
            #print(np.min(rdtn),np.max(rdtn)) 
            #print(np.min(rdt),np.max(rdt)) 
        tspec = join(tspec,tf,keys=['TARGETID'])
        td = Table([bdt,rdt,zdt,tid],names=('B_DEPTH','R_DEPTH','Z_DEPTH','TARGETID'))
        tspec = join(tspec,td,keys=['TARGETID'])
        wtype = ((tspec[tp] & 2**tarbit) > 0)
        print(str(len(tspec))+' total entries '+str(len(tspec[wtype]))+' that are requested type entries with '+str(len(np.unique(tspec[wtype]['TARGETID'])))+' unique target IDs')
        tspec = tspec[wtype]
        tspec['subset'] = night
        return tspec
    return None    
        
        
