#Tools to study and correct for trends in spectroscopic succes rate (ssr)
#Initial LRG model fitting taken from Ronpgpu Zhou's notebook

import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
from scipy.optimize import curve_fit, minimize
from astropy.coordinates import SkyCoord
from astropy import units as u

import LSS.common_tools as common

elgcol = ['SUBSET','EBV','PRIORITY','TARGETID','OII_FLUX','OII_FLUX_IVAR','ELG_LOP','ELG_VLO','TSNR2_ELG','TSNR2_LRG','PHOTSYS','MASKBITS','FIBERFLUX_G','FIBERFLUX_R','FIBERFLUX_Z','COADD_FIBERSTATUS','Z','ZWARN','DELTACHI2']


def ELG_goodobs(data,fbs_col='COADD_FIBERSTATUS'):#,dt_col='DESI_TARGET'):
    mask = data[fbs_col]==0
    print(fbs_col,np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))

    # Remove "no data" fibers
    mask &= data['ZWARN'] & 2**9==0
    print('& No data', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))

    # Apply imaging mask
    #mask &= data['lrg_mask']==0
    #print('& LRG imaging mask', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))

    data['q'] = ELG_goodz(data)#data['ZWARN']==0
    print('failure rate is '+str(np.sum(~data['q'])/len(data)))
    return data

def ELG_goodz(data,zcol='Z'):
    o2c = np.log10(data['OII_FLUX'] * np.sqrt(data['OII_FLUX_IVAR']))+0.2*np.log10(data['DELTACHI2'])
    sel = o2c > 0.9
    return sel


def LRG_goodobs(data,fbs_col='COADD_FIBERSTATUS',dt_col='DESI_TARGET'):
    mask = data[fbs_col]==0
    print(fbs_col,np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))

    # Remove "no data" fibers
    mask &= data['ZWARN'] & 2**9==0
    print('& No data', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))

    # Apply LRG mask
    #mask &= data['lrg_mask']==0
    #print('& LRG imaging mask', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))

    # Remove QSO targets
    mask &= data[dt_col] & 2**2 ==0
    print('& Remove QSO targets', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
    data = data[mask]
    data['q'] = LRG_goodz(data)#data['ZWARN']==0
    #data['q'] &= data['Z']<1.5
    #data['q'] &= data['DELTACHI2']>15  
    print('failure rate is '+str(np.sum(~data['q'])/len(data)))
    return data

def LRG_goodz(data,zcol='Z'):
    sel = data['ZWARN']==0
    sel &= data[zcol]<1.5
    sel &= data['DELTACHI2']>15  
    return sel

def get_ELG_data_full(tracer,surveys=['DA02'],versions=['test'],specrels=['guadalupe']):
    
    cats = []
    for sur,ver,sr in zip(surveys,versions,specrels):
        dir = '/global/cfs/cdirs/desi/survey/catalogs/'+sur+'/LSS/'+sr+'/LSScats/'+ver+'/'
        tfn = tracer
        #if sur == 'DA02':
        #    tfn+='zdone'
        fn = dir+tfn+'_full.dat.fits'    
        data = Table(fitsio.read(fn))
        #print(len(data))
        sel = data['ZWARN'] != 999999
        data = data[sel]
        #print(len(data))
        data['q'] = data['o2c'] > 0.9
        cats.append(data)
        print('# of spectra for fits:'+str(len(data)))
        print('# of good z for fits:'+str(np.sum(data['q'])))

    if len(cats) == 1:
        cat = cats[0]

    cat['EFFTIME_ELG'] = 8.60 * cat['TSNR2_ELG']
    cat['EFFTIME_LRG'] = 12.15 * cat['TSNR2_LRG']
    cat['zfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_Z']) - 1.211 * cat['EBV']
    cat['FIBERFLUX_Z_EC'] = cat['FIBERFLUX_Z']*10**(0.4*1.211*cat['EBV'])
    gextc = 3.214
    cat['gfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_G']) - gextc * cat['EBV']
    cat['FIBERFLUX_G_EC'] = cat['FIBERFLUX_G']*10**(0.4*gextc*cat['EBV'])
    rextc = 2.165
    cat['rfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_R']) - rextc * cat['EBV']
    cat['FIBERFLUX_R_EC'] = cat['FIBERFLUX_R']*10**(0.4*rextc*cat['EBV'])
    cat['qf'] = np.array(cat['q'], dtype=float)
    
    return cat

def get_BGS_data_full(tracer,surveys=['DA02'],versions=['test'],specrels=['guadalupe']):
    
    cats = []
    for sur,ver,sr in zip(surveys,versions,specrels):
        dir = '/global/cfs/cdirs/desi/survey/catalogs/'+sur+'/LSS/'+sr+'/LSScats/'+ver+'/'
        tfn = tracer
        #if sur == 'DA02':
        #    tfn+='zdone'
        fn = dir+tfn+'_full.dat.fits'  
        print('loading info from '+fn)  
        data = Table(fitsio.read(fn))
        #print(len(data))
        sel = data['ZWARN'] != 999999
        data = data[sel]
        #print(len(data))
        gz = data['ZWARN'] == 0
        gz &= data['DELTACHI2'] > 40
        data['q'] = gz
        cats.append(data)

    if len(cats) == 1:
        cat = cats[0]

    cat['EFFTIME_ELG'] = 8.60 * cat['TSNR2_ELG']
    cat['EFFTIME_LRG'] = 12.15 * cat['TSNR2_LRG']
    cat['EFFTIME_BGS'] = 12.15/89.8 * cat['TSNR2_BGS']
    cat['zfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_Z']) - 1.211 * cat['EBV']
    cat['FIBERFLUX_Z_EC'] = cat['FIBERFLUX_Z']*10**(0.4*1.211*cat['EBV'])
    gextc = 3.214
    cat['gfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_G']) - gextc * cat['EBV']
    cat['FIBERFLUX_G_EC'] = cat['FIBERFLUX_G']*10**(0.4*gextc*cat['EBV'])
    rextc = 2.165
    cat['rfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_R']) - rextc * cat['EBV']
    cat['FIBERFLUX_R_EC'] = cat['FIBERFLUX_R']*10**(0.4*rextc*cat['EBV'])
    cat['qf'] = np.array(cat['q'], dtype=float)
    
    return cat


def get_QSO_data_full(tracer,surveys=['DA02'],versions=['test'],specrels=['guadalupe'],cut=None,cut_condition=None):
    
    cats = []
    for sur,ver,sr in zip(surveys,versions,specrels):
        dir = '/global/cfs/cdirs/desi/survey/catalogs/'+sur+'/LSS/'+sr+'/LSScats/'+ver+'/'
        tfn = tracer
        #if sur == 'DA02':
        #    tfn+='zdone'
        fn = dir+tfn+'_full.dat.fits'    
        data = Table(fitsio.read(fn))
        print(len(data))
        sel = data['ZWARN'] != 999999
        sel &= ((data['SPECTYPE'] == 'STAR') & (data['Z_not4clus'] != 999999)) | (data['SPECTYPE'] != 'STAR')
        data = data[sel]
        
        if cut == 'galactic':
            coords = SkyCoord(ra=data['RA']*u.deg,dec=data['DEC']*u.deg)
            if cut_condition[0] == '>':
                cond = np.where(np.abs(coords.galactic.b.value) > float(cut_condition[1:]))
            elif cut_condition[0] == '<':
                cond = np.where(np.abs(coords.galactic.b.value) < float(cut_condition[1:]))
            data = data[cond]

        wz = data['Z_not4clus']*0 == 0
        wz &= data['Z_not4clus'] != 999999
        wz &= data['Z_not4clus'] != 1.e20

        print(len(data),len(wz),np.sum(wz))
        data['q'] = wz
        cats.append(data)

    if len(cats) == 1:
        cat = cats[0]

    cat['EFFTIME_ELG'] = 8.60 * cat['TSNR2_ELG']
    cat['EFFTIME_QSO'] = 8.60/0.255 * cat['TSNR2_QSO']
    cat['EFFTIME_LRG'] = 12.15 * cat['TSNR2_LRG']
    cat['zfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_Z']) - 1.211 * cat['EBV']
    cat['FIBERFLUX_Z_EC'] = cat['FIBERFLUX_Z']*10**(0.4*1.211*cat['EBV'])
    gextc = 3.214
    cat['gfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_G']) - gextc * cat['EBV']
    cat['FIBERFLUX_G_EC'] = cat['FIBERFLUX_G']*10**(0.4*gextc*cat['EBV'])
    rextc = 2.165
    cat['rfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_R']) - rextc * cat['EBV']
    cat['FIBERFLUX_R_EC'] = cat['FIBERFLUX_R']*10**(0.4*rextc*cat['EBV'])
    cat['qf'] = np.array(cat['q'], dtype=float)
    
    return cat

def get_data_full(tracer,surveys=['DA02'],versions=['test'],specrels=['guadalupe'],cut=None,cut_condition=None):
    
    cats = []
    for sur,ver,sr in zip(surveys,versions,specrels):
        dir = '/global/cfs/cdirs/desi/survey/catalogs/'+sur+'/LSS/'+sr+'/LSScats/'+ver+'/'
        tfn = tracer
        #if sur == 'DA02':
        #    tfn+='zdone'
        fn = dir+tfn+'_full.dat.fits'  
        print('loading info from '+fn)  
        data = Table(fitsio.read(fn))
        #print(len(data))
        sel = data['ZWARN'] != 999999
        if tracer[:3] == 'QSO':
            sel &= ((data['SPECTYPE'] == 'STAR') & (data['Z_not4clus'] != 999999)) | (data['SPECTYPE'] != 'STAR')
        data = data[sel]

        if cut == 'galactic':
            coords = SkyCoord(ra=data['RA']*u.deg,dec=data['DEC']*u.deg)
            if cut_condition[0] == '>':
                cond = np.where(np.abs(coords.galactic.b.value) > float(cut_condition[1:]))
            elif cut_condition[0] == '<':
                cond = np.where(np.abs(coords.galactic.b.value) < float(cut_condition[1:]))
            data = data[cond]

        #print(len(data))

        if tracer[:3] == 'BGS':
            gz = data['ZWARN'] == 0
            gz &= data['DELTACHI2'] > 40
        
        if tracer[:3] == 'LRG':
            gz = data['ZWARN']==0
            gz &= data['Z_not4clus']<1.5
            gz &= data['DELTACHI2']>15  


        if tracer[:3] == 'QSO':
            gz = data['Z_not4clus']*0 == 0
            gz &= data['Z_not4clus'] != 999999
            gz &= data['Z_not4clus'] != 1.e20
        
        if tracer[:3] == 'ELG':
            gz = data['o2c'] > 0.9

        data['q'] = gz
        print('# of spectra for fits:'+str(len(data)))
        print('# of good z for fits:'+str(np.sum(data['q'])))
        
        cats.append(data)

    if len(cats) == 1:
        cat = cats[0]

    cat['EFFTIME_ELG'] = 8.60 * cat['TSNR2_ELG']
    cat['EFFTIME_LRG'] = 12.15 * cat['TSNR2_LRG']
    cat['EFFTIME_QSO'] = 8.60/0.255 * cat['TSNR2_QSO']
    cat['EFFTIME_BGS'] = 12.15/89.8 * cat['TSNR2_BGS']
    cat['zfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_Z']) - 1.211 * cat['EBV']
    cat['FIBERFLUX_Z_EC'] = cat['FIBERFLUX_Z']*10**(0.4*1.211*cat['EBV'])
    gextc = 3.214
    cat['gfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_G']) - gextc * cat['EBV']
    cat['FIBERFLUX_G_EC'] = cat['FIBERFLUX_G']*10**(0.4*gextc*cat['EBV'])


    rextc = 2.165
    cat['rfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_R']) - rextc * cat['EBV']
    cat['FIBERFLUX_R_EC'] = cat['FIBERFLUX_R']*10**(0.4*rextc*cat['EBV'])
    cat['qf'] = np.array(cat['q'], dtype=float)
    
    return cat



def get_ELG_data(specrel='fuji',tr='ELG_LOP',maskbits=[1,11,12,13],notqso=True):
    maintids = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/'+tr+'targetsDR9v1.1.1.fits',columns=['TARGETID','DESI_TARGET','MASKBITS','NOBS_G','NOBS_R','NOBS_Z'])
    maintids = common.cutphotmask(maintids,maskbits)
    elgcatdir = '/global/cfs/cdirs/desi/users/raichoor/spectro/'+specrel
    

    sv3 = fitsio.read(elgcatdir+'/sv3-elg-fuji-tiles.fits',columns=elgcol)
    st = []
    for i in range(0,len(sv3)):
        st.append(sv3['SUBSET'][i][:4])
    st = np.array(st)
    wg = st == "thru"
    sv3 = sv3[wg]

    if tr != 'ELG':
        print('cutting SV3 to main '+tr)
        sel = sv3[tr] == True
        print('length before is '+str(len(sv3)))
        sv3 = sv3[sel]
        print('length after is '+str(len(sv3)))
    sel = sv3['PRIORITY'] > 10000
    sv3 = sv3[sel]
    print('length after cutting to priority > 10000 '+str(len(sv3)))
    sv3 = ELG_goodobs(Table(sv3))
    
    sv3 = join(sv3,maintids,keys=['TARGETID'])
    print('length after join to main targets to get DESI_TARGET and cut on maskbits values '+str(len(sv3)))
    
    
    elgcatdirg = '/global/cfs/cdirs/desi/users/raichoor/spectro/guadalupe'
    
    main = fitsio.read(elgcatdirg+'/main-elg-guadalupe-tiles.fits',columns=elgcol)
    st = []
    for i in range(0,len(main)):
        st.append(main['SUBSET'][i][:4])
    st = np.array(st)
    wg = st == "thru"
    main = main[wg]

    if tr != 'ELG':
        print('cutting main to main '+tr)
        sel = main[tr] == True
        print('length before is '+str(len(main)))
        main = main[sel]
        print('length after is '+str(len(main)))
    main = ELG_goodobs(Table(main))
    main = join(main,maintids,keys=['TARGETID'])
    print('length after join to main targets to get DESI_TARGET and cut on maskbits values '+str(len(main)))
    

    sv1 = fitsio.read(elgcatdir+'/sv1-elg-fuji-tiles.fits',columns=elgcol)
    if tr != 'ELG':
        print('cutting SV1 to main '+tr)
        sel = sv1[tr] == True
        print('length before is '+str(len(sv1)))
        sv1 = sv1[sel]
        print('length after is '+str(len(sv1)))
    sv1 = ELG_goodobs(Table(sv1))
    sv1 = join(sv1,maintids,keys=['TARGETID'])
    print('length after join to main targets to get DESI_TARGET and cut on maskbits values '+str(len(sv1)))
    


    #cat = vstack([sv1, sv3, main], join_type='inner')
    #cat = vstack([sv1, main], join_type='inner')
    cat = main
    print(len(cat))

    if notqso:
        # Remove QSO targets
        mask = cat['DESI_TARGET'] & 2**2 ==0
        print(' Remove QSO targets', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
        cat = cat[mask]


    cat['EFFTIME_ELG'] = 8.60 * cat['TSNR2_ELG']
    cat['EFFTIME_LRG'] = 12.15 * cat['TSNR2_LRG']
    cat['zfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_Z']) - 1.211 * cat['EBV']
    cat['FIBERFLUX_Z_EC'] = cat['FIBERFLUX_Z']*10**(0.4*1.211*cat['EBV'])
    gextc = 3.214
    cat['gfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_G']) - gextc * cat['EBV']
    cat['FIBERFLUX_G_EC'] = cat['FIBERFLUX_G']*10**(0.4*gextc*cat['EBV'])
    rextc = 2.165
    cat['rfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_R']) - rextc * cat['EBV']
    cat['FIBERFLUX_R_EC'] = cat['FIBERFLUX_R']*10**(0.4*rextc*cat['EBV'])
    cat['qf'] = np.array(cat['q'], dtype=float)
    
    return cat
    

def get_LRG_data(specrel='fuji'):
    
    maintids = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/LRGtargetsDR9v1.1.1.fits',columns=['TARGETID','lrg_mask'])
    sel = maintids['lrg_mask'] == 0
    maintids = maintids[sel]
    
    zcatdir = '/global/cfs/cdirs/desi/spectro/redux/'+specrel+'/zcatalog/'

    perexpall = Table(fitsio.read(zcatdir+'ztile-sv1-dark-perexp.fits'))
    sel = np.isin(perexpall['TARGETID'],maintids['TARGETID'])
    perexplrg = perexpall[sel]
    del perexpall
    perexplrg = LRG_goodobs(perexplrg,'COADD_FIBERSTATUS','SV1_DESI_TARGET')
    
    cat_1xall = Table(fitsio.read(zcatdir+'ztile-sv1-dark-1x_depth.fits'))
    sel = np.isin(cat_1xall['TARGETID'],maintids['TARGETID'])
    cat_1xlrg = cat_1xall[sel]
    del cat_1xall
    cat_1xlrg = LRG_goodobs(cat_1xlrg,'COADD_FIBERSTATUS','SV1_DESI_TARGET')   
    
    cat_deepall = Table(fitsio.read(zcatdir+'ztile-sv1-dark-cumulative.fits'))
    sel = np.isin(cat_deepall['TARGETID'],maintids['TARGETID'])
    cat_deeplrg = cat_deepall[sel]
    del cat_deepall
    cat_deeplrg = LRG_goodobs(cat_deeplrg,'COADD_FIBERSTATUS','SV1_DESI_TARGET') 
 
    cat_sv3all = Table(fitsio.read(zcatdir+'ztile-sv3-dark-cumulative.fits'))
    sel = np.isin(cat_sv3all['TARGETID'],maintids['TARGETID'])
    sel &= cat_sv3all['PRIORITY'] == 103200 #we don't want to include the failed repeats in the statistics
    cat_sv3lrg = cat_sv3all[sel]
    del cat_sv3all
    cat_sv3lrg = LRG_goodobs(cat_sv3lrg,'COADD_FIBERSTATUS','SV3_DESI_TARGET') 
     
    if specrel == 'fuji':
        specrelmain = 'guadalupe'
        zcatdirm = '/global/cfs/cdirs/desi/spectro/redux/'+specrelmain+'/zcatalog/'

    cat_mainall = Table(fitsio.read(zcatdirm+'ztile-main-dark-cumulative.fits'))
    sel = np.isin(cat_mainall['TARGETID'],maintids['TARGETID'])
    cat_mainlrg = cat_mainall[sel]
    del cat_mainall
    cat_mainlrg = LRG_goodobs(cat_mainlrg,'COADD_FIBERSTATUS','DESI_TARGET') 

    cat = vstack([perexplrg, cat_1xlrg, cat_mainlrg, cat_deeplrg, cat_sv3lrg], join_type='inner')
    print(len(cat))

    cat['EFFTIME_ELG'] = 8.60 * cat['TSNR2_ELG']
    cat['EFFTIME_LRG'] = 12.15 * cat['TSNR2_LRG']
    cat['zfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_Z']) - 1.211 * cat['EBV']
    cat['FIBERFLUX_Z_EC'] = cat['FIBERFLUX_Z']*10**(0.4*1.211*cat['EBV'])
    cat['qf'] = np.array(cat['q'], dtype=float)
    
    return cat

def fit_cons(dl,el,minv=0,step=0.01):
    c = minv
    newcost = np.sum((dl-c)**2./el**2.)
    oldcost = newcost + 1
    while newcost < oldcost:
        oc = c
        oldcost = newcost
        c += step
        newcost = np.sum((dl-c)**2./el**2.)
    
    return oldcost,c


class LRG_ssr:
    def __init__(self,specrel='fuji',efftime_min=500,efftime_max=2000,surveys=['DA02'],versions=['test'],specrels=['guadalupe']):
        #self.cat = get_LRG_data(specrel)
        self.cat = get_data_full('LRG',surveys=surveys,versions=versions,specrels=specrels)
        mask = self.cat['EFFTIME_LRG']>efftime_min
        mask &= self.cat['EFFTIME_LRG']<efftime_max
        self.cat = self.cat[mask]
        print('success rate is '+str(sum(self.cat['qf'])/len(self.cat)))

    def cost(self,q_predict):
        return np.sum((self.cat['qf']-q_predict)**2)

    def wrapper(self,params):
        q_predict = 1-self.failure_rate(self.cat['FIBERFLUX_Z_EC'], self.cat['EFFTIME_LRG'], *params)
        return self.cost(q_predict)

    def failure_rate(self,flux, efftime, a, b, c):
        sn = flux * np.sqrt(efftime)
        #return np.clip(np.exp(-(sn+a)/b)+c/flux, 0, 1)
        return np.clip(np.exp(a*sn+b)+c/flux, 0, 1)

    def add_modpre(self,data):
        res = minimize(self.wrapper, [-.1, 3., 0.02], bounds=((-200, 0), (0, 100), (0., 1)),
               method='Powell', tol=1e-6)
        pars = res.x
        print(pars)
        dflux = data['FIBERFLUX_Z']*10**(0.4*1.211*data['EBV'])#data['FIBERFLUX_Z_EC']
        deff = 12.15 * data['TSNR2_LRG']#data['EFFTIME_LRG']
        data['mod_success_rate'] = 1. -self.failure_rate(dflux,deff,*pars) 
        print(min(data['mod_success_rate']),max(data['mod_success_rate']))
        data['WEIGHT_ZFAIL'] = 1./data['mod_success_rate'] 
        return data

class BGS_ssr:
    def __init__(self,specrel='fuji',efftime_min=120,efftime_max=300,surveys=['DA02'],versions=['test'],specrels=['guadalupe']):
        
        self.cat = get_data_full('BGS_BRIGHT',surveys=surveys,versions=versions,specrels=specrels)
        mask = self.cat['EFFTIME_BGS']>efftime_min
        mask &= self.cat['EFFTIME_BGS']<efftime_max
        print('using '+str(len(self.cat)/len(self.cat[mask])))
        self.cat = self.cat[mask]
        self.selgz = self.cat['q'] == 1
        ha,bine = np.histogram(self.cat['EFFTIME_BGS'])
        hf,_ = np.histogram(self.cat['EFFTIME_BGS'][~self.selgz])
        self.nzf = hf/ha
        print(self.nzf)
        self.nzfe = np.sqrt(hf)/ha
        bc = []
        bs = bine[1]-bine[0]
        for i in range(0,len(bine)-1):
            bc.append(bine[i]+bs/2.) 
        self.bc = np.array(bc)
        self.bine = bine
        self.vis_5hist = False


    def cost(self,q_predict):
        return np.sum((self.cat['qf']-q_predict)**2)

    def wrapper(self,params):
        q_predict = 1-self.failure_rate(self.cat['FIBERFLUX_R_EC'], self.cat['EFFTIME_BGS'], *params)
        return self.cost(q_predict)

    def wrapper_hist(self,params):
        h_predict = self.failure_rate_eff(self.bc, *params)
        diff = self.nzf-h_predict
        cost = np.sum((diff/self.nzfe)**2.)
        return cost


    def failure_rate(self,flux, efftime, a, b, c):
        sn = flux * np.sqrt(efftime)
        return np.clip(np.exp(-(sn+a)/b)+c/flux, 0, .5)

    def failure_rate_eff(self, efftime, a, b, c):
        #sn = flux * np.sqrt(efftime)
        #return np.clip(np.exp(-(sn+a)/b)+c/flux, 0, 1)
        return np.clip(np.exp(-(efftime+a)/b)+c, 0, 1)

    def hist_norm(self,fluxc):
        nzfper = []
        consl = []
        
        nb = 5
        pstep = 100//5
        costt = 0
        for i in range(0,nb):
            sel = self.cat['FIBERFLUX_R_EC'] > np.percentile(self.cat['FIBERFLUX_R_EC'],i*pstep)
            sel &= self.cat['FIBERFLUX_R_EC'] < np.percentile(self.cat['FIBERFLUX_R_EC'],(i+1)*pstep)
            mf = np.median(self.cat['FIBERFLUX_R_EC'][sel])
            if self.vis_5hist:
                print(mf)
            #fper.append(mf)
            wtf = (fluxc*(self.mft-self.cat['FIBERFLUX_R_EC'])/self.mft+1)*(self.wts_fid-1)+1
            selw = wtf < 1
            wtf[selw] = 1
            ha,_ = np.histogram(self.cat['EFFTIME_BGS'][sel],bins=self.bine)
            hf,_ = np.histogram(self.cat['EFFTIME_BGS'][sel&self.selgz],weights=wtf[sel&self.selgz],bins=self.bine)
            #if self.vis_5hist:
            #    print(mf)
            #    print(np.sum(ha))
            #    print(np.sum(hf))
            dl = hf/ha
            nzfper.append(dl)
            def ccost(c):
                return np.sum((dl-c)**2./self.nzfpere[i]**2.)
            resc = minimize(ccost, np.ones(1))
            bc = resc.x
            cost = ccost(bc)
            consl.append(bc)
            costt += cost
        if self.vis_5hist:
            for i in range(0,nb):
                plt.errorbar(self.bc,nzfper[i],self.nzfpere[i])
                plt.plot(self.bc,np.ones(len(self.bc))*consl[i],'k:')
            plt.show()
            self.consl = consl
        return costt    



    def add_modpre(self,data):#,fn_root):
        #res = minimize(self.wrapper, [0, 10., 0.01], bounds=((-200, 200), (0, 100), (0., 1)),
        #       method='Powell', tol=1e-6)
        #pars = res.x
        #print(pars,self.wrapper(pars))
        res = minimize(self.wrapper_hist, [-2, 3, 0.01], bounds=((-1000, 0), (0, 1000), (0., 0.02)),
               method='Powell', tol=1e-6)
        pars = res.x
        chi2 = self.wrapper_hist(pars)
        print(pars,chi2)
        plt.errorbar(self.bc,self.nzf,self.nzfe,fmt='ko',label='data')
        mod = self.failure_rate_eff(self.bc, *pars)
        plt.plot(self.bc,mod,'k--',label='model; chi2='+str(round(chi2,3)))
        plt.ylabel('BGS_BRIGHT Z failure rate')
        plt.xlabel('BGS_BRIGHT EFFECTIVE exp time')
        plt.legend()
        #plt.savefig(fn_root+'overall_failratefit.png')        
        plt.show()
        plt.clf()

        dflux = data['FIBERFLUX_R']*10**(0.4*2.165*data['EBV'])#data['FIBERFLUX_Z_EC']
        deff = 12.15/89.8 * data['TSNR2_BGS']#data['EFFTIME_LRG']
        
        #data['mod_success_rate'] = 1. -self.failure_rate_eff(deff,*pars)   
        #assr = 1. -self.failure_rate_eff(self.cat['EFFTIME_BGS'],*pars)   
        #relssr = assr/np.max(assr) 
        #drelssr = data['mod_success_rate']/np.max(assr)#np.max(data['mod_success_rate'])
        #self.wts_fid = 1/relssr
        fail_frac_mod_data = self.failure_rate_eff(deff,*pars) #
        fail_frac_mod_cat = self.failure_rate_eff(self.cat['EFFTIME_BGS'],*pars)  #
        minfail = np.min(fail_frac_mod_cat )
        relfail = fail_frac_mod_cat/minfail
        relfail_data = fail_frac_mod_data/minfail
        self.wts_fid = (1.-minfail)/(1.-fail_frac_mod_cat)#relfail
        nzfper = []
        nzfpere = []
        fper = []
        self.mft = np.median(self.cat['FIBERFLUX_R_EC'])
        nb = 5
        pstep = 100//5
        for i in range(0,nb):
            sel = self.cat['FIBERFLUX_R_EC'] > np.percentile(self.cat['FIBERFLUX_R_EC'],i*pstep)
            sel &= self.cat['FIBERFLUX_R_EC'] < np.percentile(self.cat['FIBERFLUX_R_EC'],(i+1)*pstep)
            mf = np.median(self.cat['FIBERFLUX_R_EC'][sel])
            fper.append(mf)
            ha,_ = np.histogram(self.cat['EFFTIME_BGS'][sel],bins=self.bine)
            hf,_ = np.histogram(self.cat['EFFTIME_BGS'][sel&self.selgz],bins=self.bine)
            hfw,_ = np.histogram(self.cat['EFFTIME_BGS'][sel&self.selgz],weights=self.wts_fid[sel&self.selgz],bins=self.bine)
            nzfper.append(hf/ha)
            nzfpere.append(np.sqrt(ha-hf)/ha)
            #plt.plot(self.bc,hfw/ha)
        #plt.title('inputs')
        #plt.show()
        self.nzfpere = nzfpere    
        print(nzfpere)
        rest = minimize(self.hist_norm, np.ones(1))#, bounds=((-10, 10)),
               #method='Powell', tol=1e-6)
        fcoeff = rest.x
        self.vis_5hist = True
        print(fcoeff,self.hist_norm(fcoeff),self.consl,fper)#,self.hist_norm(0.),self.hist_norm(1.)) 
        plt.plot(fper,self.consl)
        plt.show()
        #wtf = (fcoeff*(self.mft-dflux)/self.mft+1)*(1/drelssr-1)+1
        wtf = (fcoeff*(self.mft-dflux)/self.mft+1)*((1.-minfail)/(1.-fail_frac_mod_data)-1)+1
        data['mod_success_rate'] = 1. - fail_frac_mod_data +minfail
        minfail = np.zeros(len(dflux))
        for i in range(0,nb):
            sel = dflux > np.percentile(self.cat['FIBERFLUX_R_EC'],i*pstep)
            sel &= dflux < np.percentile(self.cat['FIBERFLUX_R_EC'],(i+1)*pstep)
            minfail_bin = 1-self.consl[i][0]
            minfail[sel] = minfail_bin
            print(np.unique(minfail[sel]),1-minfail_bin)
        data['mod_success_rate'] = data['mod_success_rate']-minfail
        print(np.min(data['mod_success_rate']))
        #minfail_flux = (fcoeff*(self.mft-dflux)/self.mft+1)*
        sel = wtf < 1
        wtf[sel] = 1
        data['WEIGHT_ZFAIL'] =  wtf
        return data


        #print(len(data),np.sum(data['mod_success_rate']))
#         ha,_ = np.histogram(deff,bins=self.bine)
#         gz = data['ZWARN'] == 0
#         gz &= data['DELTACHI2'] > 40
#         hf,_ = np.histogram(deff[gz],weights=1/data[gz]['mod_success_rate'],bins=self.bine)
#         plt.errorbar(self.bc,1.-self.nzf,self.nzfe,fmt='ko')
#         plt.errorbar(self.bc,hf/ha,self.nzfe,fmt='rd')
#         
#         plt.show()
# 
# 
#         return data


class ELG_ssr:
    def __init__(self,specrel='fuji',efftime_min=450,efftime_max=1500,surveys=['DA02'],versions=['test'],specrels=['guadalupe'],reg=None):
        #self.cat = get_ELG_data_full('ELG_LOPnotqso')#,surveys=surveys,versions=versions,specrels=specrels)#get_ELG_data(specrel)
        self.cat = get_data_full('ELG_LOPnotqso',surveys=surveys,versions=versions,specrels=specrels)#get_ELG_data(specrel)

        mask = self.cat['EFFTIME_ELG']>efftime_min
        mask &= self.cat['EFFTIME_ELG']<efftime_max
        self.reg = reg
        if reg is not None:
            mask &= self.cat['PHOTSYS'] == reg
        self.cat = self.cat[mask]
        self.selgz = self.cat['q'] == 1
        ha,bine = np.histogram(self.cat['EFFTIME_ELG'])
        hf,_ = np.histogram(self.cat['EFFTIME_ELG'][~self.selgz])
        self.nzf = hf/ha
        print(self.nzf)
        self.nzfe = np.sqrt(hf)/ha
        bc = []
        bs = bine[1]-bine[0]
        for i in range(0,len(bine)-1):
            bc.append(bine[i]+bs/2.) 
        self.bc = np.array(bc)
        self.bine = bine
        self.vis_5hist = False
        self.outdir = '/global/cfs/cdirs/desi/survey/catalogs/'+surveys[0]+'/LSS/'+specrels[0]+'/LSScats/'+versions[0]+'/'
        
        
        
    def cost(self,q_predict):
        return np.sum((self.cat['qf']-q_predict)**2)

    def wrapper(self,params):
        q_predict = 1-self.failure_rate(self.cat['FIBERFLUX_G_EC'], self.cat['EFFTIME_ELG'], *params)
        return self.cost(q_predict)

    def wrapper_hist(self,params):
        h_predict = self.failure_rate_eff(self.bc, *params)
        diff = self.nzf-h_predict
        cost = np.sum((diff/self.nzfe)**2.)
        return cost

    def failure_rate(self,flux, efftime, a, b, c):
        #sn = flux * np.sqrt(efftime)
        #return np.clip(np.exp(-(sn+a)/b)+c/flux, 0, 1)
        return np.clip(np.exp(-(efftime+a)/b)+c/flux, 0, 1)

    def failure_rate_eff(self, efftime, a, b, c):
        #sn = flux * np.sqrt(efftime)
        #return np.clip(np.exp(-(sn+a)/b)+c/flux, 0, 1)
        return np.clip(np.exp(-(efftime+a)/b)+c, 0, 1)

    
    def hist_norm(self,fluxc):
        nzfper = []
        consl = []
        
        nb = 5
        pstep = 100//5
        costt = 0
        for i in range(0,nb):
            sel = self.cat['FIBERFLUX_G_EC'] > np.percentile(self.cat['FIBERFLUX_G_EC'],i*pstep)
            sel &= self.cat['FIBERFLUX_G_EC'] < np.percentile(self.cat['FIBERFLUX_G_EC'],(i+1)*pstep)
            mf = np.median(self.cat['FIBERFLUX_G_EC'][sel])
            if self.vis_5hist:
                print(mf)
            #fper.append(mf)
            wtf = (fluxc*(self.mft-self.cat['FIBERFLUX_G_EC'])/self.mft+1)*(self.wts_fid-1)+1
            selw = wtf < 1
            wtf[selw] = 1
            ha,_ = np.histogram(self.cat['EFFTIME_ELG'][sel],bins=self.bine)
            hf,_ = np.histogram(self.cat['EFFTIME_ELG'][sel&self.selgz],weights=wtf[sel&self.selgz],bins=self.bine)
            #if self.vis_5hist:
            #    print(mf)
            #    print(np.sum(ha))
            #    print(np.sum(hf))
            dl = hf/ha
            nzfper.append(dl)
            def ccost(c):
                return np.sum((dl-c)**2./self.nzfpere[i]**2.)
            resc = minimize(ccost, np.ones(1))
            bc = resc.x
            cost = ccost(bc)
            consl.append(bc)
            costt += cost
        if self.vis_5hist:
            for i in range(0,nb):
                plt.errorbar(self.bc,nzfper[i],self.nzfpere[i])
                plt.plot(self.bc,np.ones(len(self.bc))*consl[i],'k:')
            plt.ylabel('ELG_LOPnotqso Z failure rate, in fiber bins')
            plt.xlabel('ELG EFFECTIVE exp time')
            plt.legend()
            plt.savefig(self.outdir+'ELG_LOPnotqso_gfibbin_failratefit.png')        

            plt.show()
        return costt    
        
    
    def add_modpre(self,data):
        res = minimize(self.wrapper_hist, [-20, 225., 0.28], bounds=((-1000, 0), (0, 1000), (0., 1)),
               method='Powell', tol=1e-6)
        pars = res.x
        chi2 = self.wrapper_hist(pars)
        print(pars,chi2)
        rw = ''
        if self.reg is not None:
            rw = self.reg
        fo = open(self.outdir+'ELG_LOPnotqso_'+rw+'pars.txt','w')
        fo.write('#overall fit\n')
        fo.write('#a b c chi2\n')
        for par in pars:
            fo.write(str(par)+' ')
        fo.write(str(chi2)+'\n')
        plt.errorbar(self.bc,self.nzf,self.nzfe,fmt='ko',label='data')
        mod = self.failure_rate_eff(self.bc, *pars)
        plt.plot(self.bc,mod,'k--',label='model; chi2='+str(round(chi2,3)))
        plt.ylabel('ELG_LOPnotqso Z failure rate')
        plt.xlabel('ELG EFFECTIVE exp time')
        plt.legend()
        plt.savefig(self.outdir+'ELG_LOPnotqso_'+rw+'overall_failratefit.png')        
        plt.show()
        plt.clf()

        gextc = 3.214
        dflux = data['FIBERFLUX_G']*10**(0.4*gextc*data['EBV']) #data['FIBERFLUX_G_EC']
        deff = 8.60 * data['TSNR2_ELG']#data['EFFTIME_ELG']
        #data['mod_success_rate'] = 1. -self.failure_rate(dflux,deff,*pars) 
        data['mod_success_rate'] = 1. -self.failure_rate_eff(deff,*pars)   
        assr = 1. -self.failure_rate_eff(self.cat['EFFTIME_ELG'],*pars)   
        relssr = assr/np.max(assr) 
        drelssr = data['mod_success_rate']/np.max(assr)#np.max(data['mod_success_rate'])
        seld = deff > 450
        seld &= deff < 1500
        print(len(relssr),len(drelssr[seld]),np.max(assr),np.max(data[seld]['mod_success_rate']))
        self.wts_fid = 1/relssr
        nzfper = []
        nzfpere = []
        fper = []
        self.mft = np.median(self.cat['FIBERFLUX_G_EC'])
        nb = 5
        pstep = 100//5
        for i in range(0,nb):
            sel = self.cat['FIBERFLUX_G_EC'] > np.percentile(self.cat['FIBERFLUX_G_EC'],i*pstep)
            sel &= self.cat['FIBERFLUX_G_EC'] < np.percentile(self.cat['FIBERFLUX_G_EC'],(i+1)*pstep)
            mf = np.median(self.cat['FIBERFLUX_G_EC'][sel])
            fper.append(mf)
            ha,_ = np.histogram(self.cat['EFFTIME_ELG'][sel],bins=self.bine)
            hf,_ = np.histogram(self.cat['EFFTIME_ELG'][sel&self.selgz],bins=self.bine)
            hfw,_ = np.histogram(self.cat['EFFTIME_ELG'][sel&self.selgz],weights=self.wts_fid[sel&self.selgz],bins=self.bine)
            nzfper.append(hf/ha)
            nzfpere.append(np.sqrt(ha-hf)/ha)
            #plt.plot(self.bc,hfw/ha)
        #plt.title('inputs')
        #plt.show()
        self.nzfpere = nzfpere    
        rest = minimize(self.hist_norm, np.ones(1))#, bounds=((-10, 10)),
               #method='Powell', tol=1e-6)
        fcoeff = rest.x
        self.vis_5hist = True
        chi2 = self.hist_norm(fcoeff)
        print(fcoeff,chi2)#,self.hist_norm(0.),self.hist_norm(1.)) 
        fo.write('#gflux fit\n')
        fo.write('#fcoeff chi2\n')
        
        fo.write(str(fcoeff)+' ')
        fo.write(str(chi2)+'\n')
        fo.close()
        wtf = (fcoeff*(self.mft-dflux)/self.mft+1)*(1/drelssr-1)+1
        sel = wtf < 1
        wtf[sel] = 1
        data['WEIGHT_ZFAIL'] =  wtf
        return data

#         nb = 5
#         pstep = 100//5
#         costt = 0
#         
#         seld = np.ones(len(dflux),dtype='bool')
#         dflux = dflux[seld]
#         deff =deff[seld]
#         dselgz = data[seld]['o2c'] > 0.9
#         wtf = (1/drelssr[seld]-1)+1
        #print('are weight arrays equal?',np.array_equal(self.wts_fid,wtf))
#         for i in range(0,nb):
#             sel = dflux > np.percentile(dflux,i*pstep)
#             sel &= dflux < np.percentile(dflux,(i+1)*pstep)
#             mf = np.median(dflux[sel])
#             
#             
#             
#             ha,_ = np.histogram(deff[sel],bins=self.bine)
#             hf,_ = np.histogram(deff[sel&dselgz],weights=wtf[sel&dselgz],bins=self.bine)

class QSO_ssr:
    def __init__(self,specrel='fuji',efftime_min=450,efftime_max=1500,cut=None,cut_condition=None,surveys=['DA02'],versions=['test'],specrels=['guadalupe']):
        self.cat = get_data_full('QSO',cut=cut,cut_condition=cut_condition,surveys=surveys,versions=versions,specrels=specrels)#get_ELG_data(specrel)
        mask = self.cat['EFFTIME_QSO']>efftime_min
        mask &= self.cat['EFFTIME_QSO']<efftime_max
        self.cat = self.cat[mask]
        self.selgz = self.cat['q'] == 1
        ha,bine = np.histogram(self.cat['EFFTIME_QSO'])
        hf,_ = np.histogram(self.cat['EFFTIME_QSO'][~self.selgz])
        self.nzf = hf/ha
        print(self.nzf)
        self.nzfe = np.sqrt(hf)/ha
        bc = []
        bs = bine[1]-bine[0]
        for i in range(0,len(bine)-1):
            bc.append(bine[i]+bs/2.) 
        self.bc = np.array(bc)
        self.bine = bine
        self.vis_5hist = False
        self.cut = cut
        self.cut_condition = cut_condition
        
        
        
    def cost(self,q_predict):
        return np.sum((self.cat['qf']-q_predict)**2)

    def wrapper(self,params):
        q_predict = 1-self.failure_rate(self.cat['FIBERFLUX_G_EC'], self.cat['EFFTIME_QSO'], *params)
        return self.cost(q_predict)

    def wrapper_hist(self,params):
        h_predict = self.failure_rate_eff(self.bc, *params)
        diff = self.nzf-h_predict
        cost = np.sum((diff/self.nzfe)**2.)
        return cost

    def failure_rate(self,flux, efftime, a, b, c):
        #sn = flux * np.sqrt(efftime)
        #return np.clip(np.exp(-(sn+a)/b)+c/flux, 0, 1)
        return np.clip(np.exp(-(efftime+a)/b)+c/flux, 0, 1)

    def failure_rate_eff(self, efftime, a, b, c):
        #sn = flux * np.sqrt(efftime)
        #return np.clip(np.exp(-(sn+a)/b)+c/flux, 0, 1)
        return np.clip(np.exp(-(efftime+a)/b)+c, 0, 1)

    
    def hist_norm(self,fluxc):
        nzfper = []
        consl = []
        
        nb = 5
        pstep = 100//5
        costt = 0
        for i in range(0,nb):
            sel = self.cat['FIBERFLUX_G_EC'] > np.percentile(self.cat['FIBERFLUX_G_EC'],i*pstep)
            sel &= self.cat['FIBERFLUX_G_EC'] < np.percentile(self.cat['FIBERFLUX_G_EC'],(i+1)*pstep)
            mf = np.median(self.cat['FIBERFLUX_G_EC'][sel])
            if self.vis_5hist:
                print(mf)
            #fper.append(mf)
            wtf = (fluxc*(self.mft-self.cat['FIBERFLUX_G_EC'])/self.mft+1)*(self.wts_fid-1)+1
            selw = wtf < 1
            wtf[selw] = 1
            ha,_ = np.histogram(self.cat['EFFTIME_QSO'][sel],bins=self.bine)
            hf,_ = np.histogram(self.cat['EFFTIME_QSO'][sel&self.selgz],weights=wtf[sel&self.selgz],bins=self.bine)
            #if self.vis_5hist:
            #    print(mf)
            #    print(np.sum(ha))
            #    print(np.sum(hf))
            dl = hf/ha
            nzfper.append(dl)
            def ccost(c):
                return np.sum((dl-c)**2./self.nzfpere[i]**2.)
            resc = minimize(ccost, np.ones(1))
            bc = resc.x
            cost = ccost(bc)
            consl.append(bc)
            costt += cost
        if self.vis_5hist:
            for i in range(0,nb):
                plt.errorbar(self.bc,nzfper[i],self.nzfpere[i])
                plt.plot(self.bc,np.ones(len(self.bc))*consl[i],'k:')
            plt.show()
        return costt    
        
    
    def add_modpre(self,data,fn_root='',plot_color='k',plot_only=False,savefig=True):
        res = minimize(self.wrapper_hist, [-220, 130, 0.4], bounds=((-1000, 0), (0, 1000), (.25, .5)),
               method='Powell', tol=1e-6)
        pars = res.x
        chi2 = self.wrapper_hist(pars)
        print(pars,chi2)
        if self.cut == 'galactic':
            if self.cut_condition[0] == '<':
                label = 'b < %.2f deg' % float(self.cut_condition[1:])
            elif self.cut_condition[0] == '>':
                label = 'b > %.2f deg' % float(self.cut_condition[1:])
            plt.errorbar(self.bc,self.nzf,self.nzfe,fmt='ko',label=label,color=plot_color)
        else:
            plt.errorbar(self.bc,self.nzf,self.nzfe,fmt='ko',label='data',color=plot_color)
        mod = self.failure_rate_eff(self.bc, *pars)
        plt.plot(self.bc,mod,'--',label='model; chi2='+str(round(chi2,3)),color=plot_color)
        plt.ylabel('QSO Z failure rate')
        plt.xlabel('QSO EFFECTIVE exp time')
        plt.legend()
        if savefig:
            plt.savefig(fn_root+'overall_failratefit.png')
            plt.show()
        if not plot_only:
            gextc = 3.214
            rextc = 2.165
            dflux = data['FIBERFLUX_G']*10**(0.4*gextc*data['EBV']) #data['FIBERFLUX_G_EC']

            deff = 8.60/0.255 * data['TSNR2_QSO']#data['EFFTIME_ELG']
            #data['mod_success_rate'] = 1. -self.failure_rate(dflux,deff,*pars) 
            data['mod_success_rate'] = 1. -self.failure_rate_eff(deff,*pars)   
            assr = 1. -self.failure_rate_eff(self.cat['EFFTIME_QSO'],*pars)   
            relssr = assr/np.max(assr) 
            drelssr = data['mod_success_rate']/np.max(assr)#np.max(data['mod_success_rate'])
            seld = deff > 450
            seld &= deff < 1500
            print(len(relssr),len(drelssr[seld]),np.max(assr),np.max(data[seld]['mod_success_rate']))
            self.wts_fid = 1/relssr
            nzfper = []
            nzfpere = []
            fper = []
            self.mft = np.median(self.cat['FIBERFLUX_G_EC'])
            nb = 5
            pstep = 100//5
            for i in range(0,nb):
                sel = self.cat['FIBERFLUX_G_EC'] > np.percentile(self.cat['FIBERFLUX_G_EC'],i*pstep)
                sel &= self.cat['FIBERFLUX_G_EC'] < np.percentile(self.cat['FIBERFLUX_G_EC'],(i+1)*pstep)
                mf = np.median(self.cat['FIBERFLUX_G_EC'][sel])
                fper.append(mf)
                ha,_ = np.histogram(self.cat['EFFTIME_QSO'][sel],bins=self.bine)
                hf,_ = np.histogram(self.cat['EFFTIME_QSO'][sel&self.selgz],bins=self.bine)
                hfw,_ = np.histogram(self.cat['EFFTIME_QSO'][sel&self.selgz],weights=self.wts_fid[sel&self.selgz],bins=self.bine)
                nzfper.append(hf/ha)
                nzfpere.append(np.sqrt(ha-hf)/ha)
                #plt.plot(self.bc,hfw/ha)
            #plt.title('inputs')
            #plt.show()
            self.nzfpere = nzfpere    
            rest = minimize(self.hist_norm, np.ones(1))#, bounds=((-10, 10)),
                   #method='Powell', tol=1e-6)
            fcoeff = rest.x
            self.vis_5hist = True
            print(fcoeff,self.hist_norm(fcoeff))#,self.hist_norm(0.),self.hist_norm(1.)) 
            wtf = (fcoeff*(self.mft-dflux)/self.mft+1)*(1/drelssr-1)+1
            sel = wtf < 1
            wtf[sel] = 1
            data['WEIGHT_ZFAIL'] =  wtf
            return data

#             print(mf)
#             print(np.sum(ha))
#             print(np.sum(hf))
#             dl = hf/ha
#             plt.plot(self.bc,dl)
#         plt.show()
          
    
     