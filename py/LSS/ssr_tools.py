#Tools to study and correct for trends in spectroscopic succes rate (ssr)
#Initial LRG model fitting taken from Ronpgpu Zhou's notebook

import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
from scipy.optimize import curve_fit, minimize

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
        if sur == 'DA02':
            tfn+='zdone'
        fn = dir+tfn+'_full.dat.fits'    
        data = Table(fitsio.read(fn))
        data['q'] = data['o2c'] > 0.9
        cats.append(data)

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

class LRG_ssr:
    def __init__(self,specrel='fuji',efftime_min=500,efftime_max=2000):
        self.cat = get_LRG_data(specrel)
        mask = self.cat['EFFTIME_LRG']>efftime_min
        mask &= self.cat['EFFTIME_LRG']<efftime_max
        self.cat = self.cat[mask]

    def cost(self,q_predict):
        return np.sum((self.cat['qf']-q_predict)**2)

    def wrapper(self,params):
        q_predict = 1-self.failure_rate(self.cat['FIBERFLUX_Z_EC'], self.cat['EFFTIME_LRG'], *params)
        return self.cost(q_predict)

    def failure_rate(self,flux, efftime, a, b, c):
        sn = flux * np.sqrt(efftime)
        return np.clip(np.exp(-(sn+a)/b)+c/flux, 0, 1)

    def add_modpre(self,data):
        res = minimize(self.wrapper, [0, 10., 0.01], bounds=((-200, 200), (0, 100), (0., 1)),
               method='Powell', tol=1e-6)
        pars = res.x
        print(pars)
        dflux = data['FIBERFLUX_Z']*10**(0.4*1.211*data['EBV'])#data['FIBERFLUX_Z_EC']
        deff = 12.15 * data['TSNR2_LRG']#data['EFFTIME_LRG']
        data['mod_success_rate'] = 1. -self.failure_rate(dflux,deff,*pars)       
        return data

class ELG_ssr:
    def __init__(self,specrel='fuji',efftime_min=450,efftime_max=1500):
        self.cat = get_ELG_data_full('ELG_LOPnotqso')#get_ELG_data(specrel)
        mask = self.cat['EFFTIME_ELG']>efftime_min
        mask &= self.cat['EFFTIME_ELG']<efftime_max
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

    def fit_cons(dl,el,minv=0,step=0.01):
        c = minv
        newcost = np.sum((dl-c)**2./el**2.)
        oldcost = newcost + 1
        while newcost < oldcost:
            oc = c
            oldcost = newcost
            c += step
            newcost = np.sum((dl-c)**2./el**2.)
        return oldcost
    
    def hist_norm(self,fluxc):
        nzfper = []
        nzfpere = []
        fper = []
        consl = []
        mft = np.median(self.cat['FIBERFLUX_G_EC'])
        nb = 5
        pstep = 100//5
        costt = 0
        for i in range(0,nb):
            sel = self.cat['FIBERFLUX_G_EC'] > np.percentile(self.cat['FIBERFLUX_G_EC'],i*pstep)
            sel &= self.cat['FIBERFLUX_G_EC'] < np.percentile(self.cat['FIBERFLUX_G_EC'],(i+1)*pstep)
            mf = np.median(self.cat['FIBERFLUX_G_EC'][sel])
            fper.append(mf)
            wtf = fluxc*mft/mf*(self.wts_fid-1)+1
            ha,_ = np.histogram(self.cat['EFFTIME_ELG'][sel])
            hf,_ = np.histogram(self.cat['EFFTIME_ELG'][sel&self.selgz],weights=wtf)
            dl = hf[0]/ha[0]
            cost = fit_cons(dl,self.nzfpere[i])
            costt += cost
        return cost    
        
    
    def add_modpre(self,data):
        res = minimize(self.wrapper_hist, [-200, 10., 0.01], bounds=((-10000, 0), (0, 10000), (0., 1)),
               method='Powell', tol=1e-6)
        pars = res.x
        print(pars,self.wrapper_hist(pars))
        gextc = 3.214
        dflux = data['FIBERFLUX_G']*10**(0.4*gextc*data['EBV']) #data['FIBERFLUX_G_EC']
        deff = 8.60 * data['TSNR2_ELG']#data['EFFTIME_ELG']
        #data['mod_success_rate'] = 1. -self.failure_rate(dflux,deff,*pars) 
        data['mod_success_rate'] = 1. -self.failure_rate_eff(deff,*pars)   
        assr = 1. -self.failure_rate_eff(self.cat['EFFTIME_ELG'],*pars)   
        relssr = assr/np.max(assr) 
        drelssr = data['mod_success_rate']/np.max(data['mod_success_rate'])
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
            ha,_ = np.histogram(self.cat['EFFTIME_ELG'][sel])
            hf,_ = np.histogram(self.cat['EFFTIME_ELG'][sel&self.selgz])
            nzfper.append(hf[0]/ha[0])
            nzfpere.append(np.sqrt(ha[0]-hf[0])/ha[0])
        self.nzfpere = nzfpere    
        rest = minimize(self.hist_norm, 1, bounds=(-10, 10),
               method='Powell', tol=1e-6)
        fcoeff = rest.x
        print(fcoeff,self.hist_norm(fcoeff)) 
        
        data['WEIGHT_ZFAIL'] =  fluxc*mft/self.data[dflux]*(1/drelssr-1)+1
        return data
          
    
     