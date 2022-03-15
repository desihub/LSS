#Tools to study and correct for trends in spectroscopic succes rate (ssr)
#Initial LRG model fitting taken from Ronpgpu Zhou's notebook

import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
from scipy.optimize import curve_fit, minimize

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
    data['q'] = data['ZWARN']==0
    data['q'] &= data['Z']<1.5
    data['q'] &= data['DELTACHI2']>15  
    print('failure rate is '+str(np.sum(~data['q'])/len(data)))
    return data


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
    def __init__(self,specrel='fuji'):
        self.cat = get_LRG_data(specrel)

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
        dflux = data['FIBERFLUX_Z']*10**(0.4*1.211*data['EBV'])
        deff = 12.15 * data['TSNR2_LRG']
        data['mod_success_rate'] = 1. -self.failure_rate(dflux,deff,*pars)       
        return data
          
    
     