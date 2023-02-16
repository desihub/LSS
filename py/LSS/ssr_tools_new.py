#Tools to study and correct for trends in spectroscopic succes rate (ssr)
#Initial LRG model fitting taken from Ronpgpu Zhou's notebook

import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
from scipy.optimize import curve_fit, minimize
from scipy.special import erf
from astropy.coordinates import SkyCoord
from astropy import units as u

import LSS.common_tools as common

extdict ={'G':3.214,'R':2.165,'Z':1.211}

def gen_erf(val,a,b,c):
    return a*erf((b+val)/c)


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




class model_ssr:
    def __init__(self,input_data,tsnr_min=80,tsnr_max=200,tracer='ELG',reg=None,outdir='',band='G',outfn_root='test',readpars=False):
        self.cat = input_data

        mask = self.cat['TSNR2_'+tracer]>tsnr_min
        mask &= self.cat['TSNR2_'+tracer]<tsnr_max
        self.tsnr_max = tsnr_max
        self.reg = reg
        if reg is not None:
            mask &= self.cat['PHOTSYS'] == reg
        self.cat = self.cat[mask]
        print(len(self.cat))
        self.cat['FIBERFLUX_'+band+'_EC'] = self.cat['FIBERFLUX_'+band]*10**(0.4*extdict[band]*self.cat['EBV'])
        self.selgz = common.goodz_infull(tracer,self.cat,zcol='Z_not4clus')
        ha,bine = np.histogram(self.cat['TSNR2_'+tracer])
        hf,_ = np.histogram(self.cat['TSNR2_'+tracer][~self.selgz],bins=bine)
        self.nzf = hf/ha
        tot_failrate = np.sum(hf)/np.sum(ha)
        high_failrate = np.sum(hf[5:])/np.sum(ha[5:])
        print(self.nzf)
        self.nzfe = np.sqrt(hf)/ha
        bc = []
        bs = bine[1]-bine[0]
        for i in range(0,len(bine)-1):
            bc.append(bine[i]+bs/2.) 
        self.bc = np.array(bc)
        self.bine = bine
        self.vis_5hist = False
        self.outdir = outdir
        self.band = band
        self.tracer = tracer
        self.outfn_root = outfn_root
        rw = ''
        if self.reg is not None:
            rw = self.reg
        
        #fit to TSNR2
        
        if readpars:
            parsf = np.loadtxt(self.outdir+outfn_root+rw+'pars_overall.txt')
            pars = np.array([parsf[0],parsf[1],parsf[2]])
            chi2 = parsf[3]
        else:
            res = minimize(self.wrapper_hist, [-16, 10., high_failrate], bounds=((-2*tsnr_max, 2*tsnr_max), (0.001, tsnr_max), (0., tot_failrate)),method='Powell')#,
               #method='Powell', tol=1e-6)
            pars = res.x
            chi2 = self.wrapper_hist(pars)
            print(pars,chi2)
            fo = open(self.outdir+outfn_root+rw+'pars_overall.txt','w')
            fo.write('#overall fit\n')
            fo.write('#a b c chi2\n')
            for par in pars:
                fo.write(str(par)+' ')
            fo.write(str(chi2)+'\n')
            fo.close()
        self.pars = pars
        plt.errorbar(self.bc,self.nzf,self.nzfe,fmt='ko',label='data')
        mod = self.failure_rate_eff(self.bc, *pars)
        plt.plot(self.bc,mod,'k--',label='model; chi2='+str(round(chi2,3)))
        plt.ylabel(outfn_root+rw+' Z failure rate')
        plt.xlabel('TSNR2_'+tracer)
        plt.legend()
        plt.savefig(self.outdir+outfn_root+rw+'overall_failratefit.png')        
        plt.show()
        plt.clf()
        #fit to fiberflux trend
        assr = 1. -self.failure_rate_eff(self.cat['TSNR2_'+tracer],*pars)   
        relssr = assr/np.max(assr) 
        self.wts_fid = 1/relssr
        nzfper = []
        nzfpere = []
        fper = []
        self.mft = np.median(self.cat['FIBERFLUX_'+self.band+'_EC'])
        nb = 5
        pstep = 100//5
        for i in range(0,nb):
            sel = self.cat['FIBERFLUX_'+self.band+'_EC'] > np.percentile(self.cat['FIBERFLUX_'+self.band+'_EC'],i*pstep)
            sel &= self.cat['FIBERFLUX_'+self.band+'_EC'] < np.percentile(self.cat['FIBERFLUX_'+self.band+'_EC'],(i+1)*pstep)
            mf = np.median(self.cat['FIBERFLUX_'+self.band+'_EC'][sel])
            fper.append(mf)
            ha,_ = np.histogram(self.cat['TSNR2_'+tracer][sel],bins=self.bine)
            hf,_ = np.histogram(self.cat['TSNR2_'+tracer][sel&self.selgz],bins=self.bine)
            hfw,_ = np.histogram(self.cat['TSNR2_'+tracer][sel&self.selgz],weights=self.wts_fid[sel&self.selgz],bins=self.bine)
            nzfper.append(hf/ha)
            sel = ha == hf
            hf[sel] -= 1 #so that the errors aren't 0
            err = np.sqrt(ha*(1-hf/ha))/ha
            nzfpere.append(err)
        self.nzfpere = nzfpere
        #print(nzfpere)    
        
        if readpars:
            parsflux = np.loadtxt(self.outdir+outfn_root+rw+'pars_fluxfit.txt')
            fcoeff,piv = parsflux[0],parsflux[1]
            ssrvsflux = np.loadt(self.outdir+outfn_root+rw+'maxssrvsflux.txt').transpose()
            self.mfl = ssrvsflux[0]
            self.consl = ssrvsflux[1]
            
        else:
            rest = minimize(self.hist_norm, [2,self.mft],method='Powell')#np.ones(1))#, bounds=((-10, 10)),
               #method='Powell', tol=1e-6)
            fcoeff,piv = rest.x
            self.vis_5hist = True
            chi2 = self.hist_norm([fcoeff,piv])
            print(fcoeff,piv,chi2)#,self.hist_norm(0.),self.hist_norm(1.)) 
            fo = open(self.outdir+outfn_root+rw+'pars_fluxfit.txt','w')
            fo.write('#'+self.band+'flux fit\n')
            fo.write('#coeff flux_pivot chi2\n')
        
            fo.write(str(fcoeff)+' '+str(piv)+' ')
            fo.write(str(chi2)+'\n')
            fo.close()
            self.mfl = np.array(self.mfl)
            print(self.consl)
            fo = open(self.outdir+outfn_root+rw+'maxssrvsflux.txt','w')
            fo.write('#flux max_ssr\n')
            for i in range(0,len(self.mfl)):
                fo.write(str(self.mfl[i])+' '+str(self.consl[i])+'\n')
            fo.close()
            
        self.fcoeff = fcoeff
        self.piv = piv
            #print(self.mfl)
        
        #Now, we need a smooth function for maximum ssr vs. flux
        if readpars:
            parsmaxflux = np.loadtxt(self.outdir+outfn_root+rw+'pars_ssrmaxflux.txt')
            #if tracer == 'ELG':
            #    self.flux_mod = np.poly1d(parsmaxflux)
            #else:
            self.pars_ferf = parsmaxflux
            self.flux_mod = self.ssrvflux_erf
                
        else:
            fo = open(self.outdir+outfn_root+rw+'pars_ssrmaxflux.txt','w')
            fo.write('#fit parameters for maximum ssr as a function of flux\n')
        

            #if tracer == 'ELG':
            #    flux_par = np.polyfit(np.array(self.mfl),np.array(self.consl),2)
            #    print(flux_par)
            #    self.flux_mod = np.poly1d(flux_par)
            #    for par in flux_par :
            #        fo.write(str(par)+' ')
            #    fo.write('\n')    

            #else:
            #we expect asymptotic behavior for LRG and BGS
            rel_flux = self.cat['FIBERFLUX_'+self.band+'_EC']/self.piv#self.mft
            wtf = (self.fcoeff*(1-rel_flux)+1)*(self.wts_fid-1)+1
            selw = wtf < 1
            wtf[selw] = 1
            flux_max = np.percentile(self.cat['FIBERFLUX_'+self.band+'_EC'],99)
            flux_min = np.min(self.cat['FIBERFLUX_'+self.band+'_EC'])

            a = np.histogram(self.cat['FIBERFLUX_'+self.band+'_EC'][self.selgz],weights=wtf[self.selgz],bins=20,range=(flux_min,flux_max))
            b = np.histogram(self.cat['FIBERFLUX_'+self.band+'_EC'],bins=a[1])
            self.ssr_flux = a[0]/b[0]
            self.flux_vals = a[1][:-1]+(a[1][1]-a[1][0])/2

            ssrvflux = minimize(self.wrapper_ssrvflux,[self.consl[-1],self.mfl[0],self.mfl[-1]],method='Powell')
            self.pars_ferf = ssrvflux.x
            print(self.pars_ferf)
            self.flux_mod = self.ssrvflux_erf
            for par in self.pars_ferf :
                fo.write(str(par)+' ')
            fo.write('\n')    
            fo.close()
            plt.plot(self.mfl,self.consl,'rd')
            plt.plot(self.mfl,self.flux_mod(self.mfl),'r--')
            plt.plot(self.flux_vals,self.ssr_flux,'ko')
            plt.plot(self.flux_vals,self.flux_mod(self.flux_vals),'k-')
            plt.show()
           
            
        
        
    def ssrvflux_erf(self,flux):
        return self.pars_ferf[0]*erf((self.pars_ferf[1]+flux)/self.pars_ferf[2])
    
    def wrapper_ssrvflux(self,params):
        #mod = gen_erf(self.mfl,*params)
        #cost = np.sum((self.consl-mod)**2.)
        mod = gen_erf(self.flux_vals,*params)
        cost = np.sum((self.ssr_flux-mod)**2.)        
        return cost
    
    def wrapper_hist(self,params):
        h_predict = self.failure_rate_eff(self.bc, *params)
        diff = self.nzf-h_predict
        cost = np.sum((diff/self.nzfe)**2.)
        return cost


    def failure_rate_eff(self, efftime, a, b, c):
        #sn = flux * np.sqrt(efftime)
        #return np.clip(np.exp(-(sn+a)/b)+c/flux, 0, 1)
        return np.clip(np.exp(-(efftime+a)/b)+c, 0, 1)

    
    def hist_norm(self,params,outfn='test.png'):
        nzfper = []
        consl = []
        fluxc,piv_flux = params
        nb = 5
        pstep = 100//5
        costt = 0
        mfl = []
        for i in range(0,nb):
            sel = self.cat['FIBERFLUX_'+self.band+'_EC'] > np.percentile(self.cat['FIBERFLUX_'+self.band+'_EC'],i*pstep)
            sel &= self.cat['FIBERFLUX_'+self.band+'_EC'] < np.percentile(self.cat['FIBERFLUX_'+self.band+'_EC'],(i+1)*pstep)
            mf = np.median(self.cat['FIBERFLUX_'+self.band+'_EC'][sel])
            if self.vis_5hist:
                print(mf)
                mfl.append(mf)
            #fper.append(mf)
            
            rel_flux = self.cat['FIBERFLUX_'+self.band+'_EC']/piv_flux#self.mft
            wtf = (fluxc*(1-rel_flux)+1)*(self.wts_fid-1)+1
            selw = wtf < 1
            wtf[selw] = 1
            ha,_ = np.histogram(self.cat['TSNR2_'+self.tracer][sel],bins=self.bine)
            hf,_ = np.histogram(self.cat['TSNR2_'+self.tracer][sel&self.selgz],weights=wtf[sel&self.selgz],bins=self.bine)
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
            consl.append(bc[0])
            costt += cost
        if self.vis_5hist:
            for i in range(0,nb):
                plt.errorbar(self.bc,nzfper[i],self.nzfpere[i])
                plt.plot(self.bc,np.ones(len(self.bc))*consl[i],'k:')
            plt.ylabel(self.outfn_root+' Z success rate, in fiber bins')
            plt.xlabel('TSNR2_'+self.tracer)
            plt.legend()
            plt.savefig(self.outdir+outfn)        

            plt.show()
            self.consl = consl
            self.mfl = mfl
        return costt    
        
    
    def add_modpre(self,data):
        dflux = data['FIBERFLUX_'+self.band]*10**(0.4**extdict[self.band]*data['EBV']) #data['FIBERFLUX_G_EC']
        deff = data['TSNR2_'+self.tracer]#data['EFFTIME_ELG']
        #data['mod_success_rate'] = 1. -self.failure_rate(dflux,deff,*pars) 
        tssr = 1.-self.failure_rate_eff(deff,*self.pars)
        max_tssr = 1. - self.failure_rate_eff(self.tsnr_max,*self.pars)
        relssr = tssr/max_tssr
        max_ssr_flux = self.flux_mod(dflux) 
        print(np.min(max_ssr_flux),np.max(max_ssr_flux),np.mean(max_ssr_flux))
        #data['mod_success_rate'] = 1. -   
        rel_flux = dflux/self.piv
        wtf = (self.fcoeff*(1-rel_flux)+1)*(1/relssr-1)+1
        
        sel = wtf < 1
        wtf[sel] = 1
        mod = max_ssr_flux/wtf
        #data['WEIGHT_ZFAIL'] =  wtf
        return wtf,mod

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

