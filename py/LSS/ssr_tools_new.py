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
import time
from astropy.io import fits
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

def get_tsnr2z(tracer='ELG',night=20230128,expid=165078,tsnrdir='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/TSNR2z/',dz=0.001):
    '''
    get the relative template signal to noise ^2 for a given tracer type, night, and expid
    Copied from script from Julien Guy
    '''
    from desispec.io import read_frame
    import desisim.templates
    import scipy.ndimage
    ivar_table_filename = tsnrdir+"desi-"+str(night)+"-"+str(expid)+"-ivar.csv"
    if not os.path.isfile(ivar_table_filename) :
        spec_wave = None
        spec_ivar = []

        for spec in range(10) :

            wave_list=[]
            ivar_list=[]

            for arm in ["b","r","z"] :
                cam=f"{arm}{spec}"
                cframe_filename="/global/cfs/cdirs/desi/spectro/redux/daily/exposures/{}/{:08d}/cframe-{}-{:08d}.fits.gz".format(night,expid,cam,expid)
                if not os.path.isfile(cframe_filename) :
                    print("No file",cframe_filename,"?")
                    continue

                cframe = read_frame(cframe_filename)
                # directly get inverse variance from the sky fibers
                selection = np.where(cframe.fibermap["OBJTYPE"]=="SKY")[0]
                if len(selection)==0 :
                    print("No sky fibers?")
                    continue

                sw =  np.sum((cframe.ivar[selection]>0)*(cframe.mask[selection]==0),axis=0)
                median_ivar = np.sum(cframe.ivar[selection]*(cframe.mask[selection]==0),axis=0)/(sw+(sw==0))

                #plt.plot(cframe.wave,median_ivar)

                # store wavelength and ivar to merge them
                wave_list.append(cframe.wave)
                ivar_list.append(median_ivar)
            wave=np.unique(np.hstack(wave_list)) # this works because the wavelength of different cameras are aligned
            ivar=np.zeros(wave.size)
            for cam_wave,cam_ivar in zip(wave_list,ivar_list) :
                # we add the ivar of different cameras in overlapping area
                ivar += np.interp(wave,cam_wave,cam_ivar,left=0,right=0) # this works because the wavelength of different cameras are aligned


            if spec_wave is None :
                spec_wave = wave
            else :
                assert(np.all(spec_wave==wave))
            spec_ivar.append(ivar)
            #plt.plot(wave,ivar,"-")

        ivar=np.mean(np.array(spec_ivar),axis=0)
        #plt.plot(spec_wave,ivar,c="k")
        t=Table()
        t["WAVE"]=spec_wave[spec_wave<9820]
        t["IVAR"]=ivar[spec_wave<9820]
        t.write(ivar_table_filename,overwrite=True)

    t=Table.read(ivar_table_filename)
    wave=t["WAVE"]
    ivar=t["IVAR"]

    if tracer == 'ELG':
        # average an ensemble of ELG spectra at the same redshift and magnitude
        #if 1 :
        zmin=0.
        zmax=2.
        zref=1.
        wmin=wave[0]*(1+zmin)/(1+zref)
        wmax=wave[-1]*(1+zmax)/(1+zref)
        twave=np.linspace(wmin,wmax,int((wmax-wmin)/0.4))
        maker = desisim.templates.ELG(wave=twave)#, normfilter_south=normfilter_south)
        nmodel=50
        flux, other_wave, meta, objmeta = maker.make_templates(nmodel=nmodel, redshift=np.repeat(1.,nmodel), mag=np.repeat(22,nmodel), seed=12)
        tflux=np.mean(flux,axis=0)
        print(flux.shape)

        # use median filter to remove the continuum
        dw = twave[1]-twave[0]
        width = int(100./dw) #100A smoothing
        smooth_flux = scipy.ndimage.median_filter(tflux, width, mode='constant')
        dflux = tflux-smooth_flux

        # compute tsnr2 on a redshift range
        redshift=np.linspace(0,1.7,int(1.700001/dz)+1)
        tsnr2=np.zeros(redshift.shape)
        for i,z in enumerate(redshift):
            tmp=np.interp(wave,twave*(1+z)/(1+zref),dflux)
            tsnr2[i]=np.sum(ivar*tmp**2)
        #print(redshift[:20])
        return redshift,tsnr2
    else:
        print('ERROR, only ELG tracer type supported so far')




def normed_linfit(data,seltot,sel,range=(80,200),col='TSNR2_ELG',cl='k',ps='o',lab=''):
    a,bins = np.histogram(data[seltot][col],range=range)
    b,_ = np.histogram(data[sel][col],bins=bins)
    sp = bins[1]-bins[0]
    err = np.sqrt(b)/a #rough, correct formula is to use binomial
    normb = np.sum(b)/np.sum(a)
    cchi = np.sum((b/a/normb-1.)**2./(err/normb)**2.)
    bc = bins[:-1]+sp/2.
    #plt.errorbar(bc,b/a/normb,err/normb,fmt=cl+ps,label=lab+str(round(cchi,3)))
    res = np.polyfit(bc, b/a/normb,1, w=1./(err/normb))
    #chi2lin = np.sum((b/a/normb-(res[0]*bc+res[1]))**2./(err/normb)**2.)
    #plt.plot(bc,res[0]*bc+res[1],'r--',label=r'$\chi^2_{\rm lin}=$'+str(round(chi2lin,3)))
    return res

class model_ssr_zfac:
    def __init__(self,input_data,tsnr_min=80,tsnr_max=200,tracer='ELG',reg=None,outdir='',band='G',outfn_root='test',readpars=False):
        self.cat = input_data

        mask = self.cat['TSNR2_'+tracer]>tsnr_min
        mask &= self.cat['TSNR2_'+tracer]<tsnr_max
        self.tsnr_max = tsnr_max
        self.reg = reg
        if reg is not None:
            mask &= self.cat['PHOTSYS'] == reg
        self.cat = self.cat[mask]
        self.selgz = common.goodz_infull(tracer,self.cat,zcol='Z_not4clus')
        self.ztsnr,self.tsnr2z = get_tsnr2z(tracer)
        modo2 = (1500/self.tsnr2z-1)/(15/((self.ztsnr-1)))+1
        for i in range(0,len(self.ztsnr)):
            self.ztsnr[i] = round(self.ztsnr[i],3)
        self.modo2_dict = dict(zip(self.ztsnr,modo2))
        self.relzfac = self.get_relzfac(self.cat)
        if tracer == 'ELG':
            minz = 0.8
            maxz = 1.6
        self.selz = self.cat['Z_not4clus'] > minz
        self.selz &= self.cat['Z_not4clus'] < maxz
        self.maxz = maxz
        
        self.res_mod_slp = self.get_slpfunc()

        if not os.path.exists(outdir+outfn_root+'slp_wzfac.txt'):
            fo = open(outdir+outfn_root+'slp_wzfac.txt','w')
            fo.write(str(self.res_mod_slp))
            fo.close()
        if not os.path.exists(outdir+outfn_root+'slp_vszfac.txt'):
            fo = open(outdir+outfn_root+'slp_vszfac.txt','w')
            for i in range(0,len(self.slpl)):
                fo.write(str(self.zfacl[i])+' '+str(self.slpl[i])+'\n')
            fo.close()
        #fo.write('#a b c chi2\n')
        #for par in pars:
        #   fo.write(str(par)+' ')
        #fo.write(str(chi2)+'\n')

        self.ssrtot = len(self.cat[self.selgz])/len(self.cat)
        self.tracer = tracer
            
 
    def get_relzfac(self,data,zcol='Z_not4clus',zmax=1.6):           
        relzfac = []
        for z in data[zcol]:
            rz = round(z,3)
            if rz >0 and rz < zmax:
                relzfac.append(self.modo2_dict[rz])
            else:
                relzfac.append(1)
        relzfac = np.array(relzfac)
        return relzfac


    
    def get_slpfunc(self,pstep=5):
        slpl = []
        zfacl = []
        for i in range(0,100//pstep):
            minzfac = np.percentile(self.relzfac[self.selz],i*pstep)
            maxzfac = np.percentile(self.relzfac[self.selz],(i+1)*pstep)
            selzfac = self.relzfac > minzfac
            selzfac &= self.relzfac < maxzfac
            seltot = np.ones(len(self.cat),dtype='bool')
            sel = self.selz & self.selgz & selzfac
            res = normed_linfit(self.cat,seltot,sel)
            slpl.append(res[0])
            zfacl.append(np.median(self.relzfac[sel]))
        self.slpl = slpl
        self.zfacl = zfacl
        res = np.polyfit(zfacl,slpl,1)
        return res

    def add_modpre(self,data,zcol='Z'):
        relzfac = self.get_relzfac(data,zcol=zcol)
        mod_slp = self.res_mod_slp[1]+self.res_mod_slp[0]*relzfac
        sel_neg = mod_slp < 0
        mod_slp[sel_neg] = 0
        mtsnr2 = np.median(self.cat['TSNR2_'+self.tracer])
        mod_relssr = mod_slp*data['TSNR2_'+self.tracer]+1-mod_slp*mtsnr2
        mod_ssr = np.clip(self.ssrtot*mod_relssr,0,1)
        wtf = 1./mod_relssr

        return wtf,mod_ssr




class model_ssr:
    '''
    This class will fit a model based on TSNR2_<type> and FIBERFLUX_<band> for the redshift success and produce weights
    that are the inverse of the relative predicted redshift success
    '''
    def __init__(self,input_data,tsnr_min=80,tsnr_max=200,tracer='ELG',reg=None,outdir='',band='G',outfn_root='test',readpars=False,overwrite_pars_ssrmaxflux=True):
        self.cat = input_data

        mask = self.cat['TSNR2_'+tracer]>tsnr_min
        mask &= self.cat['TSNR2_'+tracer]<tsnr_max
        self.tsnr_max = tsnr_max
        self.reg = reg
        if reg is not None:
            mask &= self.cat['PHOTSYS'] == reg
        self.cat = self.cat[mask]
        if tracer == 'QSO':
            names = list(self.cat.dtype.names)
            if 'OII_FLUX' not in names:
                emline = fits.open('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/emlin_catalog.fits')[1].data
        
                ss = np.searchsorted(sorted(emline['TARGETID']), self.cat['TARGETID'])
        
                ass = np.argsort(emline['TARGETID'])
        
                oii_flux = emline['OII_FLUX'][ass][ss]
                oii_flux_ivar = emline['OII_FLUX_IVAR'][ass][ss]
        
                oiii_flux = emline['OIII_FLUX'][ass][ss]
                oiii_flux_ivar = emline['OIII_FLUX_IVAR'][ass][ss]
            else:
                oii_flux = self.cat['OII_FLUX']
                oii_flux_ivar = self.cat['OII_FLUX_IVAR']
        
                oiii_flux = self.cat['OIII_FLUX']
                oiii_flux_ivar = self.cat['OIII_FLUX_IVAR']
           
            dchi_cut =40 #was 30
            o2c_cut = 1.2 #was 0.9
            oiii_cut = 5
            selgal = (( ~common.goodz_infull(tracer,self.cat,zcol='Z_not4clus') & (self.cat['Z_RR'] > 0.01) ) & 
                ((self.cat['DELTACHI2'] > dchi_cut) | (np.log10(oii_flux * oii_flux_ivar**0.5) > o2c_cut - 0.2 * np.log10(self.cat['DELTACHI2']))
                | (oiii_flux * oiii_flux_ivar**0.5 > oiii_cut))
                )  #...now 40,1.2,5...earlier 30,0.9,5
            selstar = ( ~common.goodz_infull(tracer,self.cat,zcol='Z_not4clus') & (self.cat['Z_RR'] < 0.01))     
            self.cat=self.cat[~selgal & ~selstar]
        
        print(len(self.cat))
        self.cat['FIBERFLUX_'+band+'_EC'] = self.cat['FIBERFLUX_'+band]*10**(0.4*extdict[band]*self.cat['EBV'])
        self.selgz = common.goodz_infull(tracer,self.cat,zcol='Z_not4clus')
        ha,bine = np.histogram(self.cat['TSNR2_'+tracer])
        medt = np.median(self.cat['TSNR2_'+tracer])
        hf,_ = np.histogram(self.cat['TSNR2_'+tracer][~self.selgz],bins=bine)
        self.nzf = hf/ha
        tot_failrate = np.sum(hf)/np.sum(ha)
        high_failrate = np.sum(hf[5:])/np.sum(ha[5:])
        print(self.nzf)
        self.nzfe = np.sqrt(hf)/ha
        bc = []
        median_bins = []
        bs = bine[1]-bine[0]
        for i in range(0,len(bine)-1):
            bc.append(bine[i]+bs/2.) 
            median_bins.append(np.median(self.cat['TSNR2_'+tracer][(self.cat['TSNR2_'+tracer] >= bine[i]) & (self.cat['TSNR2_'+tracer] < bine[i+1])]))
        self.bc = np.array(bc)
        self.median_bins = np.array(median_bins)
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
            maxv = tot_failrate
            if high_failrate > tot_failrate:
                maxv = 1.01*high_failrate
            bguess = 10*medt/120.
            aguess = -16*medt/120.
            if tracer == 'QSO':
                #DR1 S parameters
                aguess = -3.1
                bguess = 6.33
                high_failrate = 0.33
            res = minimize(self.wrapper_hist, [aguess, bguess, high_failrate], bounds=((-2*tsnr_max, 2*tsnr_max), (0.001, tsnr_max), (0., maxv)),method='Powell')#,
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
        plt.clf()
        plt.errorbar(self.bc,self.nzf,self.nzfe,fmt='ko',label='data')
        mod = self.failure_rate_eff(self.bc, *pars)
        plt.plot(self.bc,mod,'k--',label='model; chi2='+str(round(chi2,3)))
        plt.ylabel(outfn_root+rw+' Z failure rate')
        plt.xlabel('TSNR2_'+tracer)
        plt.legend()
        plt.savefig(self.outdir+outfn_root+rw+'overall_failratefit.png')        
        #plt.show()
        plt.clf()
        #fit to fiberflux trend
        print('self.cat[TSNR2_tracer]',self.cat['TSNR2_'+tracer])
        print('self.cat',self.cat)
        assr = 1. -self.failure_rate_eff(self.cat['TSNR2_'+tracer],*pars)   
        print('assr',assr)
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
            fcoeff,piv,C = parsflux[0],parsflux[1],parsflux[2]
            ssrvsflux = np.loadtxt(self.outdir+outfn_root+rw+'maxssrvsflux.txt').transpose()
            self.mfl = ssrvsflux[0]
            self.consl = ssrvsflux[1]
            
        else:
            if tracer == 'LRG':
               self.fluxfittype = 'piecewise'
               print('about to do linear fit for flux dependence')
               rest = minimize(lambda params: self.hist_norm(params, fluxfittype='linear'), [2,self.mft],method='Powell')
               fcoeff_start, piv_start = rest.x
               #if reg == 'N':
               #print('reg',reg)
               print('linear fit done, now doing piecewise')
               rest = minimize(lambda params: self.hist_norm(params, fluxfittype=self.fluxfittype), [fcoeff_start, piv_start, 3],method='Powell')
               #elif reg == 'S':
               #    print('reg',reg)
               #    rest = minimize(lambda params: self.hist_norm(params, fluxfittype=self.fluxfittype), [14.82,2.94,2.09],method='Powell')
               #rest = minimize(self.hist_norm, [2,self.mft,3],method='Powell',args='linear')#np.ones(1))#, bounds=((-10, 10)),
               #method='Powell', tol=1e-6)
               fcoeff,piv,C = rest.x
               self.vis_5hist = True
               chi2 = self.hist_norm([fcoeff,piv,C],fluxfittype=self.fluxfittype,outfn=outfn_root+rw+'5histfit.png')
               print('results and chi2:')
               print(fcoeff,piv,C,chi2)#,self.hist_norm(0.),self.hist_norm(1.)) 
               fo = open(self.outdir+outfn_root+rw+'pars_fluxfit.txt','w')
               fo.write('#'+self.band+'flux fit\n')
               fo.write('#coeff flux_pivot C chi2\n')
        
               fo.write(str(fcoeff)+' '+str(piv)+' '+str(C)+' ')
               fo.write(str(chi2)+'\n')
               fo.close()
            else:
               self.fluxfittype = 'linear'
               print('about to do linear fit for flux dependence')
               rest = minimize(lambda params: self.hist_norm(params, fluxfittype=self.fluxfittype), [2,self.mft],method='Powell')
               #rest = minimize(self.hist_norm, [2,self.mft],method='Powell')#np.ones(1))#, bounds=((-10, 10)),
               #method='Powell', tol=1e-6)
               fcoeff,piv = rest.x
               self.vis_5hist = True
               chi2 = self.hist_norm([fcoeff,piv],fluxfittype=self.fluxfittype,outfn=outfn_root+rw+'5histfit.png')
               print('results and chi2:')
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
        if tracer == 'LRG':
            self.C = C
            #print(self.mfl)
        
        #Now, we need a smooth function for maximum ssr vs. flux
        if readpars == 'junk':
            parsmaxflux = np.loadtxt(self.outdir+outfn_root+rw+'pars_ssrmaxflux.txt')
            #if tracer == 'ELG':
            #    self.flux_mod = np.poly1d(parsmaxflux)
            #else:
            self.pars_ferf = parsmaxflux
            self.flux_mod = self.ssrvflux_erf
                
        else:
            if overwrite_pars_ssrmaxflux:
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
            flux_max = np.percentile(self.cat['FIBERFLUX_'+self.band+'_EC'],95)
            flux_min = np.min(self.cat['FIBERFLUX_'+self.band+'_EC'])

            a = np.histogram(self.cat['FIBERFLUX_'+self.band+'_EC'][self.selgz],weights=wtf[self.selgz],bins=20,range=(flux_min,flux_max))
            b = np.histogram(self.cat['FIBERFLUX_'+self.band+'_EC'],bins=a[1])
            self.ssr_flux = a[0]/b[0]
            self.flux_vals = a[1][:-1]+(a[1][1]-a[1][0])/2

            ssrvflux = minimize(self.wrapper_ssrvflux,[self.consl[-1],self.mfl[0],self.mfl[-1]],method='Powell')
            self.pars_ferf = ssrvflux.x
            print(self.pars_ferf)
            self.flux_mod = self.ssrvflux_erf
            if overwrite_pars_ssrmaxflux:
                for par in self.pars_ferf :
                    fo.write(str(par)+' ')
                fo.write('\n')    
                fo.close()
            plt.plot(self.mfl,self.consl,'rd')
            plt.plot(self.mfl,self.flux_mod(self.mfl),'r--')
            plt.plot(self.flux_vals,self.ssr_flux,'ko')
            plt.plot(self.flux_vals,self.flux_mod(self.flux_vals),'k-')
            plt.savefig(self.outdir+outfn_root+rw+'flux_dep.png') 
            #plt.show()
           
            
        
        
    def ssrvflux_erf(self,flux):
        return self.pars_ferf[0]*erf((self.pars_ferf[1]+flux)/self.pars_ferf[2])
    
    def wrapper_ssrvflux(self,params):
        #mod = gen_erf(self.mfl,*params)
        #cost = np.sum((self.consl-mod)**2.)
        mod = gen_erf(self.flux_vals,*params)
        cost = np.sum((self.ssr_flux-mod)**2.)        
        return cost
    
    def wrapper_hist(self,params):
        h_predict = self.failure_rate_eff(self.median_bins, *params)
        diff = self.nzf-h_predict
        cost = np.sum((diff/self.nzfe)**2.)
        return cost


    def failure_rate_eff(self, efftime, a, b, c):
        #sn = flux * np.sqrt(efftime)
        #return np.clip(np.exp(-(sn+a)/b)+c/flux, 0, 1)
        return np.clip(np.exp(-(efftime+a)/b)+c, 0, 1)

    
    def hist_norm(self,params,fluxfittype='linear',outfn='test.png'):
        if (fluxfittype != 'linear') and (fluxfittype != 'piecewise'):
            print('ERROR, fluxfittype must be either linear or piecewise')
        #print('call to hist norm, params',params)
        #print('len of mod.cat',len(self.cat['FIBERFLUX_'+self.band+'_EC']))
        t0 = time.time()
        nzfper = []
        consl = []
        if fluxfittype == 'piecewise':
            fluxc,piv_flux,C = params
            flux_break = 3.
            self.flux_break = flux_break
        elif fluxfittype == 'linear':
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
            
            if fluxfittype == 'piecewise':
                rel_flux = self.cat['FIBERFLUX_'+self.band+'_EC']/piv_flux #mod.mft
                flux_piece = fluxc * (1 - rel_flux) + 1
            
                #flux_break = 3.
                D = (fluxc + 1 - flux_break * fluxc/piv_flux - C) * 1./flux_break
                flux_piece[self.cat['FIBERFLUX_'+self.band+'_EC'] > flux_break] = C + D * self.cat['FIBERFLUX_'+self.band+'_EC'][self.cat['FIBERFLUX_'+self.band+'_EC'] > flux_break]
            
                #inv_rel_flux[mod.cat['FIBERFLUX_'+mod.band+'_EC'] > piv_flux2] = - 0.5 * (rel_flux[mod.cat['FIBERFLUX_'+mod.band+'_EC'] > piv_flux2] - piv_flux2/piv_flux) + (1 - piv_flux2/piv_flux)
                wtf = flux_piece*(self.wts_fid-1)+1
            
            elif fluxfittype == 'linear':
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
            #print('cost',cost)
        if self.vis_5hist:
            for i in range(0,nb):
                plt.errorbar(self.bc,nzfper[i],self.nzfpere[i])
                plt.plot(self.bc,np.ones(len(self.bc))*consl[i],'k:')
            plt.ylabel(self.outfn_root+' Z success rate, in fiber bins')
            plt.xlabel('TSNR2_'+self.tracer)
            plt.legend()
            plt.savefig(self.outdir+outfn)
            plt.clf()        

            #plt.show()
            self.consl = consl
            self.mfl = mfl
            
        #print('time for hist_norm',time.time()-t0)
        #print('costt',costt)
        return costt    
        
    
    def add_modpre(self,data):
        dflux = data['FIBERFLUX_'+self.band]*10**(0.4*extdict[self.band]*data['EBV']) #data['FIBERFLUX_G_EC']
        deff = data['TSNR2_'+self.tracer]#data['EFFTIME_ELG']
        #data['mod_success_rate'] = 1. -self.failure_rate(dflux,deff,*pars) 
        tssr = 1.-self.failure_rate_eff(deff,*self.pars)
        sel = tssr == 0
        tssr[sel] = .01

        max_tssr = 1. - self.failure_rate_eff(self.tsnr_max,*self.pars)
        relssr = tssr/max_tssr
        max_ssr_flux = self.flux_mod(dflux) 
        print(np.min(max_ssr_flux),np.max(max_ssr_flux),np.mean(max_ssr_flux))
        #data['mod_success_rate'] = 1. -   
        rel_flux = dflux/self.piv
        
        
        if self.fluxfittype == 'piecewise':
            flux_piece = self.fcoeff * (1 - rel_flux) + 1
            
            #flux_break = 3.
            D = (self.fcoeff + 1 - self.flux_break * self.fcoeff/self.piv - self.C) * 1./self.flux_break
            flux_piece[dflux > self.flux_break] = self.C + D * dflux[dflux > self.flux_break]
            
            #inv_rel_flux[mod.cat['FIBERFLUX_'+mod.band+'_EC'] > piv_flux2] = - 0.5 * (rel_flux[mod.cat['FIBERFLUX_'+mod.band+'_EC'] > piv_flux2] - piv_flux2/piv_flux) + (1 - piv_flux2/piv_flux)
            wtf = flux_piece*(1/relssr-1)+1
            
        elif self.fluxfittype == 'linear':
            #print('rel_flux',np.shape(rel_flux))
            #print('wts_fid',np.shape(self.wts_fid))
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

