from __future__ import print_function

import os, sys
import numpy as np
import pylab as plt
import healpy as hp
from astropy.io import fits
from astropy.table import Table
from iminuit import Minuit
#import iminuit.frontends
from scipy.optimize import minimize

class Syst:

    def __init__(self, data_we, rand_we):
        self.data_we = data_we 
        self.rand_we = rand_we
        self.data_syst = {}
        self.rand_syst = {}
        self.syst_names = []
        self.nsyst = 0
        self.ndata = data_we.size
        self.nrand = rand_we.size
        #self.factor = np.sum(self.rand_we)/np.sum(self.data_we) #AJR put this here because we still want normalization based on total sample
        #print(self.factor)

    def add_syst(self, name, data_syst, rand_syst):
        assert( data_syst.size == self.ndata ) 
        assert( rand_syst.size == self.nrand )
        assert( name not in self.syst_names )
        #-- checking bad values of systematics
        wd = np.isnan(data_syst)|np.isinf(data_syst)
        data_syst[wd] = hp.UNSEEN 
        wr = np.isnan(rand_syst)|np.isinf(rand_syst)
        rand_syst[wr] = hp.UNSEEN 

        self.syst_names.append(name)
        self.data_syst[name] = data_syst
        self.rand_syst[name] = rand_syst
        self.nsyst += 1

    def cut_outliers(self, p=1., verbose=False):
        ''' Cut galaxies and randoms with extreme values of systematics '''

        if p==0:
            return 

        #print(f'Cutting {p}% of extreme values in each map:') 
        w_data = np.ones(self.ndata, dtype=bool)
        w_rand = np.ones(self.nrand, dtype=bool) 
        for name in self.syst_names:
            data_syst = self.data_syst[name]
            rand_syst = self.rand_syst[name]
    
            w = data_syst!=hp.UNSEEN
            syst_min = np.percentile(data_syst[w], p/2) #22.4
            syst_max = np.percentile(data_syst[w], 100-p/2) #23.8

            w_data &= (data_syst >= syst_min) & \
                      (data_syst <= syst_max)
            w_rand &= (rand_syst >= syst_min) & \
                      (rand_syst <= syst_max)
            if verbose:
                print(' ', name, 'from', syst_min, 'to', syst_max)
            
        
        #-- Applying cuts and updating arrays
        for name in self.syst_names:
            self.data_syst[name] = self.data_syst[name][w_data] 
            self.rand_syst[name] = self.rand_syst[name][w_rand] 
        self.data_we = self.data_we[w_data]
        self.rand_we = self.rand_we[w_rand]
        self.w_data = w_data
        self.w_rand = w_rand
        self.ndata = self.data_we.size
        self.nrand = self.rand_we.size
        self.factor = np.sum(self.rand_we)/np.sum(self.data_we)
        if verbose:
            print('Number of galaxies before/after cutting outliers: ', 
                  w_data.size, np.sum(w_data))
            print('Number of randoms  before/after cutting outliers: ', 
                  w_rand.size, np.sum(w_rand))

    def prepare(self, nbins=10):

        nsyst = self.nsyst
        data_syst = self.data_syst
        rand_syst = self.rand_syst

        #-- compute histograms        
        edges      = {}
        centers    = {}
        h_rand     = {}
        h_index    = {}
        h_indexr    = {}
        edelta     = {}
        
        for name in data_syst:
            syst = data_syst[name]
            edg = np.linspace(syst.min()-1e-7, syst.max()+1e-7, nbins+1)
            cen = 0.5*(edg[:-1]+edg[1:]) 
            h_index[name] = np.floor((syst   -edg[0])/\
                                     (edg[-1]-edg[0]) * nbins).astype(int).T
            h_indexr[name] = np.floor((rand_syst[name]   -edg[0])/\
                                     (edg[-1]-edg[0]) * nbins).astype(int).T
            h_rand[name], _ = np.histogram(rand_syst[name], bins=edg, 
                                        weights=self.rand_we)
            print(h_rand[name])
            h_dat,_ = np.histogram(data_syst[name],bins=edg,weights=self.data_we)
            #h_dat = np.bincount(h_index[name], weights=self.data_we)
            #h_randn = np.bincount(h_indexr[name], weights=self.rand_we)  
            print(h_dat)
            #print(h_randn)  
            edges[name] = edg
            centers[name] = cen
            edelta[name] =  np.sqrt((h_dat   /h_rand[name]**2 + \
                                     h_dat**2/h_rand[name]**3 )) * self.factor \
                            + 1e10*(h_dat==0)

        self.edges = edges
        self.centers = centers
        self.h_index = h_index
        self.h_rand = h_rand
        self.edelta = edelta
        self.nbins = nbins
        
    def get_subsample(self, wd):
        wd = wd[self.w_data]
        s = Syst(self.data_we[wd], self.rand_we)
        for name in self.data_syst:
            s.add_syst(name, self.data_syst[name][wd], self.rand_syst[name])
        s.nbins = self.nbins
        s.factor = np.sum(s.rand_we)/np.sum(s.data_we)
        s.edges = self.edges
        s.centers = self.centers
        s.h_rand = self.h_rand
        s.h_index = {}
        s.edelta  = {}
        for name in s.data_syst:
            syst = s.data_syst[name]
            edg = s.edges[name]
            s.h_index[name] = np.floor((syst   -edg[0])/\
                                       (edg[-1]-edg[0]) * s.nbins).astype(int).T
            h_dat = np.bincount(s.h_index[name], weights=s.data_we,
                                minlength=s.nbins)
            s.edelta[name] =  np.sqrt((h_dat   /s.h_rand[name]**2 + \
                                       h_dat**2/s.h_rand[name]**3 )) * s.factor \
                            + 1e10*(h_dat==0)

        return s
        
    def get_model(self, pars, syst): 
        ''' Compute model from parameters and systematic values
            Input
            ------
            pars : dictionary containing parameters of fit
            syst : dictionary containing systematic values
        '''

        #-- same but using dictionary
        model = 1.+pars['constant']
        #print(pars['constant'])
        #print(pars)
        #model = np.ones(len(self.data_we))+pars['constant']
        for p in pars:
            if p == 'constant': continue
            #print(pars[p])
            edges = self.edges[p]
            edgemin, edgemax = edges[0], edges[-1]
            mp = pars[p]* (syst[p]-edgemin)/(edgemax-edgemin)
            #print(p,len(mp))
            model += mp
        return model

    def get_histograms(self, pars=None):
        data_syst = self.data_syst
        data_we = self.data_we
        
        h_rand = self.h_rand
        h_index = self.h_index

        h_data = {} 
        delta = {}

        if pars is None:
            we_model = data_we*0+1
        else:
            we_model = 1/self.get_model(pars, data_syst)

        #-- doing histograms with np.bincount, it's faster
        for name in data_syst:
            #h_dat = np.bincount(h_index[name], weights=data_we*we_model, 
            #                    minlength=self.nbins)
            #print(name,len(data_syst[name]),len(data_we),len(we_model))
            h_dat,_ = np.histogram(data_syst[name],bins=self.edges[name],weights=data_we*we_model)
            h_ran = h_rand[name]
            #-- computing overdensity and error assuming poisson
            delt = h_dat/h_ran * self.factor
            #edelt = np.sqrt((h_dat   /h_ran**2 + \
            #                 h_dat**2/h_ran**3 )) * self.factor
            h_data[name] = h_dat
            delta[name] = delt
            #edelta[name] = edelt

        self.h_data = h_data
        self.delta = delta


    def get_chi2(self, *pars):
        ''' Computes chi2 for a set of parameters 
            - for minuit fitter, pars is a tuple
            - but usually it is easy to give a dictionary
            - if no argument is give, compute chi2 for constant=1 and zero slopes
        '''
        #print('length of pars is '+str(len(pars)))
        #print(pars)
        if len(pars) == 0: 
           self.get_histograms()
        elif isinstance(pars[0], dict):
           self.get_histograms(pars=pars[0])
        else:
           pars_dict = {}
           for par_name, p in zip(self.par_names, list(pars)[0]):
               pars_dict[par_name] = p
           #print('pars_dict is '+str(pars_dict))
           self.get_histograms(pars=pars_dict)
        #self.get_histograms()
        chi2 = 0.
        for name in self.syst_names:
            chi2+= np.sum( (self.delta[name]-1)**2/self.edelta[name]**2)
        return chi2

    def fit_minuit(self, fit_maps=None, fixes=None, limits=None, priors=None):

        #-- If fit_maps is None, fit all maps 
        #-- Otherwise, define indices of maps to be fitted
        if fit_maps is None:
            fit_maps = self.syst_names
            #fit_index = np.arange(len(fit_maps), dtype=int)
        else:
            maps = self.syst_names
            for fit_map in fit_maps:
                if fit_map not in maps:
                    print(fit_map, 'not available for fitting')
                    fit_maps.remove(fit_map)
            #fit_index = []
            #fit_maps_ordered = []
            #for i in range(len(maps)):
            #    if maps[i] in fit_maps:
            #        fit_index.append(i)
            #        fit_maps_ordered.append(maps[i])
            #fit_index = np.array(fit_index, dtype=int)
            #fit_maps = np.array(fit_maps_ordered)
        #self.fit_index = fit_index
        self.fit_maps = fit_maps

#         par_names = []
        init_pars = {}
        init_errs ={}
#         par_names.append('constant')
        init_pars['constant'] = 0.
        init_errs['error_constant'] = 0.1
         
        for par in self.fit_maps:
            value = 0
            init_pars[par] = value
            init_errs['error_'+par] = abs(value)/10. if value!=0 else 0.1
#             par_names.append(par)
        par_names = [par for par in init_pars]
        #print(par_names)
        pars_values = [ init_pars[par] for par in init_pars]
        print(tuple(pars_values))
        mig = Minuit(self.get_chi2, tuple(pars_values), name=tuple(par_names))
        mig.errordef = Minuit.LEAST_SQUARES
# 
        self.fixes = fixes
        #print(init_errs)
        for par in init_pars:
            mig.errors[par] = init_errs['error_'+par]
            if fixes:
                mig.fixed[par] = fixes[par] if par in fixes else False
            if limits:
                mig.limits[par] = limits[par] if par in limits else (None, None)


#          if fixes:
#              for key in fixes.keys():
#                  init_pars[key] = fixes[key]
#                  init_pars['fix_'+key] = True 
#          if limits:
#              for key in limits.keys():
#                  init_pars['limit_'+key] = (limits[key][0], limits[key][1])
        

        self.priors = priors
        self.par_names = par_names
        

        print('Maps available for chi2:')
        print(self.syst_names)
        print('Fitting for:')
        print(self.par_names)

        #mig = Minuit(self.get_chi2,  \
#         mig = Minuit(self.get_chi2, throw_nan=False, \
#                              forced_parameters=par_names, \
#                              print_level=1, errordef=1, \
#                              **init_pars)
                             #frontend=iminuit.frontends.ConsoleFrontend(), \

        mig.tol = 1.0 
        imin = mig.migrad()
        self.mig = mig
        self.imin = imin
        #self.is_valid = imin[0]['is_valid']
        self.best_pars = mig.values 
        self.errors = mig.errors
        self.chi2min = mig.fval
        self.ndata = self.nbins*self.nsyst
        #self.npars = mig.narg
        self.npars = len(par_names)
        self.covariance = mig.covariance
        #for par in par_names:
        #    if mig.fitarg['fix_'+par]:
        #        self.npars -= 1
        self.rchi2min = self.chi2min/(self.ndata-self.npars)
        self.chi2_before = self.get_chi2()
        self.rchi2_before =  self.get_chi2()/self.ndata
        print('chi2 (before fit) = %.2f   ndata = %d                rchi2 = %.4f'%\
                (self.chi2_before, self.ndata, self.rchi2_before))
        print('chi2 (after  fit) = %.2f   ndata = %d   npars = %d   rchi2 = %.4f'%\
                (self.chi2min, self.ndata, self.npars, self.rchi2min))

    def plot_overdensity(self, pars=[None], ylim=[0.75, 1.25], 
        nbinsh=50, title=None):

        #-- setting up the windows
        nmaps = self.nsyst
        figsize = (15, 3) if nmaps > 1 else (5,3)
        f, ax = plt.subplots(1, nmaps, sharey=True, figsize=figsize)
        if nmaps == 1:
            ax = [ax] 
        if nmaps > 1: 
            f.subplots_adjust(wspace=0.05, left=0.05, right=0.98, 
                              top=0.98, bottom=0.15)
        ax[0].set_ylim(ylim)

        #-- compute histograms for before/after parameters
        nbins = self.nbins
        centers = self.centers

        for par in pars:
            self.get_histograms(pars=par)
            for i in range(nmaps):
                name   = self.syst_names[i]
                delta  = self.delta[name]
                edelta = self.edelta[name]
                chi2   = np.sum( (delta-1.)**2/edelta**2)
                label = r'$\chi^2_{r}  = %.1f/%d = %.2f$'%\
                          (chi2, nbins, chi2/nbins)
                ax[i].errorbar(centers[name], delta, edelta, fmt='.', label=label)
                ax[i].axhline( 1.0, color='k', ls='--')
                ax[i].locator_params(axis='x', nbins=5, tight=True)
                
                #-- add title and legend
                ax[i].legend(loc=0, numpoints=1, fontsize=8)
                ax[i].set_xlabel(name)

        #-- overplot histogram (normalizing to the 1/3 of the y-axis)
        for i in range(nmaps):
            name = self.syst_names[i]
            h_syst, bins = np.histogram(self.data_syst[name], bins=nbinsh)
            x = 0.5*(bins[:-1]+bins[1:])
            y = h_syst/h_syst.max()*0.3*(ylim[1]-ylim[0])+ylim[0]
            ax[i].step(x, y, where='mid', color='g')

        ax[0].set_ylabel('Density fluctuations')
        if title:
            f.subplots_adjust(top=0.9)
            plt.suptitle(title)

    def export(self, fout):

        fout = open(fout, 'w')
        nbins = self.centers[0].size
       
        delta_before, edelta_before = self.get_histograms()
        delta_after,  edelta_after  = self.get_histograms(pars=self.pars)
        chi2_before = np.sum( (delta_before-1)**2/edelta_before**2)
        chi2_after  = np.sum( (delta_after -1)**2/edelta_after**2 )
        rchi2_before = chi2_before/(self.ndata)
        rchi2_after  = chi2_after /(self.ndata-self.npars)
        
        print('#- Photometric systematic fit output by Julian Bautista', file=fout) 
        print('#- Number of systematic maps read: %d'%self.nsyst, file=fout)
        print('#- Number of systematic maps fitted: %d'%len(self.fit_maps), file=fout)
        print('#- Number of bins per systematic: %d'%nbins, file=fout)
        print('#- Maps read:', file=fout)
        print('#- '+' '.join(self.syst_names), file=fout)
        print('#- Maps fitted:', file=fout)
        print('#- '+' '.join(self.fit_maps), file=fout)
        print('#- chi2 ndata npars rchi2 (before fit)', file=fout)
        print('#- %f %d %d %f'%(chi2_before, self.ndata, 0, rchi2_before), file=fout)
        print('#- chi2 ndata npars rchi2 (after fit)', file=fout)
        print('#- %f %d %d %f'%(chi2_after, self.ndata, self.npars, rchi2_after), file=fout)
        print('#-- Parameters (const + slope per fitted systematic)', file=fout)
        print('#--', self.pars, file=fout) 

        
 
        for j in range(self.nsyst):
            sname = self.syst_names[j]
            line = '#- %s_min  %s_cen  %s_max'%\
                   (sname, sname, sname)
            line += '  delta_before  edelta_before  delta_after  edelta_after \t'
            print(line, file=fout)

            for i in range(nbins):
                smin = self.edges[j, i]
                smax = self.edges[j, i+1]
                scen = self.centers[j, i]
                den0 =   delta_before[j, i]
                eden0 = edelta_before[j, i]
                den1 =   delta_after[j, i]
                eden1 = edelta_after[j, i]
                line = '%f \t %f \t %f \t %f \t %f \t %f \t %f'%\
                       (smin, scen, smax, den0, eden0, den1, eden1)
                print(line, file=fout)

        fout.close()



def get_pix(nside, ra, dec, nest=0):
    return hp.ang2pix(nside, np.radians(-dec+90), np.radians(ra), nest=nest)

def flux_to_mag(flux, band, ebv=None):
    ''' Converts SDSS fluxes to magnitudes, correcting for extinction optionally (EBV)'''
    #-- coefs to convert from flux to magnitudes
    b = np.array([1.4, 0.9, 1.2, 1.8, 7.4])[band]*1e-10
    mag = -2.5/np.log(10)*(np.arcsinh((flux/1e9)/(2*b)) + np.log(b))
    #-- extinction coefficients for SDSS u, g, r, i, and z bands
    ext_coeff = np.array([4.239, 3.303, 2.285, 1.698, 1.263])[band]
    if not ebv is None:
        mag -= ext_coeff*ebv
    return mag



def fit_slopes_per_xbin(s, x_name, x_data, x_nbins=10, p=1., fit_maps=None, plot_delta=False,
                save_plot_delta=False, sample_name=None, sample_root=None):
    ''' From a systematic_fitter.Syst object, divide in subsamples based 
        on the value of a variable called zname with values data_z 
    '''

    #-- define x_bins excluding p% of extreme values
    x_bins = np.array([ np.percentile(x_data, p/2+(i*(100-p)/(x_nbins))) 
                       for i in range(x_nbins+1)])
    xcen = 0.5*(x_bins[1:]+x_bins[:-1]) 

    #-- List of best fit parameters
    chi2_list = []

    #-- List of fitter objects
    slist = []

    #-- Loop over x_bins
    for i in range(x_nbins):
        xmin = x_bins[i]
        xmax = x_bins[i+1]
        wx = (x_data>=xmin) & (x_data<xmax)

        print('===')
        #print(f'=== Getting subsample {i+1} of {x_nbins}, {xmin:.3f} < {x_name} < {xmax:.3f} ===')
        print('===')
        
        ss = s.get_subsample(wx)
        ss.fit_minuit(fit_maps=fit_maps)
        if plot_delta:
            ss.plot_overdensity(pars=[None, s.best_pars, ss.best_pars], ylim=[0., 2.])#, 
                 #title=f'{sample_name}: {xmin:.3f} < {x_name} < {xmax:.3f} fit in this bin #{i}')
            if save_plot_delta:
                plt.savefig('plots/syst_'+str(sample_root)+'_bin'+str(i)+'.pdf')

        chi2_list.append({'before':ss.get_chi2(), 
                     'global':ss.get_chi2(s.best_pars), 
                     'bin':ss.get_chi2(ss.best_pars)})
        slist.append(ss)

    return x_bins, chi2_list, slist

def fit_smooth_slopes_vs_x(x_bins, s_list, chi2_thres=16):
    
    xmin = x_bins[0]
    xmax = x_bins[-1]
    xcen = 0.5*(x_bins[:-1]+x_bins[1:])

    par_names = s_list[0].par_names
    npar = len(par_names)
    
    coeffs = {}
    
    print('\nFitting polynomials over slopes per bin')
    for par_name in par_names:
        y  = np.array([ss.best_pars[par_name] for ss in s_list])
        dy = np.array([ss.errors[par_name]    for ss in s_list])

        #-- Fitting slopes with polynomials with increasing order
        for order in range(4):
            #-- Fit poly
            coeff = np.polyfit(xcen, y, order, w=1/dy)
            #-- Compute model to get chi2
            ymodel = np.polyval(coeff, xcen)
            chi2 = np.sum((y-ymodel)**2/dy**2)
            if par_name in coeffs:
                #-- check if chis is smaller by at least chi2_thres units than lower order
                coeff_before = coeffs[par_name]
                ymodel_before = np.polyval(coeff_before, xcen)
                chi2_before = np.sum((y-ymodel_before)**2/dy**2)
                if (chi2 < chi2_before - chi2_thres):
                    #print('  ', par_name, ' fit with ', order, 'order poly with chi2= %.1f/%d = %.2f'%(chi2, y.size, chi2/y.size))
                    coeffs[par_name] = coeff
            else:
                coeffs[par_name] = coeff

    return coeffs


def plot_slopes_vs_x(x_bins, s_list, x_name='Z', ylim=None, title=None,
    global_pars=None, global_errors=None):
    
    xmin = x_bins[0]
    xmax = x_bins[-1]
    xcen = 0.5*(x_bins[:-1]+x_bins[1:])

    par_names = s_list[0].par_names
    npar = len(par_names)
    
    #-- some options for plotting
    figsize = (15, 3) if npar > 2 else (6,3)
    f, ax = plt.subplots(1, npar, sharey=False, figsize=figsize)
    if npar == 1:
        ax = [ax]
    if npar > 1:
        f.subplots_adjust(wspace=0.13, left=0.05, right=0.98,
                          top=0.98, bottom=0.15)
    if not ylim is None:
        ax[0].set_ylim(ylim)

    for par_name, axx in zip(par_names, ax):
        y  = np.array([ss.best_pars[par_name] for ss in s_list])
        dy = np.array([ss.errors[par_name]    for ss in s_list])
        axx.errorbar(xcen, y, dy, fmt='o', ms=3)
        #-- Fitting slopes with polynomials with increasing order
        for order in range(3):
            coeff = np.polyfit(xcen, y, order, w=1/dy)
            #-- Compute model to get chi2
            ymodel = np.polyval(coeff, xcen)
            chi2 = np.sum((y-ymodel)**2/dy**2)
            #-- Plot the model using the full range of x_data
            x = np.linspace(xmin, xmax, 30)
            ymodel = np.polyval(coeff, x)
            axx.plot(x, ymodel, label=r'$n_{\rm poly}=%d, \chi^2 = %.1f$'%(order, chi2))
        if not global_pars is None:
            x = np.median(xcen)
            y = global_pars[par_name]
            dy = global_errors[par_name]
            axx.errorbar(x, y, dy, fmt='*', ms=8, label='Global fit') 
        axx.locator_params(axis='x', nbins=4, tight=True)
        axx.legend(loc=0, fontsize=8)
        axx.set_xlabel(x_name)
        axx.set_title(par_name)

    if title:
        f.subplots_adjust(top=0.85)
        plt.suptitle(title)


def get_pars_from_coeffs(coeffs, x_values):

    pars_bin = {}
    for par_name in coeffs:
        pars_bin[par_name] = np.polyval(coeffs[par_name], x_values)
    return pars_bin

def get_chi2_xbin_smooth(s, x_bins, x_data, coeffs, chi2s,
    plot_delta=False, sample_name=None, x_name=None):
    ''' Compute chi2 for each bin, this time with the smooth poly parameters
    '''
    for i in range(x_bins.size-1):
        xmin = x_bins[i]
        xmax = x_bins[i+1]
        wdx = (x_data>=xmin) & (x_data<xmax) & (s.w_data)
        ss = s.get_subsample(wdx)
        pars_xbin = get_pars_from_coeffs(coeffs, x_data[wdx])
        chi2s[i]['bin_smooth'] = ss.get_chi2(pars_xbin)
        if plot_delta:
            ss.plot_overdensity(pars=[None, s.best_pars, pars_xbin], ylim=[0, 2])#, 
                #title=f'{sample_name}: {xmin:.3f} < {x_name} < {xmax:.3f} fit in this bin')

def plot_chi2_vs_x(x_bins, chi2s, x_name='redshift', title=None):
   
    xcen = 0.5*(x_bins[:-1]+x_bins[1:])
    plt.figure()
    for chi2name in chi2s[0]:
        c2 = np.array([chi2[chi2name] for chi2 in chi2s])
        plt.plot(xcen, c2, 'o-', label=chi2name)
    plt.legend(loc=0)
    plt.ylabel(r'$\chi^2$ per %s bin'%x_name)
    plt.xlabel(x_name)
    if title:
        plt.title(title)

def ra(x):
    return x-360*(x>300)

#-- Plot extreme values of weights in the sky
def plot_weights_sky(data_ra, data_dec, data_weightsys, title=None):
    plt.figure(figsize=(12, 7))
    plt.scatter(ra(data_ra), data_dec, c=data_weightsys, vmin=0.5, vmax=1.5, lw=0, s=2, cmap='jet', label=None)
    wext = (data_weightsys<0.5)|(data_weightsys > 2.)
    if sum(wext)>0:
        plt.plot(ra(data_ra[wext]), data_dec[wext], 'ko', ms=4, 
            label=r'$w_{\rm sys} < 0.5 \ {\rm or} \ w_{\rm sys} > 2.$ : %d galaxies'%sum(wext))
    plt.xlabel('RA [deg]')
    plt.ylabel('DEC [deg]')
    c = plt.colorbar()
    c.set_label('WEIGHT_SYSTOT')
    if title:
        plt.title(title)
    plt.legend(loc=0)
    plt.tight_layout()
