# plot the successful rate (SSR) of tracers w.r.t to observation conditions
# suses "full" catalogs as inputs
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import sys
import argparse

import fitsio
from astropy.table import Table

parser = argparse.ArgumentParser()
parser.add_argument('--basedir', help='where to find catalogs', type=str, default='/global/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--tracers", help="only ELG_LOPnotqso is available",default='ELG_LOPnotqso')
parser.add_argument("--zmin", help="the redshift lower limit",default=0.01, type=float)
parser.add_argument("--zmax", help="the redshift upper limit",default=1.6, type=float)

args = parser.parse_args()


indir = args.basedir+'/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
outdir = indir+'plots/ssr/'

# create the susscessful rate vs observation figure
if args.data == 'LSS':
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print('made '+outdir)

def SSR(fullcata, quantity, selection_quality, selection_goodz, weights, fiberbins=None):
    # the binning range of the observing condition
    if (quantity == 'TSNR2_ELG'):
        smin= 80
        smax= 200
        a,bins  = np.histogram(fullcata[quantity][selection_quality],range=(smin,smax))
    else:
        smin= np.nanpercentile(fullcata[quantity][selection_quality],(1,99))[0]
        smax= np.nanpercentile(fullcata[quantity][selection_quality],(1,99))[1]
        a,bins  = np.histogram(fullcata[quantity][selection_quality],range=(smin,smax))
    # all ELG targets, good ELG targets, binning and error
    b,_     = np.histogram(fullcata[quantity][selection_goodz],bins=bins,weights=weights)            
    BIN     = (bins[:-1]+bins[1:])/2
    err     = np.sqrt(b*(1-b/a))/a 
    return [a,b,BIN,err]

def SSR_chi2(goodz, allz, err):
    standard = np.sum(goodz)/np.sum(allz)
    return np.sum((goodz/allz-standard)**2/err**2)

# list all tracers
tps = [args.tracers]
for tp in tps:
    # only ELG-related SSR is included, need help for other tracers
    if tp.find('ELG') == -1:
        print('the sccessful rate of {} is still developping'.format(tp))
        os.exit()
    else:
        # read the full catalogue 
        full = Table(fitsio.read(indir+tp+'_full.dat.fits'))
        # add new deducted observing conditions
        full['COADD_EXPTIME-EBVFAC2'] = full['COADD_EXPTIME']/10**(2*2.165*full['EBV']/2.5)
        full['FIBERFLUX_G/MW_TRANSMISSION_G'] = full['FIBERFLUX_G']/full['MW_TRANSMISSION_G']
        full['FIBERASSIGN'] = np.sqrt(full['FIBERASSIGN_X']**2+full['FIBERASSIGN_Y']**2)
        quantities = ['TSNR2_ELG','COADD_EXPTIME-EBVFAC2','FIBERASSIGN'] #,'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z','FIBERFLUX_G/MW_TRANSMISSION_G','EBV'

        # sample selections for ELG-related samples
        ## 'N' = 'BASS+MzLS' photometric information, 'S' = 'DeCaLS'
        seln = full['PHOTSYS'] == 'N' 
        ## the selection of valid samples
        sel_obs = full['ZWARN'] != 999999
        sel_obs &= full['TSNR2_ELG'] > 80
        ## the selection of good redshift
        gz = full['o2c'] > 0.9
        print('In the {} full catalogue, there are {} objects in total, {} valid objects and {} objects with good redshift'.format(tp,len(full),len(full[sel_obs]), len(full[gz&sel_obs])))
        print('The overall successful rate of {} is {:.1f}%'.format(tp,len(full[gz&sel_obs])/len(full[sel_obs])*100))

        # plot the successful rate
        nrows     = len(quantities)//4 if len(quantities)%4 ==0 else len(quantities)//4+1
        plt.rc('font', family='serif', size=12)
        fig   = plt.figure(figsize=(4*5,nrows*5))
        spec  = gridspec.GridSpec(nrows=nrows,ncols=4,left = 0.05,right = 0.98,bottom=0.1,top = 0.98,wspace=0.2)#,hspace=0.15,wspace=0)
        ax    = np.empty((nrows,4), dtype=type(plt.axes))
        for q in range(len(quantities)):
            i,j     = q//4,q%4
            ax[i,j] = fig.add_subplot(spec[i,j])
            # define the observing condition: quantity and zmin,zmax
            quantity = quantities[q]
            zmin     = args.zmin
            zmax     = args.zmax
            # the selection of redshift
            selz    = (zmin<full['Z_not4clus'])&(full['Z_not4clus']<zmax)

            # the histogram of valid samples w.r.t quantities
            # and that of the weighted good-redshift samples w.r.t quantities in PHOTOSYS "N"
            ## the histogram range
            for split in ['N','S']:
                if split == 'N':
                    selection = sel_obs&seln
                    fmt='ro'
                    fmt_model='r'
                elif split == 'S':
                    selection = sel_obs&~seln
                    fmt='bo'
                    fmt_model='b'
                
                # select the good z at a given redshift range, selected by selq
                selection_gz = selection&gz&selz
                
                ## the histogram of valid samples w.r.t quantities
                ## and that of the weighted good-redshift samples w.r.t quantities 
                ALL, GOOD, BIN, err = SSR(full, quantity, selection, selection_gz, weights=full['WEIGHT_ZFAIL'][selection_gz])
                meanssr = np.sum(GOOD)/np.sum(ALL)
                ax[i,j].errorbar(BIN,GOOD/ALL/meanssr,err/meanssr,label=split+r': ZFAIL $\chi^2/dof={:.1f}/{}$'.format(SSR_chi2(GOOD,ALL,err),len(ALL)),fmt=fmt)
                plt.xlabel(f'{quantity} at {zmin}<z<{zmax}')

                # the SSR corrected by model
                ALL, GOOD_uncorr, BIN, err_uncorr = SSR(full, quantity, selection, selection_gz, weights=np.ones(np.sum(selection_gz)))
                meanssr_uncorr = np.sum(GOOD_uncorr)/np.sum(ALL)
                plt.fill_between(BIN,(GOOD_uncorr/ALL+err_uncorr)/meanssr_uncorr,(GOOD_uncorr/ALL-err_uncorr)/meanssr_uncorr,color=fmt_model,alpha=0.2,label='_hidden')
                ax[i,j].plot(BIN,GOOD_uncorr/ALL/meanssr_uncorr,label=split+r': unweighted $\chi^2/dof={:.1f}/{}$'.format(SSR_chi2(GOOD_uncorr,ALL,err_uncorr),len(ALL)),color=fmt_model,alpha=0.5)

            plt.grid(True)        
            plt.legend()
            plt.ylabel('{} z success rate'.format(tp))
        
        plt.savefig(outdir+'{}_success_rate_z{}z{}.png'.format(tp,args.zmin,args.zmax))        
        plt.close('all')

