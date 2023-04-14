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
parser.add_argument("--verspec",help="version for redshifts",default='himalayas')
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--tracers", help="only ELG_LOPnotqso is available",default='ELG_LOPnotqso')
parser.add_argument("--type", help="observing-conditions or redshift-bins",default='redshift-bins')
parser.add_argument("--zmin", help="the redshift lower limit",default=0.6, type=float)
parser.add_argument("--zmax", help="the redshift upper limit",default=1.6, type=float)
parser.add_argument("--quantity", help="the observing condition: TSNR2_ELG, FIBERASSIGN, COADD_EXPTIME-EBVFAC2",default='TSNR2_ELG')

args = parser.parse_args()


indir = args.basedir+'/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
outdir = indir+'plots/ssr/'

# create the susscessful rate vs observation figure
if args.data == 'LSS':
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print('made '+outdir)

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
        if args.type == 'observing-conditions':
            quantities = ['TSNR2_ELG','COADD_EXPTIME-EBVFAC2','FIBERASSIGN','PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z','FIBERFLUX_G/MW_TRANSMISSION_G','EBV']
        elif args.type == 'redshift-bins':
            quantities = [args.zmin+0.1*binwidth for binwidth in range(int((args.zmax-args.zmin)/0.1))]

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
        nrows = len(quantities)//4+1
        plt.rc('font', family='serif', size=12)
        fig   = plt.figure(figsize=(4*6,nrows*5))
        spec  = gridspec.GridSpec(nrows=nrows,ncols=4,left = 0.05,right = 0.98,bottom=0.05,top = 0.96,wspace=0.2)#,hspace=0.15,wspace=0)
        ax    = np.empty((nrows,4), dtype=type(plt.axes))
        for q in range(len(quantities)):
            i,j     = q//4,q%4
            ax[i,j] = fig.add_subplot(spec[i,j])
            # define the observing condition: quantity and zmin,zmax
            if args.type == 'observing-conditions':
                quantity = quantities[q]
                zmin     = args.zmin
                zmax     = args.zmax
            elif args.type == 'redshift-bins':
                quantity = args.quantity
                zmin     = quantities[q]
                zmax     = quantities[q]+0.1
            # the selection of redshift
            selz    = (zmin<full['Z_not4clus'])&(full['Z_not4clus']<zmax)

            # the histogram of valid samples w.r.t quantities
            # and that of the weighted good-redshift samples w.r.t quantities in PHOTOSYS "N"
            ## the histogram range
            if quantity == 'TSNR2_ELG':
                smin= 80
                smax= 200
            else:
                smin= np.nanpercentile(full[quantity][sel_obs&selz&seln],(1,99))[0]
                smax= np.nanpercentile(full[quantity][sel_obs&selz&seln],(1,99))[1]

            ## SSR_weighted = N_{weighted good-redshift}/N_{valid samples}
            a,bins  = np.histogram(full[quantity][sel_obs&seln],range=(smin,smax))
            b,_     = np.histogram(full[quantity][sel_obs&selz&gz&seln],bins=bins,weights=full['WEIGHT_ZFAIL'][sel_obs&selz&gz&seln])
            ## binomial error of SSR_weighted
            err     = np.sqrt(b*(1-b/a))/a 
            plt.errorbar((bins[:-1]+bins[1:])/2,b/a,err,fmt='rd',label='{} N'.format(args.verspec))
            
            ## SSR_unweighted = N_{good-redshift}/N_{valid samples}
            b,_     = np.histogram(full[quantity][sel_obs&selz&gz&seln],bins=bins)
            err     = np.sqrt(b*(1-b/a))/a 
            plt.fill_between((bins[:-1]+bins[1:])/2,b/a+err,b/a-err,color='r',alpha=0.2,label='{} N unweighted'.format(args.verspec))
            plt.plot((bins[:-1]+bins[1:])/2,b/a,'r--',label='_hidden')
            
            # SSR_weighted and SSR_unweighted in PHOTOSYS "S"
            if quantity == 'TSNR2_ELG':
                smin= 80
                smax= 200
            else:
                smin= np.nanpercentile(full[quantity][sel_obs&selz&~seln],(1,99))[0]
                smax= np.nanpercentile(full[quantity][sel_obs&selz&~seln],(1,99))[1]
            a,bins  = np.histogram(full[quantity][sel_obs&~seln],range=(smin,smax))
            b,_     = np.histogram(full[quantity][sel_obs&selz&gz&~seln],bins=bins,weights=full['WEIGHT_ZFAIL'][sel_obs&selz&gz&~seln])            
            err     = np.sqrt(b*(1-b/a))/a 
            plt.errorbar((bins[:-1]+bins[1:])/2,b/a,err,fmt='bo',label='{} S'.format(args.verspec))
            
            b,_     = np.histogram(full[quantity][sel_obs&selz&gz&~seln],bins=bins)
            err     = np.sqrt(b*(1-b/a))/a 
            plt.fill_between((bins[:-1]+bins[1:])/2,b/a+err,b/a-err,color='b',alpha=0.2,label='{} S unweighted'.format(args.verspec))
            plt.plot((bins[:-1]+bins[1:])/2,b/a,'b--',label='_hidden')

            plt.legend()
            plt.ylabel('{} z success rate'.format(tp))
            plt.xlabel(quantity+' at {:.1f}<z<{:.1f}'.format(zmin,zmax))
            plt.grid()
        
        if args.type == 'observing-conditions':
            plt.savefig(outdir+'{}_success_rate_{}_z{}z{}.png'.format(tp,args.type,args.zmin,args.zmax))        
        elif args.type == 'redshift-bins':
            plt.savefig(outdir+'{}_success_rate_{}_{}.png'.format(tp,quantity,args.type))        
        plt.close('all')

