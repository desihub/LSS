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
import LSS.common_tools as common


parser = argparse.ArgumentParser()
parser.add_argument('--basedir', help='where to find catalogs', type=str, default='/global/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--tracers", help="only ELG_LOPnotqso is available",default='all')

args = parser.parse_args()


indir = args.basedir+'/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
############
outdir = indir+'plots/ssr/'
############

# create the susscessful rate vs observation figure
if args.data == 'LSS':
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print('made '+outdir)

def SSR(fullcata, quantity, selection_quality, selection_goodz, weights, fiberbins=None):
    # the binning range of the observing condition
    if (quantity == 'FIBER'):
        a,bins  = np.histogram(fullcata[quantity][selection_quality],bins=fiberbins)
    else:
        if   (quantity == 'TSNR2_ELG'):
            smin= 80
            smax= 200
        elif (quantity == 'TSNR2_BGS'):
            smin= 880
            smax= 2200
        else:
            smin= np.nanpercentile(fullcata[quantity][selection_quality],(1,99))[0]
            smax= np.nanpercentile(fullcata[quantity][selection_quality],(1,99))[1]
        a,bins  = np.histogram(fullcata[quantity][selection_quality],range=(smin,smax))
    # all  targets, good  targets, binning and error
    b,_     = np.histogram(fullcata[quantity][selection_goodz],bins=bins,weights=weights)            
    BIN     = (bins[:-1]+bins[1:])/2
    err     = np.sqrt(b*(1-b/a))/a 
    return [a,b,BIN,err]

def SSR_chi2(goodz, allz, err):
    standard = np.sum(goodz)/np.sum(allz)
    return np.sum((goodz/allz-standard)**2/err**2)

# list all tracers
tps = [args.tracers]
if args.tracers == 'all':
    tps = ['BGS_BRIGHT']#,'ELG_LOPnotqso','QSO','LRG']

if args.survey == 'SV3' and args.tracers == 'all':
    tps = ['QSO','LRG','BGS_ANY','BGS_BRIGHT','ELG','ELG_HIP','ELG_HIPnotqso','ELGnotqso']
    zdw = ''
    if args.data != 'LSS':
        tps = ['QSO','LRG','ELG']

for tp in tps:
    # redshift range of different tracers
    if tp[:3] == 'QSO':
        zmin = 0.8
        zmax = 3.5
    elif tp[:3] == 'ELG':
        zmin = 0.01
        zmax = 1.8
        flux = 'G'
    elif tp[:3] == 'LRG':
        zmin = -0.1
        zmax = 1.5
    elif tp[:3] == 'BGS':
        zmin = 0.01
        zmax = 0.5
        flux = 'Z'
    # read the full catalogue 
    full = Table(fitsio.read(indir+tp+'_full.dat.fits'))
    # add new deducted observing conditions
    full["EBV"].name = "FIBER_EBV"
    # add per-tile ebv
    e = Table.read("/global/cfs/cdirs/desi/spectro/redux/daily/exposures-daily.csv")
    full["TILEID_EBV"] = np.nan
    for tileid in np.unique(full["TILEID"]):
        sel = full["TILEID"] == tileid
        full["TILEID_EBV"][sel] = e["EBV"][e["TILEID"] == tileid][0]
    full['COADD_EXPTIME-EBVFAC2'] = full['COADD_EXPTIME']/10**(2*2.165*full['TILEID_EBV']/2.5)
    full['FIBERASSIGN'] = np.sqrt(full['FIBERASSIGN_X']**2+full['FIBERASSIGN_Y']**2)
    #full[f'FIBERFLUX_{flux}/MW_TRANSMISSION_{flux}'] = full[f'FIBERFLUX_{flux}']/full[f'MW_TRANSMISSION_{flux}']
    quantities = [f'TSNR2_{tp[:3]}','COADD_EXPTIME-EBVFAC2','FIBERASSIGN'] #,'PSFDEPTH_G', 'PSFDEPTH_R', 'PSFDEPTH_Z','FIBERFLUX_G/MW_TRANSMISSION_G','EBV'

    # sample selections 
    ## 'N' = 'BASS+MzLS' photometric information, 'S' = 'DeCaLS'
    seln = full['PHOTSYS'] == 'N' 
    ## the selection of valid samples
    sel_obs = full['ZWARN'] != 999999
    sel_obs&= full['ZWARN']*0 == 0
    if tp[:3] == 'QSO':
	    sel_obs &= full['PRIORITY'] == 3400 #repeats throw things off  
    # the selection of redshift
    selz    = (zmin<full['Z_not4clus'])&(full['Z_not4clus']<zmax)
    ## the selection of good redshift
    gz = common.goodz_infull(tp[:3],full)

    print('In the {} full catalogue, there are {} objects in total, {} valid objects and {} objects with good redshift'.format(tp,len(full),len(full[sel_obs]), len(full[gz&sel_obs])))
    print('The overall successful rate of {} is {:.1f}%'.format(tp,len(full[gz&sel_obs])/len(full[sel_obs])*100))

    # plot the successful rate
    ncols = 3
    nrows = len(quantities)//ncols if len(quantities)%ncols ==0 else len(quantities)//ncols+1
    plt.rc('font', family='serif', size=12)
    fig   = plt.figure(figsize=(ncols*5,nrows*5))
    spec  = gridspec.GridSpec(nrows=nrows,ncols=ncols,left = 0.05,right = 0.98,bottom=0.1,top = 0.98,wspace=0.2)#,hspace=0.15,wspace=0)
    ax    = np.empty((nrows,ncols), dtype=type(plt.axes))
    for q in range(len(quantities)):
        i,j     = q//ncols,q%ncols
        ax[i,j] = fig.add_subplot(spec[i,j])
        quantity = quantities[q]

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

            ax[i,j].axhline(1,c='k')

        plt.grid(True)        
        plt.legend()
        plt.ylabel('{} z success rate'.format(tp))
    
    plt.savefig(outdir+'{}_success_rate_z{}z{}_{}.png'.format(tp,zmin,zmax,args.version))        
    plt.close('all')

    # obtain the fibre-wise SSR
    dv        = 0.05
    photos    = ['BASS/MzLS','DECaLS']
    # the fibreIDs
    dl        = np.loadtxt(f'/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v0.1/{tp}_zsuccess_fromfull.txt').transpose()
    FIB       = np.append(dl[0],999999)
    # the fiber label (ID) and positions
    fibloc    = Table.read(os.environ['DESIMODEL']+'/data/focalplane/desi-focalplane_2021-06-27T20:00:34+00:00.ecsv')
    fibx_dict = dict(zip(fibloc['FIBER'], fibloc['OFFSET_X']))
    fiby_dict = dict(zip(fibloc['FIBER'], fibloc['OFFSET_Y']))
    # the throughput map
    dflux5    = np.loadtxt('/global/cfs/cdirs/desi/survey/fluxcalibration/flux-ratio-vs-xy-4500A-5500A.csv',delimiter=',').transpose()
    dflux6    = np.loadtxt('/global/cfs/cdirs/desi/survey/fluxcalibration/flux-ratio-vs-xy-6000A-7300A.csv',delimiter=',').transpose()
    dflux8    = np.loadtxt('/global/cfs/cdirs/desi/survey/fluxcalibration/flux-ratio-vs-xy-8500A-9800A.csv',delimiter=',').transpose()
    # grid the throughput map to the mesh of fiber positions
    sp        = dflux5[1][1]-dflux5[1][0]
    minv      = min(dflux5[0])-sp/2
    indx      = 49*((fibloc['OFFSET_X']-minv)//sp)
    indy  = (fibloc['OFFSET_Y']-minv)//sp
    index = (indx+indy).astype(int)
    f5val = dflux5[2][index]
    f6val = dflux6[2][index]
    f8val = dflux8[2][index]
    # the gridded throughput map mated with fiber positions
    fibf5_dict = dict(zip(fibloc['FIBER'], f5val)) 
    fibf6_dict = dict(zip(fibloc['FIBER'], f6val))
    fibf8_dict = dict(zip(fibloc['FIBER'], f8val))
    # matching the throughput value with the 
    xll = []
    yll = []
    f5l = []
    f6l = []
    f8l = []
    for fib in dl[0]:
        xll.append(fibx_dict[fib])
        yll.append(fiby_dict[fib])
        f5l.append(fibf5_dict[fib])
        f6l.append(fibf6_dict[fib])
        f8l.append(fibf8_dict[fib])
    
    # plot the fibre-wise SSR
    for cp,split in enumerate(['N','S']):
        if split == 'N':
            selection = sel_obs&seln
        elif split == 'S':
            selection = sel_obs&~seln
        selection_gz = selection&selz&gz

        ALL, GOOD, BIN, err = SSR(full, 'FIBER', selection, selection_gz, weights=full['WEIGHT_ZFAIL'][selection_gz], fiberbins=FIB)
        ssrmodel = GOOD/ALL            
        # fibrewise SSR and the correction
        plt.scatter(xll,yll,c=ssrmodel/np.nanmean(ssrmodel),s=2,vmin=1-dv,vmax=1+dv)
        plt.colorbar()
        plt.title(f'fibrewise SSR on {photos[cp]}')
        plt.savefig(outdir+'{}_focalplane_success_rate_z{}z{}_{}_{}.png'.format(tp,zmin,zmax,split,args.version))        
        plt.close()
