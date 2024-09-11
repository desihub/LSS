# plot the successful rate (SSR) of tracers w.r.t to observation conditions
# suses "full" catalogs as inputs
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os
import sys
import argparse

import fitsio
from astropy.table import Table,unique,join
import LSS.common_tools as common
from LSS.globals import main


parser = argparse.ArgumentParser()
parser.add_argument('--basedir', help='where to find catalogs', type=str, default='/global/cfs/cdirs/desi/survey/catalogs/')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='Y1')
parser.add_argument("--data",help="LSS or mock directory",default='LSS')
parser.add_argument("--verspec",help="version for redshifts",default='iron')
parser.add_argument("--version", help="catalog version",default='test')
parser.add_argument("--tracers", help="only ELG_LOPnotqso is available",default='all')
parser.add_argument("--zmin", help="minimum redshift",type=float,default=None)
parser.add_argument("--zmax", help="maximum redshift",type=float,default=None)
parser.add_argument("--focalplane_SSR_LSS", help="add WEIGHT_focal to the full data or not",action='store_true',default=False)
parser.add_argument("--fullonly", help="use full data instead of full_HPmapcut",action='store_true',default=False)


args = parser.parse_args()

def apply_imaging_veto(ff,reccircmasks,ebits):
    if reccircmasks is not None:
        for maskfn in reccircmasks:
            mask = common.maskcircandrec(ff,maskfn)
            ff = ff[~mask]

    if ebits is not None:
        print('number before imaging mask '+str(len(ff)))
        if ebits == 'lrg_mask':
            sel = ff['lrg_mask'] == 0
            ff = ff[sel]
        else:
            ff = common.cutphotmask(ff,ebits)
        print('number after imaging mask '+str(len(ff)))
    return ff

if args.data == 'mock':
    indir = args.basedir
    if args.basedir[-1] == '/':
        mock_number = int(args.basedir.split('/')[-2].split('mock')[1])
    else:
        mock_number =  int(args.basedir.split('/')[-1].split('mock')[1])

    outdir = 'mock%i/' % mock_number
    os.system('mkdir %s' % outdir)
else:
    indir = args.basedir+'/'+args.survey+'/'+args.data+'/'+args.verspec+'/LSScats/'+args.version+'/'
    if args.fullonly:
        outdir = indir+'plots/ssr/'
    else:
        outdir = indir+'plots/ssr_HPmapcut/'

# create the susscessful rate vs observation figure
if args.data == 'LSS':
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        print('made '+outdir)

def SSR(fullcata, quantity, selection_quality, selection_goodz, weights, fiberbins=None, binsbins=None):
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
        if np.any(binsbins) == None:
            a,bins  = np.histogram(fullcata[quantity][selection_quality],range=(smin,smax))
        else:
            a,bins  = np.histogram(fullcata[quantity][selection_quality],range=(smin,smax),bins=binsbins)
    # all  targets, good  targets, binning and error
    if np.any(binsbins) == None:
        b,_     = np.histogram(fullcata[quantity][selection_goodz],bins=bins,weights=weights)            
    else:
        b,_     = np.histogram(fullcata[quantity][selection_goodz],bins=binsbins,weights=weights)            

    BIN     = (bins[:-1]+bins[1:])/2
    err     = np.sqrt(b*(1-b/a))/a 
    return [a,b,BIN,err,bins]

def SSR_chi2(goodz, allz, err):
    standard = np.sum(goodz)/np.sum(allz)
    return np.sum((goodz/allz-standard)**2/err**2)

# list all tracers
tps = [args.tracers]
if args.tracers == 'all':
    tps = ['BGS_BRIGHT','ELG_LOPnotqso','QSO','LRG']

if args.survey == 'SV3' and args.tracers == 'all':
    tps = ['QSO','LRG','BGS_ANY','BGS_BRIGHT','ELG','ELG_HIP','ELG_HIPnotqso','ELGnotqso']
    zdw = ''
    if args.data != 'LSS':
        tps = ['QSO','LRG','ELG']

for tp in tps:
    # redshift range of different tracers
    if (args.zmin is not None)&(args.zmax is not None):
        zmin = args.zmin
        zmax = args.zmax
    else:
        if tp[:3] == 'QSO':
            zmin = 0.8
            zmax = 2.1#3.5
        elif tp[:3] == 'ELG':
            zmin = 0.8#0.01
            zmax = 1.6
        elif tp[:3] == 'LRG':
            zmin = 0.4#0.01
            zmax = 1.1#0.4
        elif tp[:3] == 'BGS':
            zmin = 0.1#0.01
            zmax = 0.4#0.5
    # SSR variation amplitude
    if tp[:3] == 'QSO':
        dv   = 0.08
    elif tp[:3] == 'ELG':
        flux = 'G'
        dv   = 0.05
    elif tp[:3] == 'LRG':
        dv   = 0.02
    elif tp[:3] == 'BGS':
        flux = 'Z'
        dv   = 0.02
    # read the full catalogue 
    if args.data == 'LSS':
        if args.fullonly:
            full = Table(fitsio.read(indir+tp+'_full.dat.fits'))
        else:
            full = Table(fitsio.read(indir+tp+'_full_HPmapcut.dat.fits'))
    elif args.data == 'mock':
        full = Table(fitsio.read(indir+'ffa_full_' + tp+'.fits'))
    # add new deducted observing conditions
    #print('full columns',full.columns)
    #print('orig len of full',len(full))
    #full = unique(full,keys=['TARGETID'])
    #print('len of full after cutting non-unique',len(full))
    if args.data == 'mock':
        if tp[:3] == 'LRG':
            lrgmask = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/'+'forFA'+str(mock_number)+'_matched_input_full_lrg_imask.fits')
            full = join(full,lrgmask,keys=['TARGETID'])
            #print(len(mock_data_tr))
        mainp = main(tp[:3],'iron','Y1')
        ebits = mainp.ebits
        reccircmasks = mainp.reccircmasks
        full = apply_imaging_veto(full,reccircmasks,ebits)

    if args.data == 'LSS':
        full["EBV"].name = "FIBER_EBV"
    # add per-tile ebv
    e = Table.read("/global/cfs/cdirs/desi/spectro/redux/daily/exposures-daily.csv")
    full["TILEID_EBV"] = np.nan
    for tileid in np.unique(full["TILEID"]):
        sel = full["TILEID"] == tileid
        full["TILEID_EBV"][sel] = e["EBV"][e["TILEID"] == tileid][0]
    if args.data == 'LSS':
        full['COADD_EXPTIME-EBVFAC2'] = full['COADD_EXPTIME']/10**(2*2.165*full['TILEID_EBV']/2.5)
        full['FIBERASSIGN'] = np.sqrt(full['FIBERASSIGN_X']**2+full['FIBERASSIGN_Y']**2)
        quantities = [f'TSNR2_{tp[:3]}','COADD_EXPTIME-EBVFAC2','FIBERASSIGN']
    elif args.data == 'mock':
        quantities = [f'TSNR2_{tp[:3]}']
    # sample selections 
    ## 'N' = 'BASS+MzLS' photometric information, 'S' = 'DECam'
    seln = full['PHOTSYS'] == 'N' 
    ## the selection of valid samples
    sel_obs = full['ZWARN'] != 999999
    sel_obs&= full['ZWARN']*0 == 0
    if args.data=='mock':
        sel_obs&=full['WEIGHT_IIP'] != 1e20
        
    if tp[:3] == 'QSO':
	    sel_obs &= full['PRIORITY'] == 3400 #repeats throw things off  
    # the selection of redshift
    #print('full columns',full.columns)
    if args.data == 'LSS':
        selz    = (zmin<full['Z_not4clus'])&(full['Z_not4clus']<zmax)
    elif args.data == 'mock':
        selz    = (zmin<full['RSDZ'])&(full['RSDZ']<zmax)
    ## the selection of good redshift
    
    '''mock_data_tr = unique(mock_data_tr,keys=['TARGETID'])
    print('length after cutting to unique targetid',len(mock_data_tr))
    selobs = mock_data_tr['WEIGHT_IIP'] != 1e20
    mock_data_tr = mock_data_tr[selobs]
    print('length after cutting to "observed" targets',len(mock_data_tr))
    mock_data_tr.rename_column('RSDZ', 'Z')
    mock_data_tr['WEIGHT_COMP'] = mock_data_tr['WEIGHT_IIP']
    if 'imaging' in args.veto:
        if tracer == 'LRG':
            lrgmask = fitsio.read(args.base_dir+'forFA'+str(args.realization)+'_matched_input_full_lrg_imask.fits')
            mock_data_tr = join(mock_data_tr,lrgmask,keys=['TARGETID'])
            print(len(mock_data_tr))
        ebits = mainp.ebits
        reccircmasks = mainp.reccircmasks
        mock_data_tr = apply_imaging_veto(mock_data_tr,reccircmasks,ebits)'''
    
    if args.data == 'LSS':
        gz = common.goodz_infull(tp[:3],full)
    elif args.data == 'mock':
        full['WEIGHT_ZFAIL'] = np.ones_like(full['RSDZ'])
        np.random.seed(mock_number)
        # Determine good redshifts by randomly downsampling
        # The random selection isn't quite right because failure rate depends on TSNR2. 
        # So the errorbars will be underestimated at low TSNR2. 
        # For LRG, the overall success rate is 99% (similar for BGS),
        #  and success rate drops about 1% at very low TSNR2. 
        # So this is about a 40% effect on the low TSNR2 errorbar 
        # (just taking the scalings from the binomial error). 
        # This inaccuracy will be worse if I also consider binning in flux and TSNR since 
        # the success rate will be even worse at low TSNR and low flux. 
        # But I'll start with the random failure rate just to have a simple model.
        if tp[:3] == 'LRG':
            success_rate = 0.99
        elif tp[:3] == 'QSO':
            success_rate = 0.665
        gz = np.zeros(len(full))
        gz[np.random.choice(np.arange(len(full)), size=int(success_rate * len(full)),replace=False)] = 1
        gz = gz.astype('bool')
        #gz = np.random.choice(np.arange(len(full)), size=int(success_rate * len(full)),replace=False)

    print('In the {} full catalogue, there are {} objects in total, {} valid objects and {} objects with good redshift'.format(tp,len(full),len(full[sel_obs]), len(full[gz&sel_obs])))
    print('The overall successful rate of {} is {:.1f}%'.format(tp,len(full[gz&sel_obs])/len(full[sel_obs])*100))

    # plot the successful rate
    ncols = 3
    nrows = len(quantities)//ncols if len(quantities)%ncols ==0 else len(quantities)//ncols+1
    plt.rc('font', family='serif', size=12)
    fig   = plt.figure(figsize=(ncols*5,nrows*5))
    spec  = gridspec.GridSpec(nrows=nrows,ncols=ncols,left = 0.05,right = 0.99,bottom=0.1,top = 0.98,wspace=0.25)#,hspace=0.15,wspace=0)
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
            weight_type = ': ZFAIL'
            if args.data == 'mock':
                bins = np.loadtxt('/global/cfs/cdirs/desi/survey/catalogs//Y1/LSS/iron/LSScats/test/plots/ssr/'+'{}_TSNR2_success_rate_z{}z{}_{}_{}_bins.txt'.format(tp,zmin,zmax,split,args.version))
                ALL, GOOD, BIN, err, bins = SSR(full, quantity, selection, selection_gz, weights=full['WEIGHT_ZFAIL'][selection_gz], binsbins=bins)
            elif args.data == 'LSS':
                if (not 'WEIGHT_focal' in full.colnames):
                    ALL, GOOD, BIN, err, bins = SSR(full, quantity, selection, selection_gz, weights=full['WEIGHT_ZFAIL'][selection_gz])
                else:
                    ALL, GOOD, BIN, err, bins = SSR(full, quantity, selection, selection_gz, weights=full['WEIGHT_ZFAIL'][selection_gz]*full['WEIGHT_focal'][selection_gz])
                    weight_type = r': ZFAIL*$\epsilon_{\rm focal}$'
            meanssr = np.sum(GOOD)/np.sum(ALL)
            ax[i,j].errorbar(BIN,GOOD/ALL/meanssr,err/meanssr,label=split+weight_type+r', $\chi^2/dof={:.1f}/{}$'.format(SSR_chi2(GOOD,ALL,err),len(ALL)),fmt=fmt)
            #print('GOOD/ALL/meanssr',GOOD/ALL/meanssr)
            plt.xlabel(f'{quantity} at {zmin}<z<{zmax}')
            if q == 0:
                np.savetxt(outdir+'{}_TSNR2_success_rate_z{}z{}_{}_{}.txt'.format(tp,zmin,zmax,split,args.version),np.array([BIN, GOOD/ALL/meanssr,err/meanssr]).T)   
                if args.data == 'LSS':
                    np.savetxt(outdir+'{}_TSNR2_success_rate_z{}z{}_{}_{}_bins.txt'.format(tp,zmin,zmax,split,args.version),bins)
                #elif args.data == 'mock':
                
            # the SSR corrected by model
            ALL, GOOD_uncorr, BIN, err_uncorr,bins = SSR(full, quantity, selection, selection_gz, weights=np.ones(np.sum(selection_gz)))
            meanssr_uncorr = np.sum(GOOD_uncorr)/np.sum(ALL)
            plt.fill_between(BIN,(GOOD_uncorr/ALL+err_uncorr)/meanssr_uncorr,(GOOD_uncorr/ALL-err_uncorr)/meanssr_uncorr,color=fmt_model,alpha=0.2,label='_hidden')
            ax[i,j].plot(BIN,GOOD_uncorr/ALL/meanssr_uncorr,label=split+r': unweighted, $\chi^2/dof={:.1f}/{}$'.format(SSR_chi2(GOOD_uncorr,ALL,err_uncorr),len(ALL)),color=fmt_model,alpha=0.5)

            ax[i,j].axhline(1,c='k')
        handles, labels = plt.gca().get_legend_handles_labels()
        order = [2,0,3,1]
        plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order],frameon=False)

        plt.grid(True)        
        plt.ylabel('{} z success rate'.format(tp))
    
    plt.savefig(outdir+'{}_success_rate_z{}z{}_{}.png'.format(tp,zmin,zmax,args.version))        
    plt.close('all')

    # obtain the fibre-wise SSR
    photos    = ['BASS/MzLS','DECam']
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
    
    # plot the SSR and chi2 on the focal plane (fibre-wise SSR)
    for ptype in ['noZFAIL','SSR','chi2','chi2hist']:
        if ptype != 'chi2hist':
            right = 0.93
        else:
            right = 0.99
        fig = plt.figure(figsize=(9,4))
        spec = gridspec.GridSpec(nrows=1,ncols=2, left = 0.1,right = right,bottom=0.12,top = 0.93, wspace=0,width_ratios=[0.85,1])
        ax = np.empty((1,2), dtype=type(plt.axes))
        plt.rc('font', family='serif', size=12)    
        for cp,split in enumerate(['N','S']):
            ax[0,cp] = fig.add_subplot(spec[0,cp])
            if split == 'N':
                selection = sel_obs&seln
            elif split == 'S':
                selection = sel_obs&~seln
            selection_gz = selection&selz&gz
            
            if ptype == 'noZFAIL':
                ALL, GOOD, BIN, err, _ = SSR(full, 'FIBER', selection, selection_gz, weights=np.ones(np.sum(selection_gz)), fiberbins=FIB)
                ptypetl = 'no ZFAIL weight'
            else:
                ALL, GOOD, BIN, err, _ = SSR(full, 'FIBER', selection, selection_gz, weights=full['WEIGHT_ZFAIL'][selection_gz], fiberbins=FIB)
                ptypetl = 'ZFAIL weight'                
            ssrmodel   = GOOD/ALL     
            ssrmean    = np.sum(GOOD)/np.sum(ALL)   
            err[err==0]= 1  
            chi2s      = (ssrmodel-ssrmean)/err

            if ptype != 'chi2hist':
                if (ptype == 'SSR')|(ptype == 'noZFAIL'):
                    # fibrewise SSR 
                    vmin    = 1-dv
                    vmax    = 1+dv
                    value   = ssrmodel/ssrmean
                    cblabel = 'rescaled SSR'
                    # mock weight, don't need to repeat with plots
                    if split == 'N':
                        ssr_wtN = 1./(ssrmodel/np.nanmean(ssrmodel))
                        ssr_wtN[np.isnan(ssr_wtN)] = 1.
                    elif split == 'S':
                        ssr_wtS = 1./(ssrmodel/np.nanmean(ssrmodel))
                        ssr_wtS[np.isnan(ssr_wtS)] = 1.
                elif ptype == 'chi2':
                    # fibrewise SSR chi2
                    dv      = 2
                    vmin    = -dv
                    vmax    = +dv
                    value   = chi2s
                    cblabel = r'SSR $\chi^2$'
                hb = ax[0,cp].scatter(xll,yll,c=value,s=2,vmin=vmin,vmax=vmax)
                if cp == 1:
                    cb = fig.colorbar(hb, ax=ax[0,1])
                    cb.set_label(cblabel,fontsize=12)
                    plt.text(-150,410,f'{photos[cp]}',fontsize=15,weight='bold')
                    plt.yticks(alpha=0)
                else:
                    plt.text(-190,410,f'{photos[cp]}',fontsize=15,weight='bold')
                    plt.ylabel('Y (mm)')
                plt.xlabel('X (mm)')
                plt.xlim(-470,470)
                plt.ylim(-420,470)
            else:
                # the histogram of fibrewise SSR chi2
                plt.hist(chi2s[np.isfinite(chi2s)],density=True,label=f'{tp} in {split}')
                plt.xlabel('chi2')
                if cp ==0:
                    plt.ylabel('normalised counts')
                plt.legend()
            plt.title('{} chi2 = {:.1f}/{}'.format(ptypetl,np.sum(chi2s[np.isfinite(chi2s)]**2),np.sum(np.isfinite(chi2s))),fontsize=10)

        if ptype == 'SSR':
            plt.savefig(outdir+'{}_focalplane_success_rate_z{}z{}_{}.png'.format(tp,zmin,zmax,args.version))        
        else:
            plt.savefig(outdir+'{}_focalplane_success_rate_{}_z{}z{}_{}.png'.format(tp,ptype,zmin,zmax,args.version))                    

        plt.close('all')

    #print('ssr_wt',ssr_wt)
    #print('BIN',list(BIN))
    #print('FIB',list(FIB))
    #print('FIBER',full['FIBER'])
    if args.data == 'LSS':
        if args.fullonly:
            full = Table(fitsio.read(indir+tp+'_full.dat.fits'))
        else:
            full = Table(fitsio.read(indir+tp+'_full_HPmapcut.dat.fits'))
    elif args.data == 'mock':
        full = Table(fitsio.read(indir+'ffa_full_' + tp+'.fits'))
 
    if (not 'WEIGHT_focal' in full.colnames)&(args.focalplane_SSR_LSS):
        full['WEIGHT_focal'] = np.ones_like(full['WEIGHT_ZFAIL'])
        for i in range(len(full['FIBER'])):
            #print('fiber',full['FIBER'][i])
            #print('fib == fiber',np.where(FIB==full['FIBER'][i]))
            #print('len ssr_wt',len(ssr_wt))
            #print('len fib',len(FIB))
            #print('len BIN',len(BIN))
            
            if full['FIBER'][i] != 999999:
                if full['PHOTSYS'][i] == 'N':
                    if len(ssr_wtN[np.where(FIB==full['FIBER'][i])]) > 0:
                    
                        print("ssr wt[FIB==full['FIBER'][i]]",ssr_wtN[np.where(FIB==full['FIBER'][i])])
                        full['WEIGHT_focal'][i] = ssr_wtN[np.where(FIB==full['FIBER'][i])]
                        print(i)
                elif full['PHOTSYS'][i] == 'S':
                    if len(ssr_wtS[np.where(FIB==full['FIBER'][i])]) > 0:
                    
                        print("ssr wt[FIB==full['FIBER'][i]]",ssr_wtS[np.where(FIB==full['FIBER'][i])])
                        full['WEIGHT_focal'][i] = ssr_wtS[np.where(FIB==full['FIBER'][i])]
                        print(i)

        full.write(indir + tp+'_full.dat.2.fits',overwrite=True)
    
