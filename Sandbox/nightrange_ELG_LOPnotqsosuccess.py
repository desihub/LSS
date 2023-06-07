#right now, requires source /project/projectdirs/desi/software/desi_environment.sh master
from astropy.table import Table
import numpy as np
import os
import argparse
import fitsio
from desitarget.targetmask import zwarn_mask,desi_mask

from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser()
parser.add_argument("--min_night")#, help="use this if you want to specify the night, rather than just use the last one",default=None)
parser.add_argument("--max_night")#, help="use this if you want to specify the night, rather than just use the last one",default=None)
parser.add_argument("--plottsnr2",default='y')
parser.add_argument("--plotnz",default='y')
parser.add_argument("--vis",default='n',help="whether to display plots when you run")
parser.add_argument("--outdir",default='/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/plots/tests/')
args = parser.parse_args()

#one list for each petal for total targets
gz = np.zeros(10)
tz = np.zeros(10)
tsnrlsg = {x: [] for x in range(0,10)}
tsnrls = {x: [] for x in range(0,10)}

nzls = {x: [] for x in range(0,10)}
nzla = []

ss = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')

nights = np.unique(ss['LASTNIGHT'])
sel = nights >= int(args.min_night)
sel &= nights <= int(args.max_night)

nights = nights[sel]

bit = desi_mask['ELG_LOP']

for night in nights:# range(int(args.min_night),int(args.max_night)+1):
    month = str(night)[:6]
    #get the right tileids
    exps = Table.read('/global/cfs/cdirs/desi/spectro/redux/daily/exposure_tables/'+month+'/exposure_table_'+str(night)+'.csv')
    print('number of exposures found:')
    print(len(exps))
    #cut to dark tiles
    sel = exps['FAPRGRM']=='dark'
    print('number that are dark time:')
    print(len(exps[sel]))

    exps = exps[sel]



    #get the list of tileids observed on the last night
    tidl = np.unique(exps['TILEID'])

    #get total exposure time for tiles 
    exptl = np.zeros(len(tidl))
    for ii in range(0, len(tidl)):
        w = exps['TILEID'] == tidl[ii]
        expt = np.sum(exps[w]['EFFTIME_ETC'])
        exptl[ii] = expt


    sel = exptl > 850
    tidl = tidl[sel]

    print('number bright tiles that have EFFTIME_ETC/goal > 0.85 during the night:')
    print(len(tidl))


    print('looking at LRG redshift results from the night '+str(night))
    print('the tileids are:')
    print(tidl)



    zdir = '/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/'

    for tid in tidl:
        for pt in range(0,10):
        
            zmtlff = zdir+str(tid)+'/'+str(night)+'/zmtl-'+str(pt)+'-'+str(tid)+'-thru'+str(night)+'.fits'
            rrf = zdir+str(tid)+'/'+str(night)+'/redrock-'+str(pt)+'-'+str(tid)+'-thru'+str(night)+'.fits'
            emf = zdir+str(tid)+'/'+str(night)+'/emline-'+str(pt)+'-'+str(tid)+'-thru'+str(night)+'.fits'

            if os.path.isfile(zmtlff):
                zmtlf = fitsio.read(zmtlff)
                rr = fitsio.read(rrf,ext='TSNR2')
                em = fitsio.read(emf)
                nodata = zmtlf["ZWARN"] & zwarn_mask["NODATA"] != 0
                num_nod = np.sum(nodata)
                print('looking at petal '+str(pt)+' on tile '+str(tid))
                print('number with no data '+str(num_nod))
                badqa = zmtlf["ZWARN"] & zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA") != 0
                num_badqa = np.sum(badqa)
                print('number with bad qa '+str(num_badqa))
                nomtl = nodata | badqa
                wfqa = ~nomtl
                wlrg = ((zmtlf['DESI_TARGET'] & bit) > 0)
                wlrg &= ((zmtlf['DESI_TARGET'] & 4) == 0)
                zlrg = zmtlf[wfqa&wlrg]
                if len(zlrg) > 0:
                    #drz = (10**(3 - 3.5*zmtlf['Z']))
                    #mask_bad = (drz>30) & (zmtlf['DELTACHI2']<30)
                    #mask_bad |= (drz<30) & (zmtlf['DELTACHI2']<drz)
                    #mask_bad |= (zmtlf['DELTACHI2']<10)
                    #wz = zmtlf['ZWARN'] == 0
                    #wz &= zmtlf['Z']<1.4
                    #wz &= (~mask_bad)
                    #mask_bad = (zmtlf['DELTACHI2']<15)
                    #wz = zmtlf['ZWARN'] == 0
                    #wz &= zmtlf['Z']<1.5
                    #wz &= (~mask_bad)
                    o2c = np.log10(em['OII_FLUX'] * np.sqrt(em['OII_FLUX_IVAR']))+0.2*np.log10(zmtlf['DELTACHI2'])
                    wz = o2c > 0.9
                    wzwarn = wz#zmtlf['ZWARN'] == 0
                    gzlrg = zmtlf[wzwarn&wlrg&wfqa]
                    print('The fraction of good ELG_LOPnotqso is '+str(len(gzlrg)/len(zlrg))+' for '+str(len(zlrg))+' considered spectra')
                    gz[pt] += len(gzlrg)
                    tz[pt] += len(zlrg)
                    nzls[pt].append(zmtlf[wzwarn&wlrg]['Z'])
                    nzla.append(zmtlf[wzwarn&wlrg]['Z'])
                    tsnrlsg[pt].append(rr[wzwarn&wlrg&wfqa]['TSNR2_ELG'])
                    tsnrls[pt].append(rr[wfqa&wlrg]['TSNR2_ELG'])

                else:
                    print('no good elg data')  
            else:
                print(zmtlff+' not found') 
        

print('the total number of ELG_LOPnotqso considered per petal for the nights is:')
print(tz)
tzs = gz/tz
print('the total fraction of good ELG_LOPnotqso z per petal for the nights is:')
print(tzs)

if args.plotnz == 'y':
    from matplotlib import pyplot as plt
    figs = []
    all = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/ELG_LOPnotqso_full.dat.fits')
    sel = all['o2c'] > 0.9
    sel &= all['ZWARN'] != 999999
    #sel &= all['Z_not4clus'] < 1.5
    all = all[sel]
    nza = np.concatenate(nzla)
    for pt in range(0,10):
        #plt.clf()
        if len(nzls[pt]) > 0:
            fig = plt.figure()
            nzp = np.concatenate(nzls[pt])
            a = plt.hist(nzp,range=(0.1,1.6),bins=75,density=True,label='petal '+str(pt),histtype='step')
            plt.hist(nza,bins=a[1],density=True,histtype='step',label='all petals for selected nights')
            plt.hist(all['Z_not4clus'],bins=a[1],density=True,histtype='step',label='all archived in daily',color='k')
            plt.title('ELG_LOPnotqso for nights '+args.min_night+' through '+args.max_night)
            plt.xlabel('Z')
            plt.legend(loc='lower center')
            figs.append(fig)
            #plt.savefig(args.outdir+'LRG'+args.min_night+args.max_night+'_'+str(pt)+'.png')
            #if args.vis == 'y':
            #    plt.show()
    with PdfPages(args.outdir+'ELG_LOPnotqso'+args.min_night+args.max_night+'_nzpetal.pdf') as pdf:
        for fig in figs:
            pdf.savefig(fig)
            plt.close()

if args.plottsnr2 == 'y':
    from matplotlib import pyplot as plt
    
    for pt in range(0,10):
        if len(tsnrlsg[pt]) > 0:
            
            gz = np.concatenate(tsnrlsg[pt])
            az = np.concatenate(tsnrls[pt])
            a = np.histogram(gz)
            b = np.histogram(az,bins=a[1])
            bc = a[1][:-1]+(a[1][1]-a[1][0])/2.
            err = np.sqrt(b[0]-a[0])/b[0]
            plt.errorbar(bc,a[0]/b[0],err,label='petal '+str(pt))
    plt.legend()
    plt.xlabel('TSNR2_ELG')
    plt.ylabel('redshift success rate')
    plt.savefig(args.outdir+'ELG_LOPnotqso'+args.min_night+args.max_night+'_vstsnr2.png')
    if args.vis == 'y':
        plt.show()


    

