#right now, requires source /project/projectdirs/desi/software/desi_environment.sh master
from astropy.table import Table
import numpy as np
import os
import argparse
import fitsio
from desitarget.targetmask import zwarn_mask

parser = argparse.ArgumentParser()
parser.add_argument("--night", help="use this if you want to specify the night, rather than just use the last one",default=None)
parser.add_argument("--plotnz",default='n')
parser.add_argument("--redux",default='daily')
parser.add_argument("--vis",default='n',help="whether to display plots when you run")
parser.add_argument("--outdir",default='/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/plots/tests/')
args = parser.parse_args()


month = args.night[:6]
#get the right tileids
exps = Table.read('/global/cfs/cdirs/desi/spectro/redux/daily/exposure_tables/'+month+'/exposure_table_'+args.night+'.csv')
print('number of exposures found:')
print(len(exps))
#cut to dark tiles
sel = exps['FAPRGRM']=='bright'
print('number that are bright time:')
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


#sel &= exps['EFFTIME_ETC'] > 850 #select only tiles that should be near completion
sel = exptl/60 > .85
tidl = tidl[sel]

print('number bright tiles that have EFFTIME_ETC/goal > 0.85 during the night:')
print(len(tidl))


print('looking at LRG redshift results from the night '+str(args.night))
print('the tileids are:')
print(tidl)


#one list for each petal for total targets
gz = np.zeros(10)
tz = np.zeros(10)

zdir = '/global/cfs/cdirs/desi/spectro/redux/'+args.redux+'/tiles/cumulative/'

nzls = {x: [] for x in range(0,10)}
nzla = []
tilegzl = []
for tid in tidl:
    tile_tot = 0
    tile_good = 0
    for pt in range(0,10):
        
        zmtlff = zdir+str(tid)+'/'+args.night+'/zmtl-'+str(pt)+'-'+str(tid)+'-thru'+args.night+'.fits'
        if os.path.isfile(zmtlff):
            zmtlf = fitsio.read(zmtlff)
            nodata = zmtlf["ZWARN"] & zwarn_mask["NODATA"] != 0
            num_nod = np.sum(nodata)
            print('looking at petal '+str(pt)+' on tile '+str(tid))
            print('number with no data '+str(num_nod))
            badqa = zmtlf["ZWARN"] & zwarn_mask.mask("BAD_SPECQA|BAD_PETALQA") != 0
            num_badqa = np.sum(badqa)
            print('number with bad qa '+str(num_badqa))
            nomtl = nodata | badqa
            wfqa = ~nomtl
            wlrg = (zmtlf['BGS_TARGET'] ) > 0
            zlrg = zmtlf[wfqa&wlrg]
            if len(zlrg) > 0:
                #drz = (10**(3 - 3.5*zmtlf['Z']))
                #mask_bad = (drz>30) & (zmtlf['DELTACHI2']<30)
                #mask_bad |= (drz<30) & (zmtlf['DELTACHI2']<drz)
                #mask_bad |= (zmtlf['DELTACHI2']<10)
                #wz = zmtlf['ZWARN'] == 0
                #wz &= zmtlf['Z']<1.4
                #wz &= (~mask_bad)
                mask_bad = (zmtlf['DELTACHI2']<40)
                wz = zmtlf['ZWARN'] == 0
                wz &= zmtlf['Z']<0.7
                wz &= (~mask_bad)

                wzwarn = wz#zmtlf['ZWARN'] == 0
                gzlrg = zmtlf[wzwarn&wlrg]
                tile_tot += len(zlrg)
                tile_good += len(gzlrg)
                print('The fraction of good BGS is '+str(len(gzlrg)/len(zlrg))+' for '+str(len(zlrg))+' considered spectra')
                gz[pt] += len(gzlrg)
                tz[pt] += len(zlrg)
                nzls[pt].append(zmtlf[wzwarn&wlrg]['Z'])
                nzla.append(zmtlf[wzwarn&wlrg]['Z'])
            else:
                print('no good lrg data')  
        else:
            print(zmtlff+' not found') 
        
    print('the success rate for tile '+str(tid)+' is '+str(tile_good/tile_tot))
    if tile_tot > 0:
        tilegzl.append(tile_good/tile_tot)
    else:
        tilegzl.append('no data')
print('the total number of BGS considered per petal for the night is:')
print(tz)
tzs = gz/tz
print('the total fraction of good BGS z per petal for the night is:')
print(tzs)
print('the list of tiles and their success rates is')
print(tidl)
print(tilegzl)
if args.plotnz == 'y':
    from matplotlib import pyplot as plt
    all = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/LSScats/test/BGS_ANY_full.dat.fits')
    sel = all['ZWARN'] == 0
    sel &= all['DELTACHI2'] > 40
    sel &= all['Z_not4clus'] < 0.7
    all = all[sel]
    nza = np.concatenate(nzla)
    for pt in range(0,10):
        plt.clf()
        if len(nzls[pt]) > 0:
            nzp = np.concatenate(nzls[pt])
            a = plt.hist(nzp,range=(0.01,0.7),bins=28,density=True,label=args.night+' petal '+str(pt),histtype='step')
            plt.hist(nza,bins=a[1],density=True,histtype='step',label=args.night)
            plt.hist(all['Z_not4clus'],bins=a[1],density=True,histtype='step',label='all archived in daily',color='k')
            plt.title('BGS_ANY')
            plt.xlabel('Z')
            plt.legend(loc='upper right')
            plt.savefig(args.outdir+'BGS_ANY'+args.night+'_'+str(pt)+'.png')
            if args.vis == 'y':
                plt.show()
