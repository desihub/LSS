#right now, requires source /project/projectdirs/desi/software/desi_environment.sh master
from astropy.table import Table
import numpy as np
import os
import argparse
import fitsio
from desitarget.targetmask import zwarn_mask

parser = argparse.ArgumentParser()
parser.add_argument("--night", help="use this if you want to specify the night, rather than just use the last one",default=None)
parser.add_argument("--tileid", help="tileid",default=None)
parser.add_argument("--plotnz",default='n')
args = parser.parse_args()





#one list for each petal for total targets
gz = np.zeros(10)
tz = np.zeros(10)

zdir = '/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/'

nzls = {x: [] for x in range(0,10)}
nzla = []
if args.tileid is not None:
    tid = int(args.tileid)

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
		wlrg = (zmtlf['DESI_TARGET'] & 1) > 0
		zlrg = zmtlf[wfqa&wlrg]
		if len(zlrg) > 0:
			#drz = (10**(3 - 3.5*zmtlf['Z']))
			#mask_bad = (drz>30) & (zmtlf['DELTACHI2']<30)
			#mask_bad |= (drz<30) & (zmtlf['DELTACHI2']<drz)
			mask_bad = (zmtlf['DELTACHI2']<15)
			wz = zmtlf['ZWARN'] == 0
			wz &= zmtlf['Z']<1.5
			wz &= (~mask_bad)

			wzwarn = wz#zmtlf['ZWARN'] == 0
			gzlrg = zmtlf[wzwarn&wlrg]
			print('The fraction of good LRGs is '+str(len(gzlrg)/len(zlrg))+' for '+str(len(zlrg))+' considered spectra')
			gz[pt] += len(gzlrg)
			tz[pt] += len(zlrg)
			nzls[pt].append(zmtlf[wzwarn&wlrg]['Z'])
			nzla.append(zmtlf[wzwarn&wlrg]['Z'])
		else:
			print('no good lrg data')  
	else:
		print(zmtlff+' not found') 
        

print('the total number of LRG considered per petal for the night is:')
print(tz)
tzs = gz/tz
print('the total fraction of good LRG z per petal for the night is:')
print(tzs)

if args.plotnz == 'y':
    from matplotlib import pyplot as plt
    nza = np.concatenate(nzla)
    for pt in range(0,10):
        nzp = np.concatenate(nzls[pt])
        a = plt.hist(nzp,range=(0.01,1.4),bins=28,density=True)
        plt.hist(nza,bins=a[1],density=True,histtype='step')
        plt.title('petal '+str(pt))
        plt.xlabel('Z')
        plt.show()