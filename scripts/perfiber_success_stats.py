import numpy as np
#!pip install astropy
#!pip install fitsio
from scipy import stats
from scipy.stats import norm
import fitsio
import glob
import os
import matplotlib.pyplot as plt
import statistics
import argparse
import astropy
from astropy.table import Table,join
from astropy.time import Time
from astropy.io import fits

import LSS.common_tools as common


parser = argparse.ArgumentParser()
#parser.add_argument("--type", help="tracer type to be selected")
basedir='/global/cfs/cdirs/desi/survey/catalogs'
parser.add_argument("--basedir", help="base directory for input/output",default=basedir)
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='DA02')
parser.add_argument("--verspec",help="version for redshifts",default='guadalupe')
parser.add_argument("--tracer",help="tracer type (e.g., LRG)",default='LRG')

args = parser.parse_args()
basedir = args.basedir
version = args.version
survey  = args.survey
specver = args.verspec
tp = args.tracer



#ff = fitsio.read(filepathLF)
#hdul = fits.open(filepathLF)
#ff2 = fitsio.read(filepathBGS)
#hdul = fits.open(filepathBGS)

if survey != 'SV3':
    if tp == 'LRG':
        bit = 1 #for selecting LRG
    if tp == 'ELG':
        bit = 2
    if tp == 'QSO':
        bit = 4
    if tp == 'BGS_ANY':
        zf = basedir+'/'+survey+'/LSS/'+specver+'/datcomb_bright_tarspecwdup_zdone.fits'
        dz = Table(fitsio.read(zf))

        desitarg = 'BGS_TARGET'
        wtype = dz[desitarg] > 0#((dz[desitarg] & bit) > 0)
    else:    
        zf = basedir+'/'+survey+'/LSS/'+specver+'/datcomb_dark_tarspecwdup_zdone.fits'
        dz = Table(fitsio.read(zf))
        desitarg = 'DESI_TARGET'
        wtype = ((dz[desitarg] & bit) > 0)
    print(len(dz[wtype]))
    #dz = dz[wtype&wg]
    dz = dz[wtype]

    dz = common.cut_specdat(dz)
    from LSS.globals import main
    pars = main(args.type,args.verspec)


else:
    sys.exit('not written for SV3 yet')
    zf = basedir+'/'+survey+'/LSS/'+specver+'/datcomb_dark_tarspecwdup_Alltiles.fits'
    dz = Table(fitsio.read(zf))
    desitarg = 'SV3_DESI_TARGET'
    bit = 1 #for selecting LRG
    wtype = ((dz[desitarg] & bit) > 0)
    print(len(dz[wtype]))
    #dz = dz[wtype&wg]
    dz = dz[wtype]
    wz = dz['ZWARN'] != 999999 #this is what the null column becomes
    wz &= dz['ZWARN']*0 == 0 #just in case of nans
    wz &= dz['COADD_FIBERSTATUS'] == 0
    ff = dz[wz]

    zf = basedir+'/'+survey+'/LSS/'+specver+'/datcomb_bright_tarspecwdup_Alltiles.fits'
    dz = Table(fitsio.read(zf))
    desitarg = 'SV3_BGS_TARGET'
    wtype = dz[desitarg] > 0#((dz[desitarg] & bit) > 0)
    print(len(dz[wtype]))
    #dz = dz[wtype&wg]
    dz = dz[wtype]
    wz = dz['ZWARN'] != 999999 #this is what the null column becomes
    wz &= dz['ZWARN']*0 == 0 #just in case of nans
    wz &= dz['COADD_FIBERSTATUS'] == 0

    ff2 = dz[wz]

z_tot = dz['ZWARN'] != 999999
z_tot &= dz['ZWARN']*0 == 0


if tp == 'LRG':
    z_suc= dz['ZWARN']==0
    z_suc &= dz['DELTACHI2']>15
    z_suc &= dz['Z']<1.5

if tp == 'ELG':
    o2f = fitsio.read(pars.elgzf,columns=['TARGETID','LOCATION','TILEID','OII_FLUX','OII_FLUX_IVAR'])
    dz = join(dz,o2f,keys=['TARGETID','TILEID','LOCATION'])
    o2c = np.log10(dz['OII_FLUX'] * np.sqrt(dz['OII_FLUX_IVAR']))+0.2*np.log10(dz['DELTACHI2'])
    z_suc = o2c > 0.9

if tp == 'QSO':
    arz = fitsio.read(pars.qsozf,columns=['TARGETID','LOCATION','TILEID','Z','Z_QN'])

    arz['TILEID'] = arz['TILEID'].astype(int)
    dz = join(dz,arz,keys=['TARGETID','TILEID','LOCATION'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
    #dz['Z'].name = 'Z_RR' #rename the original redrock redshifts
    #dz['Z_QF'].name = 'Z' #the redshifts from the quasar file should be used instead

    z_suc = dz['Z_QF'].masked == False


if tp == 'BGS_ANY':    
    z_suc = dz['ZWARN']==0
    z_suc &= dz['DELTACHI2']>40

#print(len(ff[z_suc]),len(ff[z_tot]))
print("zsuccess rate for "+tp,len(ff[z_suc])/len(ff[z_tot]))
sys.exit()
cat1 = Table(ff[z_tot])

full=Table()
full['FIBER'] = np.arange(5000)

fiberstats = Table()
fiberstats['FIBER'], fiberstats['n_tot'] = np.unique(ff['FIBER'][z_tot], return_counts=True)
#fiberstats.sort('n_tot')

tt = Table()
tt['FIBER'], tt['n_suc'] = np.unique(ff['FIBER'][z_suc], return_counts=True)

fiberstats1 = join(fiberstats, tt, keys='FIBER', join_type='outer').filled(0)
fiberstats1 = join(fiberstats1,full, keys='FIBER',join_type='outer').filled(0)
#fiberstats1['frac_suc'] = fiberstats1['n_suc']/fiberstats1['n_tot']


z_tot = ff2['ZWARN'] != 999999
z_tot &= ff2['ZWARN']*0 == 0
#print(len(ff2[z_suc]),len(ff2[z_tot]))
print("zsuccess rate for BGS=",len(ff2[z_suc])/len(ff2[z_tot]))
cat2 = Table(ff2[z_tot])

fiberstats2 = Table()
fiberstats2['FIBER'], fiberstats2['n_tot'] = np.unique(ff2['FIBER'][z_tot], return_counts=True)
#fiberstats.sort('n_tot')

tt2 = Table()
tt2['FIBER'], tt2['n_suc'] = np.unique(ff2['FIBER'][z_suc], return_counts=True)
fiberstats2 = join(fiberstats2, tt2, keys='FIBER', join_type='outer').filled(0)
fiberstats2 = join(fiberstats2,full, keys='FIBER',join_type='outer').filled(0)
#fiberstats2['frac_suc'] = fiberstats2['n_suc']/fiberstats2['n_tot']


fstats_comb = Table()
fstats_comb['Fiber']=np.arange(5000)
fstats_comb['n_tot']=np.arange(5000)
fstats_comb['n_suc']=np.arange(5000)
for fiber in fstats_comb['Fiber']:
    m1=fiberstats1['FIBER']==fiber
    m2=fiberstats2['FIBER']==fiber
    fstats_comb['n_tot'][fiber] = fiberstats1['n_tot'][m1]+fiberstats2['n_tot'][m2]
    fstats_comb['n_suc'][fiber] = fiberstats1['n_suc'][m1]+fiberstats2['n_suc'][m2]

mask0= fstats_comb['n_tot']>1
fstats_comb=fstats_comb[mask0]
fstats_comb['frac_suc']=fstats_comb['n_suc']/fstats_comb['n_tot']
#fstats_comb 

error_floor = True

n, p = fstats_comb['n_tot'].copy(), fstats_comb['frac_suc'].copy()
if error_floor:
    p1 = np.maximum(1-p, 1/n)  # error floor
else:
    p1 = p
fstats_comb['frac_suc_err'] = np.clip(np.sqrt(n * p * (1-p))/n, np.sqrt(n * p1 * (1-p1))/n, 1)

#print("Removed fibers for having only 1 obs:\n",fstats_comb['FIBER'][ntotmask])
mean = np.sum(fstats_comb['n_suc'])/np.sum(fstats_comb['n_tot'])
fstats_comb['check'] =(mean - fstats_comb['frac_suc'])/fstats_comb['frac_suc_err']
fstats_comb.sort('frac_suc')
#fstats_comb



#mean = np.sum(fstats_comb['n_suc'])/np.sum(fstats_comb['n_tot'])
n = 3
maskcheck = fstats_comb['check']>n
print(fstats_comb)
#np.savetxt(basedir+'/'+survey+'/LSS/'+specver+'/LSScats/'+version+"/lrg+bgs_"+str(n)+"sig_bad_fibers.txt",fstats_comb[maskcheck]['Fiber'],fmt='%i')
fn = basedir+'/'+survey+'/LSS/'+specver+"/lrg+bgs_"+str(n)+"sig_bad_fibers.txt"
np.savetxt(fn,fstats_comb[maskcheck]['Fiber'],fmt='%i')
print('saved results to '+fn)