import numpy as np
#!pip install astropy
#!pip install fitsio
from scipy import stats
from scipy.stats import norm
import fitsio
import glob
import os
import sys
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
parser.add_argument("--mkfiles",help="whether or not to make the files",default='n')
#parser.add_argument("--tracer",help="tracer type (e.g., LRG)",default='LRG')

args = parser.parse_args()
basedir = args.basedir
survey  = args.survey
specver = args.verspec
#tp = args.tracer



#ff = fitsio.read(filepathLF)
#hdul = fits.open(filepathLF)
#ff2 = fitsio.read(filepathBGS)
#hdul = fits.open(filepathBGS)

tracers = ['QSO','LRG','ELG','BGS_ANY']

if args.mkfiles == 'y':
    for tp in tracers:
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
                if tp == 'ELG':
                    wtype &= ((dz[desitarg] & 4) == 0) #remove QSO
            print(len(dz[wtype]))
            #dz = dz[wtype&wg]
            dz = dz[wtype]

            dz = common.cut_specdat(dz)
            from LSS.globals import main
            pars = main(tp,args.verspec)


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
            if specver == 'guadalupe':
                qsozf = '/global/cfs/cdirs/desi/users/edmondc/QSO_catalog/guadalupe/QSO_cat_guadalupe_cumulative.fits'
            arz = fitsio.read(qsozf,columns=['TARGETID','LOCATION','TILEID','Z','Z_QN'])

            arz['TILEID'] = arz['TILEID'].astype(int)
            dz = join(dz,arz,keys=['TARGETID','TILEID','LOCATION'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF'])
            #dz['Z'].name = 'Z_RR' #rename the original redrock redshifts
            #dz['Z_QF'].name = 'Z' #the redshifts from the quasar file should be used instead

            z_suc = dz['Z_QF'].mask == False


        if tp == 'BGS_ANY':    
            z_suc = dz['ZWARN']==0
            z_suc &= dz['DELTACHI2']>40

        #print(len(ff[z_suc]),len(ff[z_tot]))
        print("zsuccess rate for "+tp,len(dz[z_suc&z_tot])/len(dz[z_tot]))
        fibl,n_tot = np.unique(dz[z_tot]['FIBER'],return_counts=True)
        fiblg,n_g = np.unique(dz[z_suc&z_tot]['FIBER'],return_counts=True)
        fib_test = np.isin(fibl,fiblg)
        z_tot &= np.isin(dz['FIBER'],fibl[fib_test])
        fibl,n_tot = np.unique(dz[z_tot]['FIBER'],return_counts=True)

        if np.array_equal(fibl,fiblg):
            gfrac = n_g/n_tot
        else:
            sys.exit('need to put something in for mismatch fiber lists')

        fn = basedir+'/'+survey+'/LSS/'+specver+"/"+tp+'_zsuccess.txt'
        fo = open(fn,'w')
        for ii in range(len(fibl)):
            fo.write(str(fibl[ii])+' '+str(n_g[ii]/n_tot[ii])+' '+str(n_g[ii])+' '+str(n_tot[ii])+'\n')
        fo.close()
 
if args.survey != 'SV3':
    from LSS.globals import main
    pars = main(tp,args.verspec)   

def plot_all_petal(petal):
    for tp in tracers:
        fn = basedir+'/'+survey+'/LSS/'+specver+"/"+tp+'_zsuccess.txt'
        d = np.loadtxt(fn).transpose()
        fibl = d[0]
        f_succ = d[1]
        n_g= d[2]
        n_tot = d[3]
        
        err = np.sqrt(n_g*(1-f_succ))/n_tot

        fmin = petal*500
        fmax = (petal+1)*500
        sel = fibl >= fmin
        sel &= fibl < fmax
        bfib = pars.badfib
        sel_bfib = np.isin(fibl[sel],bfib)

        if tp == 'LRG':
            plt.errorbar(fibl[sel],f_succ[sel],err[sel],fmt='.r',label='LRG')
            plt.plot(fibl[sel][sel_bfib],f_succ[sel][sel_bfib],'kx',label='masked')
        if tp == 'ELG':
            plt.errorbar(fibl[sel]+.25,f_succ[sel],err[sel],fmt='.b',label='ELG')
            plt.plot(fibl[sel][sel_bfib],f_succ[sel][sel_bfib],'kx')
        if tp == 'QSO':
            plt.errorbar(fibl[sel]-0.25,f_succ[sel],err[sel],fmt='.g',label='QSO')
            plt.plot(fibl[sel][sel_bfib],f_succ[sel][sel_bfib],'kx')
        if tp == 'BGS_ANY':
            plt.errorbar(fibl[sel]+0.5,f_succ[sel],err[sel],fmt='.',label='BGS',color='brown')
            plt.plot(fibl[sel][sel_bfib],f_succ[sel][sel_bfib],'kx')
    plt.grid()
    plt.legend()
    plt.xlabel('FIBER')
    plt.ylabel('redshift success rate')
    plt.title(args.verspec+' petal '+str(petal))
    plt.show()

for i in range(0,10):
    plot_all_petal(i)