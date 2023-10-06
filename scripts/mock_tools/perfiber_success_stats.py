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

basedir='/global/cfs/cdirs/desi/survey/catalogs'
parser.add_argument("--tracer", help="tracer type to be selected",default='all')
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

if args.tracer == 'all':
    tracers = ['QSO','LRG','ELG','BGS_ANY']
else:
    tracers = [args.tracer]

if args.mkfiles == 'y':
    for tp in tracers:
        if survey == 'DA02':
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
        
            

        elif survey == 'SV3':
            #ys.exit('not written for SV3 yet')
            if tp != 'BGS_ANY':
                zf = basedir+'/'+survey+'/LSS/'+specver+'/datcomb_dark_tarspecwdup_Alltiles.fits'
                dz = Table(fitsio.read(zf))
                desitarg = 'SV3_DESI_TARGET'
                if tp == 'LRG':
                    bit = 1 #for selecting LRG
                if tp == 'ELG':
                    bit = 2
                if tp == 'QSO':
                    bit = 4
                wtype = ((dz[desitarg] & bit) > 0)
                if tp == 'ELG':
                    wtype &= ((dz[desitarg] & 4) == 0) #remove QSO
            else:
                zf = basedir+'/'+survey+'/LSS/'+specver+'/datcomb_bright_tarspecwdup_Alltiles.fits'
                dz = Table(fitsio.read(zf))
                desitarg = 'SV3_BGS_TARGET'
                wtype = dz[desitarg] > 0#((dz[desitarg] & bit) > 0)
            
            print(len(dz[wtype]))
            #dz = dz[wtype&wg]
            dz = dz[wtype]
            wz = dz['COADD_FIBERSTATUS'] == 0
            dz = dz[wz]

        else: 
            zf = basedir+'/'+survey+'/LSS/'+specver+'/datcomb_'+tp+'_tarspecwdup_zdone.fits'
            dz = Table(fitsio.read(zf))
            if tp == 'ELG':
                wtype = ((dz['DESI_TARGET'] & 4) == 0) #remove QSO
                dz = dz[wtype]
            dz = common.cut_specdat(dz)
            from LSS.globals import main
            pars = main(tp,args.verspec)


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
            qsozf = pars.qsozf
            if specver == 'guadalupe':
                qsozf = '/global/cfs/cdirs/desi/users/edmondc/QSO_catalog/guadalupe/QSO_cat_guadalupe_cumulative.fits'
            arz = Table(fitsio.read(qsozf))
            arz.keep_columns(['TARGETID','LOCATION','TILEID','Z','Z_QN'])
            arz['TILEID'] = arz['TILEID'].astype(int)

            #arz = fitsio.read(qsozf,columns=['TARGETID','LOCATION','TILEID','Z','Z_QN'])

            #arz['TILEID'] = arz['TILEID'].astype(int)
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
 

