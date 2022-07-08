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
parser.add_argument("--verspec_new",help="version for redshifts",default='newQSOtemp')
parser.add_argument("--tracer",help="tracer type(s) (e.g., LRG)",default='all')

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



for tp in tracers:
	notqso = ''
	if survey == 'DA02':
		if tp == 'LRG':
			bit = 1 #for selecting LRG
		if tp == 'ELG':
			bit = 2
			notqso = 'notqso'
		if tp == 'QSO':
			bit = 4
		if tp == 'BGS_ANY':
			zf = basedir+'/'+survey+'/LSS/'+specver+'/datcomb_bright_tarspecwdup_zdone.fits'
			zf_new = basedir+'/'+survey+'/LSS/'+args.verspec_new+'/datcomb_bright_spec_zdone.fits'
			dz = Table(fitsio.read(zf))
			

			desitarg = 'BGS_TARGET'
			wtype = dz[desitarg] > 0#((dz[desitarg] & bit) > 0)
		else:    
			zf = basedir+'/'+survey+'/LSS/'+specver+'/datcomb_dark_tarspecwdup_zdone.fits'
			zf_new = basedir+'/'+survey+'/LSS/'+args.verspec_new+'/datcomb_dark_spec_zdone.fits'
			dz = Table(fitsio.read(zf))
			desitarg = 'DESI_TARGET'
			wtype = ((dz[desitarg] & bit) > 0)
			if tp == 'ELG':
				wtype &= ((dz[desitarg] & 4) == 0) #remove QSO
		print(len(dz[wtype]))
		#dz = dz[wtype&wg]
		dz = dz[wtype]

		dz = common.cut_specdat(dz)
		dz_new = Table(fitsio.read(zf_new))
		dz_new.keep_columns(['Z','ZWARN','DELTACHI2','TARGETID','TILEID','LOCATION'])
		print(len(dz))
		dz = join(dz,dz_new,keys=['TARGETID','TILEID','LOCATION'],table_names=['fid','new'])
		print(str(len(dz))+' should agree with above')
		
		 
		from LSS.globals import main
		pars = main(tp,args.verspec)
	
	elif survey == 'main':
		sys.exit(survey+' not supported yet')
		zf = basedir+'/'+survey+'/LSS/'+specver+'/datcomb_'+tp+'_tarspecwdup_zdone.fits'
		dz = Table(fitsio.read(zf))
		if tp == 'ELG':
			wtype = ((dz['DESI_TARGET'] & 4) == 0) #remove QSO
			dz = dz[wtype]
		dz = common.cut_specdat(dz)
		from LSS.globals import main
		pars = main(tp,args.verspec)
		

	elif survey == 'SV3':
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

	z_tot = dz['ZWARN_fid'] != 999999
	z_tot &= dz['ZWARN_fid']*0 == 0
	z_new = dz['ZWARN_new'] != 999999
	z_new &= dz['ZWARN_new']*0 == 0
	print('number with z to consider fid,new')
	print(len(dz[z_tot]),len(dz[z_new]))


	if tp == 'LRG':
		z_suc= dz['ZWARN_fid']==0
		z_suc &= dz['DELTACHI2_fid']>15
		z_suc &= dz['Z_fid']<1.5
		z_sucnew= dz['ZWARN_new']==0
		z_sucnew &= dz['DELTACHI2_new']>15
		z_sucnew &= dz['Z_new']<1.5
		zmin = 0.4
		zmax = 1.1

	if tp == 'ELG':
		o2f = fitsio.read(pars.elgzf,columns=['TARGETID','LOCATION','TILEID','OII_FLUX','OII_FLUX_IVAR'])
		dz = join(dz,o2f,keys=['TARGETID','TILEID','LOCATION'])
		o2c = np.log10(dz['OII_FLUX'] * np.sqrt(dz['OII_FLUX_IVAR']))+0.2*np.log10(dz['DELTACHI2_fid'])
		z_suc = o2c > 0.9
		o2f_new = fitsio.read(basedir+'/'+survey+'/LSS/'+args.verspec_new+'/emlin_catalog.fits' ,columns=['TARGETID','LOCATION','TILEID','OII_FLUX','OII_FLUX_IVAR'])
		dz = join(dz,o2f_new,keys=['TARGETID','TILEID','LOCATION'],table_names=['fid','new'])
		o2c_new = np.log10(dz['OII_FLUX_new'] * np.sqrt(dz['OII_FLUX_IVAR_new']))+0.2*np.log10(dz['DELTACHI2_new'])
		z_sucnew = o2c_new > 0.9
		zmin = 0.6
		zmax = 1.6

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

		z_suc = dz['Z'].mask == False #previous Z column should have become Z_fid
		qsozf_new = basedir+'/'+survey+'/LSS/'+args.verspec_new+'/QSO_catalog.fits'
		arz = Table(fitsio.read(qsozf_new))
		arz.keep_columns(['TARGETID','LOCATION','TILEID','Z','Z_QN'])
		arz['TILEID'] = arz['TILEID'].astype(int)
		dz = join(dz,arz,keys=['TARGETID','TILEID','LOCATION'],join_type='left',uniq_col_name='{col_name}{table_name}',table_names=['','_QF_new'])
		#print(dz.dtype.names)
		z_sucnew = dz['Z_QF_new'].mask == False
		zmin = 0.8
		zmax = 3.5


	if tp == 'BGS_ANY':    
		z_suc = dz['ZWARN']==0
		z_suc &= dz['DELTACHI2']>40
		z_sucnew = dz['ZWARN_new']==0
		z_sucnew &= dz['DELTACHI2_new']>40
		zmin = 0.01
		zmax = 0.6

	#print(len(ff[z_suc]),len(ff[z_tot]))
	print("fiducial zsuccess rate for "+tp,len(dz[z_suc&z_tot])/len(dz[z_tot]))
	print("new zsuccess rate for "+tp,len(dz[z_sucnew&z_new])/len(dz[z_new]))
	print("fraction with zsuccess in both "+tp,len(dz[z_sucnew&z_new&z_suc])/len(dz[z_new]))
	
	if tp != 'QSO':
	    plt.hist(dz['Z_fid'],histtype='step',label='fiducial',range=(zmin,zmax),bins=50)
	    plt.hist(dz['Z_new'],histtype='step',label='new',range=(zmin,zmax),bins=50)
	    plt.legend()
	    plt.xlabel('redshift')
	    plt.ylabel('# of good z in bin')
	    plt.savefig(basedir+'/'+survey+'/LSS/'+args.verspec_new+'/'+tp+notqso+'_zcompGuad.png')
	    plt.show()
	else:
	    plt.hist(dz['Z'],histtype='step',label='fiducial',range=(zmin,zmax),bins=50)
	    plt.hist(dz['Z_QF_new'],histtype='step',label='new',range=(zmin,zmax),bins=50)
	    plt.legend()
	    plt.xlabel('redshift')
	    plt.ylabel('# of good z in bin')
	    plt.savefig(basedir+'/'+survey+'/LSS/'+args.verspec_new+'/'+tp+notqso+'_zcompGuad.png')
	    plt.show()
	
	

