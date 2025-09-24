import numpy as np
import os
import sys
from astropy.io import fits
import h5py

import fitsio
from astropy.table import Table,unique,join,vstack
import LSS.common_tools as common
from LSS.globals import main

# Script to append quasar contaminants to the high-fidelity mocks
# Starts by identifying contaminants in the full file of the data
# and then adds them to the mock

append_stars = True
append_unclassified = True
append_galaxies = True


basedir = '/global/cfs/cdirs/desi/survey/catalogs/'
survey = 'DA2'
data = 'LSS'
verspec = 'loa-v1'
version = 'v1.1'
mock_path = '/global/cfs/cdirs/desi/mocks/cai/holi/v3.00/seed0201/holi_QSO_v3.00_GCcomb_clustering.dat.h5'
path_out = '/global/cfs/cdirs/desi/users/akrolew/'
name_out = 'holi_QSO_v3.00_GCcomb_clustering_append_stars_unclassified_galaxies.h5'
tp = 'QSO'

#path_star = '/global/cfs/cdirs/desi/users/akrolew/%s_%s_%s_%s_stars.fits' % (tp, survey, verspec, version)
#path_unclassified = '/global/cfs/cdirs/desi/users/akrolew/%s_%s_%s_%s_unclassified.fits' % (tp, survey, verspec, version)
#path_galaxies = '/global/cfs/cdirs/desi/users/akrolew/%s_%s_%s_%s_galaxies.fits' % (tp, survey, verspec, version)


indir = basedir+'/'+survey+'/'+data+'/'+verspec+'/LSScats/'+version+'/'


full = Table(fitsio.read(indir+tp+'_full_HPmapcut.dat.fits'))
## the selection of valid samples
sel_obs = full['ZWARN'] != 999999
sel_obs&= full['ZWARN']*0 == 0

selection = sel_obs

gz = common.goodz_infull(tp[:3],full)
selection_gz = selection&gz

emline = fits.open(basedir+'/'+survey+'/'+data+'/'+verspec+'/emlin_catalog.fits')[1].data

ss = np.searchsorted(sorted(emline['TARGETID']), full['TARGETID'])

argss = np.argsort(emline['TARGETID'])

oii_flux = emline['OII_FLUX'][argss][ss]
oii_flux_ivar = emline['OII_FLUX_IVAR'][argss][ss]

oiii_flux = emline['OIII_FLUX'][argss][ss]
oiii_flux_ivar = emline['OIII_FLUX_IVAR'][argss][ss]

dchi_cut = 30
o2c_cut = 0.9
oiii_cut = 5
selgal = ((selection  & ~gz & (full['Z_RR'] > 0.01)  & (full['Z_RR'] < 1.625)) & 
	((full['DELTACHI2'] > dchi_cut) | (np.log10(oii_flux * oii_flux_ivar**0.5) > o2c_cut - 0.2 * np.log10(full['DELTACHI2']))
	| (oiii_flux * oiii_flux_ivar**0.5 > oiii_cut))
	)  #...now 40,1.2,5...earlier 30,0.9,5
	
selstar = (selection & ~gz & (full['Z_RR'] < 0.01))
                
                
print('Fraction of stars',np.sum(selstar) / np.sum(selection))
print('Fraction of galaxies',np.sum(selgal) / np.sum(selection))
print('Fraction of qsos',np.sum(selection&gz)/np.sum(selection))
print('Fraction of junk',1-( np.sum(selstar) + np.sum(selgal) + np.sum(selection&gz) )/np.sum(selection))

#sim_data = fits.open(mock_path)[1].data
#with h5py.File(mock_path, 'r') as f:
#	data = f['my_dataset'][:]
#	print("Data from 'my_dataset':", data)
sim_data = h5py.File(mock_path,'r')
sim_ra = sim_data['RA'][:]
sim_dec = sim_data['DEC'][:]
sim_nx = sim_data['NX'][:]
sim_z = sim_data['Z'][:]
np.random.seed(123)
if append_stars:
	stars = full[selstar]
	weight = 1/full[selstar]['FRACZ_TILELOCID'] * 1/full[selstar]['FRAC_TLOBS_TILES']
	weight[np.isinf(weight)] = 0
	
	int_weight = np.round(weight)
	
	argsort_weight_diff = np.argsort(weight - int_weight)
	
	for i in range(len(argsort_weight_diff)):
		if np.sum(int_weight) < np.sum(weight):
			print(np.sum(int_weight))
			int_weight[argsort_weight_diff[i]] += 1
		else:
			break
	
	
	stars_ra = np.array([])
	stars_dec = np.array([])
	for i in range(len(weight)):
		if int_weight[i] > 1:
			offset = np.random.uniform(0, 0.025, int(int_weight[i]-1))
			offset = np.concatenate((np.array([0]), offset))
			phase = np.random.uniform(0, 2*np.pi, int(int_weight[i]-1))
			phase = np.concatenate((np.array([0]), phase))
			delta_ra = offset * np.cos(phase) / np.cos(stars['DEC'][i])
			delta_dec = offset * np.sin(phase)
			stars_ra = np.concatenate((stars_ra, stars['RA'][i] + delta_ra))
			stars_dec = np.concatenate((stars_dec, stars['DEC'][i] + delta_dec))
		else:
			stars_ra = np.concatenate((stars_ra, [stars['RA'][i]]))
			stars_dec = np.concatenate((stars_dec, [stars['DEC'][i]]))

	sim_ra = np.concatenate((sim_ra, stars_ra))
	sim_dec = np.concatenate((sim_dec, stars_dec))
	sim_nx = np.concatenate((sim_nx, np.zeros_like(stars_ra))
	sim_z = np.concatenate((sim_z, np.zeros_like(stars_ra))
	
	# star_out = Table(np.array([-9 * np.ones(len(stars_ra)).astype('int'),
	# 	-9 * np.ones(len(stars_ra)).astype('int'),
	# 	stars_ra,
	# 	stars_dec,
	# 	np.zeros_like(stars_ra),
	# 	np.zeros_like(stars_dec),
	# 	np.array(['STAR'] * len(stars_ra)),
	# 	np.zeros_like(stars_ra),
	# 	np.zeros_like(stars_ra)]).T) #names=('GALAXYID',
	# 	#'PID','RA','DEC','Z_NORSD','Z','TRACER_TYPE',
	# 	#'NX','WEIGHT_FKP'))
	# sim_data = vstack((Table(sim_data), star_out))
	#star_out.write(path_star)
if append_unclassified:
	unclassified = full[selection & ~selstar & ~selgal & ~(selection &gz)]
	weight = 1/full[selection & ~selstar & ~selgal & ~(selection &gz)]['FRACZ_TILELOCID'] * 1/full[selection & ~selstar & ~selgal & ~(selection &gz)]['FRAC_TLOBS_TILES']
	weight[np.isinf(weight)] = 0
	
	int_weight = np.round(weight)
	
	argsort_weight_diff = np.argsort(weight - int_weight)
	
	for i in range(len(argsort_weight_diff)):
		if np.sum(int_weight) < np.sum(weight):
			print(np.sum(int_weight))
			int_weight[argsort_weight_diff[i]] += 1
		else:
			break
	
	
	unclassified_ra = np.array([])
	unclassified_dec = np.array([])
	for i in range(len(weight)):
		if int_weight[i] > 1:
			offset = np.random.uniform(0, 0.025, int(int_weight[i]-1))
			offset = np.concatenate((np.array([0]), offset))
			phase = np.random.uniform(0, 2*np.pi, int(int_weight[i]-1))
			phase = np.concatenate((np.array([0]), phase))
			delta_ra = offset * np.cos(phase) / np.cos(unclassified['DEC'][i])
			delta_dec = offset * np.sin(phase)
			unclassified_ra = np.concatenate((unclassified_ra, unclassified['RA'][i] + delta_ra))
			unclassified_dec = np.concatenate((unclassified_dec, unclassified['DEC'][i] + delta_dec))
		else:
			unclassified_ra = np.concatenate((unclassified_ra, [unclassified['RA'][i]]))
			unclassified_dec = np.concatenate((unclassified_dec, [unclassified['DEC'][i]]))

	sim_ra = np.concatenate((sim_ra, unclassified_ra))
	sim_dec = np.concatenate((sim_dec, unclassified_dec))
	sim_nx = np.concatenate((sim_nx, np.zeros_like(unclassified_ra))
	sim_z = np.concatenate((sim_z, np.zeros_like(unclassified_ra))


	# unclassified_out = Table(np.array([-9 * np.ones(len(unclassified_ra)).astype('int'),
	# 	-9 * np.ones(len(unclassified_ra)).astype('int'),
	# 	unclassified_ra,
	# 	unclassified_dec,
	# 	np.zeros_like(unclassified_ra),
	# 	np.zeros_like(unclassified_dec),
	# 	np.array(['UNKNOWN'] * len(unclassified_ra)),
	# 	np.zeros_like(unclassified_ra),
	# 	np.zeros_like(unclassified_ra)]).T) #names=('GALAXYID',
	# 	#'PID','RA','DEC','Z_NORSD','Z','TRACER_TYPE',
	# 	#'NX','WEIGHT_FKP'))
	# #sim_data = vstack((Table(sim_data), unclassified_out))
	# unclassified_out.write(path_unclassified)

if append_galaxies:
	galaxies = full[selgal]
	weight = 1/full[selgal]['FRACZ_TILELOCID'] * 1/full[selgal]['FRAC_TLOBS_TILES']
	weight[np.isinf(weight)] = 0
	
	int_weight = np.round(weight)
	
	argsort_weight_diff = np.argsort(weight - int_weight)
	
	for i in range(len(argsort_weight_diff)):
		if np.sum(int_weight) < np.sum(weight):
			print(np.sum(int_weight))
			int_weight[argsort_weight_diff[i]] += 1
		else:
			break
	
	
	galaxies_ra = np.array([])
	galaxies_dec = np.array([])
	for i in range(len(weight)):
		if int_weight[i] > 1:
			offset = np.random.uniform(0, 0.025, int(int_weight[i]-1))
			offset = np.concatenate((np.array([0]), offset))
			phase = np.random.uniform(0, 2*np.pi, int(int_weight[i]-1))
			phase = np.concatenate((np.array([0]), phase))
			delta_ra = offset * np.cos(phase) / np.cos(galaxies['DEC'][i])
			delta_dec = offset * np.sin(phase)
			galaxies_ra = np.concatenate((galaxies_ra, galaxies['RA'][i] + delta_ra))
			galaxies_dec = np.concatenate((galaxies_dec, galaxies['DEC'][i] + delta_dec))
		else:
			galaxies_ra = np.concatenate((galaxies_ra, [galaxies['RA'][i]]))
			galaxies_dec = np.concatenate((galaxies_dec, [galaxies['DEC'][i]]))

	sim_ra = np.concatenate((sim_ra, galaxies_ra))
	sim_dec = np.concatenate((sim_dec, galaxies_dec))
	sim_nx = np.concatenate((sim_nx, np.zeros_like(galaxies_ra))
	sim_z = np.concatenate((sim_z, np.zeros_like(galaxies_ra))


	# galaxies_out = Table(np.array([-9 * np.ones(len(galaxies_ra)).astype('int'),
	# 	-9 * np.ones(len(galaxies_ra)).astype('int'),
	# 	galaxies_ra,
	# 	galaxies_dec,
	# 	np.zeros_like(galaxies_ra),
	# 	np.zeros_like(galaxies_dec),
	# 	np.array(['GALAXY'] * len(galaxies_ra)),
	# 	np.zeros_like(galaxies_ra),
	# 	np.zeros_like(galaxies_ra)]).T) #names=('GALAXYID',
	# 	#'PID','RA','DEC','Z_NORSD','Z','TRACER_TYPE',
	# 	#'NX','WEIGHT_FKP'))
	# #sim_data = vstack((Table(sim_data), galaxies_out))
	# galaxies_out.write(path_galaxies)

#sim_data.write(path_out + name_out)

with h5py.File(path_out + name_out 'w') as f:
	# Create a dataset named 'my_dataset'
	# The data argument directly writes the NumPy array to the dataset
	f.create_dataset('RA', data=sim_ra)
	f.create_dataset('DEC',data=sim_dec)
	f.create_dataset('Z',data=sim_z)
	f.create_dataset('NX',data=sim_nx)