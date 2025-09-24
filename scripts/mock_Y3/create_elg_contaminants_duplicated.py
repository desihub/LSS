import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.io import fits
import time
from astropy.coordinates import SkyCoord
from astropy import units as u
import yaml
import sys
import astropy
#from regressis import footprint
#from regressis.weight import PhotoWeight
import os
import fitsio
from astropy.table import Table,unique,join,vstack
import LSS.common_tools as common
from LSS.globals import main
import h5py
#from functions_AT_new import *

basedir = '/global/cfs/cdirs/desi/survey/catalogs/'
survey = 'DA2'
data = 'LSS'
verspec = 'loa-v1'
version = 'v1.1'
mock_path = '/global/cfs/cdirs/desi/mocks/cai/holi/v3.00/seed0201/holi_ELG_v3.00_GCcomb_clustering.dat.h5'
path_out = '/global/cfs/cdirs/desi/users/akrolew/'
name_out = 'holi_ELG_v3.00_GCcomb_clustering_append_failures.h5'
tp = 'ELGnotqso'
isotropic_number_density = 97. # Density of isotropic component, already accounted for
np.random.seed(42)

zmin = 1.1
zmax = 1.6

def cut_to_region_and_zrange(fname_NGC, fname_SGC):
	
	data_file1 = fits.open(fname_NGC)[1].data
	data_file2 = fits.open(fname_SGC)[1].data

	#if (region == 'NGC') or (region == 'S_NGC'):
	#	ra = data_file1['RA']
	#	dec = data_file1['DEC']
	#	z = data_file1['Z']
	#	weight = data_file1['WEIGHT']
	#	weight_sys = data_file1['WEIGHT_SYS']
	#	photsys = data_file1['PHOTSYS']
	#	targetid = data_file1['TARGETID']
	#elif (region == 'SGC') or (region == 'S_SGC-noDES') or (region == 'DES'):
	#	ra = data_file2['RA']
	#	dec = data_file2['DEC']
	#	z = data_file2['Z']
	#	weight = data_file2['WEIGHT']
	#	weight_sys = data_file2['WEIGHT_SYS']
	#	photsys = data_file2['PHOTSYS']
	#	targetid = data_file2['TARGETID']
	#else:
	ra = np.concatenate((data_file1['RA'],data_file2['RA']))
	dec = np.concatenate((data_file1['DEC'],data_file2['DEC']))
	z = np.concatenate((data_file1['Z'],data_file2['Z']))
	weight = np.concatenate((data_file1['WEIGHT'],data_file2['WEIGHT']))
	weight_sys = np.concatenate((data_file1['WEIGHT_SYS'],data_file2['WEIGHT_SYS']))
	photsys = np.concatenate((data_file1['PHOTSYS'],data_file2['PHOTSYS']))
	targetid = np.concatenate((data_file1['TARGETID'],data_file2['TARGETID']))

	# Redshift cut
	len_before_cut = np.shape(ra)[-1]
	ra = ra[(z >= zmin) & (z <= zmax)]
	dec = dec[(z >= zmin) & (z <= zmax)]
	weight = weight[(z >= zmin) & (z <= zmax)]
	weight_sys = weight_sys[(z >= zmin) & (z <= zmax)]
	photsys = photsys[(z >= zmin) & (z <= zmax)]
	targetid = targetid[(z >= zmin) & (z <= zmax)]
	z = z[(z >= zmin) & (z <= zmax)]
	len_after_cut = np.shape(ra)[-1]
	norm_fac = len_after_cut / len_before_cut



	return ra, dec, weight, weight_sys, z, photsys, targetid, norm_fac

indir = basedir+'/'+survey+'/'+data+'/'+verspec+'/LSScats/'+version+'/'
lss_dir = indir + 'nonKP/'
tracer = tp
for i in range(18):
	ra, dec, weight, weight_sys, z, photsys, targetid, _ = cut_to_region_and_zrange(lss_dir + tracer +'_NGC_%i_clustering.ran.fits' % i, lss_dir + tracer +'_SGC_%i_clustering.ran.fits' % i)

	#if ran_weight_opt == 'default':
	#    weight = weight / weight_sys
	#elif ran_weight_opt == 'default_addLIN':
	#    print("Warining! Randoms right now don't have linear weights!!")
	#    pass
	#elif ran_weight_opt == 'default_addRF':
	#    pass
	#elif ran_weight_opt == 'ones':
	#    weight = np.ones_like(weight)
	

	if i == 0:
		rand_ra = ra
		rand_dec = dec
		rand_w = weight
	else:
		rand_ra = np.concatenate( (rand_ra, ra) )
		rand_dec = np.concatenate( (rand_dec, dec) )
		rand_w = np.concatenate( (rand_w, weight) )
	print(i)



full = Table(fitsio.read(indir+tp+'_full_HPmapcut.dat.fits'))
## the selection of valid samples
sel_obs = full['ZWARN'] != 999999
sel_obs&= full['ZWARN']*0 == 0

selection = sel_obs

gz = common.goodz_infull(tp[:3],full)
selection_gz = selection&gz

data = full[selection&~gz]

ra = data['RA']
dec = data['DEC']



coords_data = SkyCoord(ra=ra*u.deg,dec=dec*u.deg)
coords_data = coords_data.transform_to(astropy.coordinates.Galactic())
coords_ran = SkyCoord(ra=rand_ra*u.deg,dec=rand_dec*u.deg)
coords_ran = coords_ran.transform_to(astropy.coordinates.Galactic())

weight = 1/data['FRACZ_TILELOCID'] * 1/data['FRAC_TLOBS_TILES']
weight[np.isinf(weight)] = 0

data_pix = hp.ang2pix(nside, coords_data.l.deg, coords_data.b.deg, lonlat=True)
data_map = np.bincount(data_pix, minlength=12*nside**2, weights=weight)


ran_pix = hp.ang2pix(nside, coords_ran.l.deg,coords_ran.b.deg, lonlat=True)
ran_map = np.bincount(ran_pix,weights=np.ones_like(rand_w),minlength=12*nside**2)
ran_comp = ran_map / (2500 * 18 * 41253 / (12 * nside**2))

number_density = data_map / (ran_comp * 41253./(12*nside**2.))


sel_fraction = isotropic_number_density / number_density

random_number = np.random.uniform(len(data))
non_isotropic_sel = np.where(random_number > sel_fraction[data_pix])

unclassified = data[non_isotropic_sel]
weight = weight[non_isotropic_sel]

sim_data = h5py.File(mock_path,'r')
sim_ra = sim_data['RA'][:]
sim_dec = sim_data['DEC'][:]
sim_nx = sim_data['NX'][:]
sim_z = sim_data['Z'][:]



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
sim_nx = np.concatenate((sim_nx, np.zeros_like(unclassified_ra)))
sim_z = np.concatenate((sim_z, np.zeros_like(unclassified_ra)))
	
with h5py.File(path_out + name_out,'w') as f:
	# Create a dataset named 'my_dataset'
	# The data argument directly writes the NumPy array to the dataset
	f.create_dataset('RA', data=sim_ra)
	f.create_dataset('DEC',data=sim_dec)
	f.create_dataset('Z',data=sim_z)
	f.create_dataset('NX',data=sim_nx)