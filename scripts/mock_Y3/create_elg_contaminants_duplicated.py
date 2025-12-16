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
import errno

def test_dir(value):
    if not os.path.exists(value):
        try:
            os.makedirs(value, 0o755)
            print('made ' + value)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise


basedir = '/global/cfs/cdirs/desi/survey/catalogs/'
survey = 'DA2'
data = 'LSS'
verspec = 'loa-v1'
version = 'v2'
tp = 'ELGnotqso'
#mock_path = '/global/cfs/cdirs/desi/mocks/cai/holi/v3.00/seed0201/holi_ELG_v3.00_GCcomb_clustering.dat.h5'
path_out = f'/global/cfs/projectdirs/desi/mocks/cai/contaminants/{survey}/{verspec}/{version}/{tp}/noveto/'
name_out = 'contaminants_rea{NUMREA}.fits'
test_dir(path_out)

#path_out = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/ELG_v5/z0.950/'
#name_out = 'contaminants_forFA0_Y3_noimagingmask_applied.fits'
isotropic_number_density = 97. # Density of isotropic component, already accounted for
np.random.seed(42)

zmin = 0.0
zmax = 5.0

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


print(tracer)
full = Table(fitsio.read(indir+tp+'_full_noveto.dat.fits'))

##full = Table(fitsio.read(indir+tp+'_full_HPmapcut.dat.fits'))

if 'FRAC_TLOBS_TILES' not in full.columns:
    print('calculating FRAC_TLOBS_TILES')

    compa = []
    fractl = []
    tll = []
    ti = 0
    full.sort('TILES')
    nts = len(np.unique(full['TILES']))
    tlsl = full['TILES']
    tlslu = np.unique(tlsl)
    laa = full['LOCATION_ASSIGNED']
    lta = full['TILELOCID_ASSIGNED']
        #print('TILELOCID_ASSIGNED',np.unique(ff['TILELOCID_ASSIGNED'],return_counts=True),len(ff))

        # for tls in np.unique(dz['TILES']): #this is really slow now, need to figure out a better way
    i = 0
    tot = 0
    atot = 0
    tltot = 0
    while i < len(full):
        tls = []
        tlis = []
        nli = 0 #initialize total available per tile group
        nai = 0 #initialize total assigned
        nti = 0 #initialize total at location where something of the same type was assigned

        while tlsl[i] == tlslu[ti]:
            nli += 1
            nai += laa[i] #laa is true/false assigned
            nti += lta[i] #lta is true/false something of the same type was assigned
            i += 1
            if i == len(full):
                break

        if ti % 100000 == 0:
            print('at tiles ' + str(ti) + ' of ' + str(nts))

        tot += nli
        atot += nai
        tltot += nti
        cp = nai / nli #
        fract = nti/nli
            # print(tls,cp,no,nt)
        compa.append(cp)
        fractl.append(fract)
        tll.append(tlslu[ti])
        ti += 1
        #print(tot,atot,tltot)

    comp_dicta = dict(zip(tll, compa))
    fract_dicta = dict(zip(tll, fractl))

    fcompa = []
    fracta = []
    for tl in full['TILES']:
        fcompa.append(comp_dicta[tl])
        fracta.append(fract_dicta[tl])
    full['COMP_TILE'] = np.array(fcompa)
    full['FRAC_TLOBS_TILES'] = np.array(fracta)

## the selection of valid samples
sel_obs = full['ZWARN'] != 999999
sel_obs&= full['ZWARN']*0 == 0

selection = sel_obs

gz = common.goodz_infull(tp[:3],full, zcol='Z')
selection_gz = selection&gz

data = full[selection&~gz]
ra = data['RA']
dec = data['DEC']



coords_data = SkyCoord(ra=ra*u.deg,dec=dec*u.deg)
coords_data = coords_data.transform_to(astropy.coordinates.Galactic())
coords_ran = SkyCoord(ra=rand_ra*u.deg,dec=rand_dec*u.deg)
coords_ran = coords_ran.transform_to(astropy.coordinates.Galactic())

weighto = 1/data['FRACZ_TILELOCID'] * 1/data['FRAC_TLOBS_TILES']
weighto[np.isinf(weighto)] = 0

nside = 32

data_pix = hp.ang2pix(nside, coords_data.l.deg, coords_data.b.deg, lonlat=True)
data_map = np.bincount(data_pix, minlength=12*nside**2, weights=weighto)


ran_pix = hp.ang2pix(nside, coords_ran.l.deg,coords_ran.b.deg, lonlat=True)
ran_map = np.bincount(ran_pix,weights=np.ones_like(rand_w),minlength=12*nside**2)
ran_comp = ran_map / (2500 * 18 * 41253 / (12 * nside**2))

number_density = data_map / (ran_comp * 41253./(12*nside**2.))


sel_fraction = isotropic_number_density / number_density



for rea in range(8, 100):
    print('NUMREA', rea)
    random_number = np.random.uniform(size=len(data))
    non_isotropic_sel = np.where(random_number > sel_fraction[data_pix])

    unclassified = data[non_isotropic_sel]

#vals,bins =np.histogram(unclassified['Z_not4clus'],bins=100)
#bin_centers = bins[:-1] + np.diff(bins)/2
#np.savetxt('zhist_unclassified.txt',np.array([bin_centers,vals]).T)
#vals,bins =np.histogram(data['Z_not4clus'],bins=100)
#bin_centers = bins[:-1] + np.diff(bins)/2
#np.savetxt('zhist_elg.txt',np.array([bin_centers,vals]).T)

#exit()

    weight = weighto[non_isotropic_sel]

    '''sim_data = h5py.File(mock_path,'r')
    sim_ra = sim_data['RA'][:]
    sim_dec = sim_data['DEC'][:]
    #sim_nx = sim_data['NX'][:]
    sim_z = sim_data['Z'][:]
    '''


    int_weight = np.round(weight)

    argsort_weight_diff = np.argsort(weight - int_weight)

    for i in range(len(argsort_weight_diff)):
        if np.sum(int_weight) < np.sum(weight):
            #print(np.sum(int_weight))
            int_weight[argsort_weight_diff[i]] += 1
        else:
            break


    unclassified_ra = np.array([])
    unclassified_dec = np.array([])
    unclassified_z = np.array([])
    unclassified_desitarget = np.array([])

    for i in range(len(weight)):
        if int_weight[i] > 1:
            offset = np.random.uniform(0, 0.025, int(int_weight[i]-1))
            offset = np.concatenate((np.array([0]), offset))
            phase = np.random.uniform(0, 2*np.pi, int(int_weight[i]-1))
            phase = np.concatenate((np.array([0]), phase))
            delta_ra = offset * np.cos(phase) / np.cos(unclassified['DEC'][i])
            delta_dec = offset * np.sin(phase)
            #print(unclassified['RA'][i], unclassified['DEC'][i], unclassified['Z_not4clus'][i], unclassified['DESI_TARGET'][i])
            unclassified_ra = np.concatenate((unclassified_ra, unclassified['RA'][i] + delta_ra))
            unclassified_dec = np.concatenate((unclassified_dec, unclassified['DEC'][i] + delta_dec))
            unclassified_z = np.concatenate((unclassified_z, np.atleast_1d(unclassified['Z'][i])))
            unclassified_desitarget = np.concatenate((unclassified_desitarget, np.atleast_1d(unclassified['DESI_TARGET'][i])))

        else:
            unclassified_ra = np.concatenate((unclassified_ra, [unclassified['RA'][i]]))
            unclassified_dec = np.concatenate((unclassified_dec, [unclassified['DEC'][i]]))
            unclassified_z = np.concatenate((unclassified_z, [unclassified['Z'][i]]))
            unclassified_desitarget = np.concatenate((unclassified_desitarget, [unclassified['DESI_TARGET'][i]]))

    sim_ra = unclassified_ra #np.concatenate((sim_ra, unclassified_ra))
    sim_dec = unclassified_dec # np.concatenate((sim_dec, unclassified_dec))
    #sim_nx = np.concatenate((sim_nx, np.zeros_like(unclassified_ra)))

    sim_z = unclassified_z # np.zeros_like(unclassified_ra) # np.concatenate((sim_z, np.zeros_like(unclassified_ra)))

    savefits=True
    if savefits:
        col1 = fits.Column(name='RA', array=sim_ra, format='D')
        col2 = fits.Column(name='DEC', array=sim_dec, format='D')
        col3 = fits.Column(name='Z', array=sim_z, format='D')
        col4 = fits.Column(name='DESI_TARGET', array=unclassified_desitarget, format='K')
        cols = fits.ColDefs([col1, col2, col3, col4])
        tbhdu = fits.BinTableHDU.from_columns(cols)
        tbhdu.writeto(path_out + name_out.format(NUMREA=rea), overwrite=True)
    else:

        with h5py.File(path_out + name_out,'w') as f:
            # Create a dataset named 'my_dataset'
        # The data argument directly writes the NumPy array to the dataset
            f.create_dataset('RA', data=sim_ra)
            f.create_dataset('DEC',data=sim_dec)
            f.create_dataset('Z',data=sim_z)
    #        f.create_dataset('NX',data=sim_nx)
