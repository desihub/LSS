from astropy.table import Table, vstack
from numpy.random import Generator, PCG64
from astropy.io import fits
import LSS.common_tools as common
import numpy as np
#z_lop = Table.read('/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph003/CutSky/ELG_v5/z1.175/forclustering/cutsky_abacusHF_DR2_ELG_LOP_z1p175_zcut_0p8to1p6_clustering.dat.fits')['Z']

#z_lrg = Table.read('/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph000/CutSky/LRG/z0.500/forclustering/cutsky_abacusHF_DR2_LRG_z0p500_zcut_0p4to1p1_clustering.dat.fits')['Z']

#z_qso = Table.read('/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph007/CutSky/QSO/z1.400/forclustering/cutsky_abacusHF_DR2_QSO_z1p400_zcut_0p8to3p5_clustering.dat.fits')['Z']

def return_north(ra, dec):
    '''
    given a table that already includes RA,DEC, add PHOTSYS column denoting whether
    the data is in the DECaLS ('S') or BASS/MzLS ('N') photometric region
    '''
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(ra* u.deg, dec* u.deg,frame='icrs')
    gc = c.transform_to('galactic')
    sel_ngc = gc.b > 0

    seln = dec > 32.375

    sel = seln&sel_ngc
    return sel


lrg_snap = ['z0.500'] #, 'z0.725', 'z0.950']
elg_snap = ['z0.950'] #, 'z1.175', 'z1.475']

for i in range(18):#,18):#,18):

    ##ranf = f'/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/randoms/rands_intiles_DARK_with_imagingmask_{i}.fits'
    ranf = f'/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/randoms/rands_intiles_DARK_NO_imagingmask_{i}.fits'
    #ranf = f'/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/randoms/rands_intiles_DARK_nomask_{i}.fits'
    random_file = Table.read(ranf)

    rng = Generator(PCG64())
    lrg_s = 'z0.500' #rng.choice(lrg_snap)
    lrg_rea = str(rng.choice(range(25))).zfill(3)
    
    z_lrg = Table.read(f"/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph{lrg_rea}/CutSky/LRG/{lrg_s}/forclustering/cutsky_abacusHF_DR2_LRG_{lrg_s.replace('.','p')}_zcut_0p4to1p1_clustering.dat.fits")['Z']

    rng = Generator(PCG64())
    sizerans = len(random_file)
    inds_z = rng.choice(len(z_lrg), sizerans)
    random_file['Z_LRG'] = z_lrg[inds_z]
    random_file['WEIGHT'] = np.ones(sizerans)

    
    
    rng = Generator(PCG64())
    qso_rea = str(rng.choice(range(25))).zfill(3)
    z_qso = Table.read(f"/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph{qso_rea}/CutSky/QSO/z1.400/forclustering/cutsky_abacusHF_DR2_QSO_z1p400_zcut_0p8to3p5_clustering.dat.fits")['Z']

    rng = Generator(PCG64())
    inds_z = rng.choice(len(z_qso), sizerans)
    random_file['Z_QSO'] = z_qso[inds_z]


    sel_ran = return_north(random_file['RA'], random_file['DEC'])

    ran_north = random_file[sel_ran]
    ran_south = random_file[~sel_ran]
    sizerans_north = len(ran_north)
    sizerans_south = len(ran_south)

    rng = Generator(PCG64())
    lop_s = 'z0.950' #TEMP rng.choice(elg_snap)
    lop_rea = str(rng.choice(range(25))).zfill(3)
    #elg_rea = str(rng.choice(range(25))).zfill(3)
    
    '''
    #TEMP lop_rea = str(rng.choice(range(25))).zfill(3)
    elg = Table.read(f"/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph{elg_rea}/CutSky/ELG_v5/{lop_s}/forclustering/cutsky_abacusHF_DR2_ELG_{lop_s.replace('.','p')}_zcut_0p8to1p6_clustering.dat.fits")
    #lop = Table.read(f"/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph{lop_rea}/CutSky/ELG_v5/{lop_s}/forclustering/cutsky_allzs_abacusHF_DR2_ELG_{lop_s.replace('.','p')}_zcut_0p8to1p6_clustering.dat.fits")

    ra_elg, dec_elg = elg['RA'], elg['DEC']
    sel_north = return_north(ra_elg, dec_elg)

    elg_north = elg[sel_north]
    elg_south = elg[~sel_north]

    print('size data', len(elg_north), len(elg_south))



    rng_N = Generator(PCG64())
    rng_S = Generator(PCG64())
    inds_z_N = rng_N.choice(len(elg_north), sizerans_north)
    inds_z_S = rng_S.choice(len(elg_south), sizerans_south)

    #ran_north['Z_ELG_LOP'] = lop_north['Z'][inds_z_N]
    ran_north['Z_ELG'] = elg_north['Z'][inds_z_N]
    ran_south['Z_ELG'] = elg_south['Z'][inds_z_S]
    #ran_south['Z_ELG_LOP'] = lop_south['Z'][inds_z_S]

    print('size randoms', len(ran_north), len(ran_south))
    print('randoms to data', float(len(ran_north))/len(elg_north), float(len(ran_south))/len(elg_south))


    newgen = Generator(PCG64())
    print(ran_north.dtype)
    N = len(ran_north)

    inds = newgen.choice(np.arange(N), int(len(elg_north)*float(len(ran_south))/len(elg_south)), replace=False)
    mask = np.zeros(N, dtype=bool)
    mask[inds] = True
    
    ran_north["ELG_MASK"] = mask

    N = len(ran_south)
    mask = np.ones(N, dtype=bool)
    ran_south["ELG_MASK"] = mask

    #ran_north = Table(inds)

    print('size randoms after normalization', len(ran_north), len(ran_south))
    print('randoms to data after normalization', float(len(ran_north))/len(elg_north), float(len(ran_south))/len(elg_south))
    #I have to subsample ran_north such float(len(ran_north))/len(lop_north) = 1.4975425785620076
    #Therefore len(ran_north) = 1.4975425785620076*len(lop_north)
    #print(ran_north.dtype)

    #print(ran_north, ran_south)

    random_file = vstack([ran_north, ran_south])

    sel_ran = return_north(random_file['RA'], random_file['DEC'])

    ran_north = random_file[sel_ran]
    ran_south = random_file[~sel_ran]
    sizerans_north = len(ran_north)
    sizerans_south = len(ran_south)

    '''

    lop = Table.read(f"/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/AbacusSummit_base_c000_ph{lop_rea}/CutSky/ELG_v5/{lop_s}/forclustering/cutsky_abacusHF_DR2_ELG_LOP_{lop_s.replace('.','p')}_zcut_0p8to1p6_clustering.dat.fits")

    ra_lop, dec_lop = lop['RA'], lop['DEC']
    sel_north = return_north(ra_lop, dec_lop)

    lop_north = lop[sel_north]
    lop_south = lop[~sel_north]

    print('size data', len(lop_north), len(lop_south))



    rng_N = Generator(PCG64())
    rng_S = Generator(PCG64())
    inds_z_N = rng_N.choice(len(lop_north), sizerans_north)
    inds_z_S = rng_S.choice(len(lop_south), sizerans_south)

    #ran_north['Z_ELG_LOP'] = lop_north['Z'][inds_z_N]
    ran_north['Z_ELG_LOP'] = lop_north['Z'][inds_z_N]
    ran_south['Z_ELG_LOP'] = lop_south['Z'][inds_z_S]
    #ran_south['Z_ELG_LOP'] = lop_south['Z'][inds_z_S]

    print('size randoms', len(ran_north), len(ran_south))
    print('randoms to data', float(len(ran_north))/len(lop_north), float(len(ran_south))/len(lop_south))


    newgen = Generator(PCG64())
    print(ran_north.dtype)
    N = len(ran_north)

    inds = newgen.choice(np.arange(N), int(len(lop_north)*float(len(ran_south))/len(lop_south)), replace=False)
    mask = np.zeros(N, dtype=bool)
    mask[inds] = True
    
    ran_north["ELG_LOP_MASK"] = mask

    N = len(ran_south)
    mask = np.ones(N, dtype=bool)
    ran_south["ELG_LOP_MASK"] = mask

    #ran_north = Table(inds)

    print('size randoms after normalization', len(ran_north), len(ran_south))
    print('randoms to data after normalization', float(len(ran_north))/len(lop_north), float(len(ran_south))/len(lop_south))
    #I have to subsample ran_north such float(len(ran_north))/len(lop_north) = 1.4975425785620076
    #Therefore len(ran_north) = 1.4975425785620076*len(lop_north)
    #print(ran_north.dtype)

 
    random_file = vstack([ran_north, ran_south])

    #ranf = f'/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/randoms/rands_intiles_DARK_{i}_v2.fits'
    ranf = f'/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0/randoms/rands_intiles_DARK_{i}_NO_imagingmask_withz.fits'
    
    common.write_LSS_scratchcp(random_file, ranf, extname = 'RANDOMS')
    fits.setval(ranf, 'ID_LRG', value = str(lrg_rea), ext = 1)
    #fits.setval(ranf, 'ID_ELG', value = str(elg_rea), ext = 1)
    fits.setval(ranf, 'ID_ELG_LOP', value = str(lop_rea), ext = 1)
    fits.setval(ranf, 'ID_QSO', value = str(qso_rea), ext = 1)

    fits.setval(ranf, 'LRG_SNAP', value = str(lrg_s), ext = 1)
    #fits.setval(ranf, 'ELG_SNAP', value = str(elg_s), ext = 1)
    fits.setval(ranf, 'ELG_LOP_SNAP', value = str(lop_s), ext = 1)
    fits.setval(ranf, 'QSO_SNAP', value = 'z.1400', ext = 1)

    


