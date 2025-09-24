#make sure add the LSS repo to your python path
from astropy.io import fits # Access to FITS (Flexible Image Transport System) files.
from astropy.table import Table, hstack, vstack, Column # A class to represent tables of heterogeneous data.
import numpy as np
import os, sys
import argparse
import random
from numpy.random import Generator, PCG64
import LSS.common_tools as common
from LSS.globals import main
import errno

def create_dirs(value):
    if not os.path.exists(value):
        try:
            os.makedirs(value, 0o755)
            print('made ' + value)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise



parser = argparse.ArgumentParser()
parser.add_argument("--input_mockpath", help = "full directory path to input mocks", default = '/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v1.0')
parser.add_argument("--input_mockfile", help = "mock file name", default = 'forFA0_Y3_noimagingmask_applied.fits')
parser.add_argument("--output_fullpathfn", help = "output mock file and full path", default = None)
parser.add_argument("--nrans", help = "number of randoms", default = None)
parser.add_argument("--tracer", help = "LRG, ELG, ELG_LOP or QSO", default = 'LRG')
parser.add_argument("--snapshot", help = "LRG: z0.500, z0.725, z0.950. ELG: z0.950, z1.175, z1.475. QSO: z1.400", default = 'z0.500')
parser.add_argument("--abacus_realization", help = "0,1...", type=int, default = 0)
parser.add_argument("--ran_seed", help = "seed for randoms; make sure this is different if running many in parallel", default = None)
parser.add_argument("--redshift_range", help = "redshift range, as separated by coma (0.5,1.4)", default = None)

args = parser.parse_args()

if args.ran_seed is not None:
    random.seed(int(args.ran_seed))

abacus_realization = str(int(args.abacus_realization)).zfill(3)

if args.redshift_range is not None:
    redshift_range = [float(x) for x in args.redshift_range.split(",")]
else:
    redshift_range = None


if args.tracer == 'ELG' or args.tracer == 'ELG_LOP': 
    tracer = 'ELG_v5'
else:
    tracer = args.tracer

inputh_path = os.path.join(args.input_mockpath, f'AbacusSummit_base_c000_ph{abacus_realization}/CutSky/{tracer}/{args.snapshot}')

if args.output_fullpathfn is None:
    output_path = os.path.join(inputh_path, 'forclustering')
else:
    output_path = os.path.join(args.output_fullpathfn, f'AbacusSummit_base_c000_ph{abacus_realization}/CutSky/{tracer}/{args.snapshot}/forclustering')

#if os.path.isdir(output_path):
#    print(output_path, 'exist')
#    exit()

create_dirs(output_path)

if args.nrans is not None:
    random_numbers = random.sample(range(18), int(args.nrans))  # range(18) gives numbers from 0 to 17
    print(random_numbers)
    ranlist = ','.join(str(x) for x in random_numbers)
    print(ranlist)

def mask_abacusHF(nz = 0, foot = None, nz_lop = 0):
    if foot == 'Y1':
        Y5 = 0
        Y1 = 1
        Y3 = 0
    elif foot == 'Y3':
        Y5 = 0
        Y1 = 0
        Y3 = 1
    else:
        Y5 = 1
        Y1 = 0
        Y3 = 0

    return nz * (2**0) + Y5 * (2**1) + nz_lop * (2**2) + Y1 * (2**3) + Y3 * (2**5) 


data = Table.read(os.path.join(inputh_path, args.input_mockfile))

# Adding the WEIGHT column
if 'WEIGHT' not in data.colnames:
    data['WEIGHT'] = np.ones(data['RA'].shape[0])


type_ = args.tracer
if type_ == 'ELG_LOP':
    byte_selection = mask_abacusHF(nz = 1, foot = 'Y3', nz_lop = 1)
else:
    byte_selection = mask_abacusHF(foot = 'Y3')
    #TEMPbyte_selection = mask_abacusHF(nz = 1, foot = 'Y3')

mask_selection = (data['STATUS'] & byte_selection) == byte_selection

if redshift_range is not None:
    mask_selection &= (data['RSDZ'] >= redshift_range[0]) & (data['RSDZ'] <= redshift_range[1])
    redshift_tag = 'zcut_' + args.redshift_range.replace('.', 'p').replace(',', 'to')
else:
    redshift_tag = 'nozcut'

targets = data[mask_selection]

snap = args.snapshot.replace('.','p')
name_output = f'cutsky_allzs_abacusHF_DR2_{type_}_{snap}_{redshift_tag}_clustering.dat.fits'
output_file = os.path.join(output_path, name_output)
print('will save to', output_file)

targets.rename_column('RSDZ', 'Z')

common.write_LSS_scratchcp(targets, output_file, extname = 'TARGETS')
fits.setval(output_file, 'EXTNAME', value = 'TARGETS', ext = 1)
fits.setval(output_file, 'OBSCON', value = 'DARK', ext = 1)




'''
if args.nrans is not None:
    fits.setval(output_file, 'RANDOMS', value = ranlist, ext = 1)
    fits.setval(output_file, 'RANSEED', value = str(args.ran_seed), ext = 1)

    ranf = '/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/rands_intiles_DARK_nomask_{rannum}.fits'

    Z = targets['Z']
    weight = targets['WEIGHT']


    for j, rannum in enumerate(random_numbers):
        print('using random', rannum)
        temp_ran = ranf.format(rannum = rannum)
        output_ran = os.path.join(output_path, f'cutsky_abacusHF_DR2_{type_}_{snap}_{redshift_tag}_{j}_clustering.ran.fits')
        random_file = Table.read(temp_ran)

        rng = Generator(PCG64())
        sizerans = len(random_file)
        inds_z = rng.choice(len(Z), sizerans)
        random_file['Z'] = Z[inds_z]
        random_file['WEIGHT'] = weight[inds_z]
        common.write_LSS_scratchcp(random_file, output_ran, extname = 'RANDOMS')
        fits.setval(output_ran, 'RANID', value = str(rannum), ext = 1)

'''
