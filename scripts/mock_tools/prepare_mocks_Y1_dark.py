from astropy.io import fits # Access to FITS (Flexible Image Transport System) files.
from astropy.table import Table, hstack, vstack, Column # A class to represent tables of heterogeneous data.
import fitsio
import numpy as np
import os
import argparse
import sys
import json
from desitarget.targetmask import obsconditions
from desimodel.footprint import is_point_in_desi

import LSS.common_tools as common
from LSS.imaging import get_pixel_bitmasknobs as bitmask #get_nobsandmask
from LSS.main.cattools import count_tiles_better
from LSS.globals import main


def create_dir(value):
    if not os.path.exists(value):
        try:
            os.makedirs(value, 0o755)
            print('Check directories', value)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

def mask_firstgen(main=0, nz=0, Y5=0, sv3=0):
    return main * (2**3) + sv3 * (2**2) + Y5 * (2**1) + nz * (2**0)

def mask_secondgen(nz=0, foot=None, nz_lop=0):
    if foot == 'Y1':
        Y5 = 0
        Y1 = 1
    elif foot == 'Y5':
        Y5 = 1
        Y1 = 0
    else:
        Y5 = 0
        Y1 = 0
    return nz * (2**0) + Y5 * (2**1) + nz_lop * (2**2) + Y1 * (2**3)




if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 

parser = argparse.ArgumentParser()
parser.add_argument("--mockver", help="type of mock to use",default=None)
parser.add_argument("--mockpath", help="Location of mock file(s)",default='/global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky')
parser.add_argument("--mockfile", help="formattable name of mock file(s). e.g. cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits. TYPE will be replaced with tracer type. PH will be replaced with realization number for simulation of mock.",default='cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits')
parser.add_argument("--realmin", help="number for the realization",default=0,type=int)
parser.add_argument("--realmax", help="number for the realization",default=1,type=int)
parser.add_argument("--prog", help="dark or bright",default='dark')
parser.add_argument("--base_output", help="base directory for output",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/')
parser.add_argument("--apply_mask", help="apply the same mask as applied to desi targets?",default='y')
parser.add_argument("--downsampling", help="downsample to Y1 target density in SecondGen Abacus mocks?",default='n')
parser.add_argument("--isProduction", help="Say yes if you want to save in main production directory",default='n')
parser.add_argument("--overwrite", help="Overwrite. if it is in production, this always will be no. You must delete by hand first", default=0, type=bool)
parser.add_argument("--split_snapshot", help="apply different snapshots to different redshift ranges?",default='n')
parser.add_argument("--new_version", help="If production, and this is a new version, set to name, for example, AbacusSummit_v3",default=None)

args = parser.parse_args()
tiletab = Table.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/tiles-{PROG}.fits'.format(PROG = args.prog.upper()))

if args.prog == 'dark':
    types = ['ELG', 'LRG', 'QSO']
    priority = {'ELG':3000, 'LRG':3200, 'QSO':3400}
    mainp = main(tp = 'QSO', specver = 'iron')
    desitar = {'ELG':34, 'LRG':1, 'QSO':4}
    numobs = {'ELG':2, 'LRG':2, 'QSO':4}
    
    if args.split_snapshot == 'y':
        zs = {'ELG':{'z0.950':[0.,1.1], 'z1.325':[1.1,99.]}, 'LRG':{'z0.500':[0.,0.6], 'z0.800':[0.6,99.]}, 'QSO':{'z1.400':[0.,99.]}}
    else:
        zs = {'ELG':'z1.100', 'LRG':'z0.800', 'QSO':'z1.400'}

    
    if args.mockver == 'ab_secondgen':
        desitar = {'ELG':2**1, 'LRG':2**0, 'QSO':2**2}
        downsampling = {'ELG':0.7345658717688022, 'LRG':0.708798313382828, 'QSO':0.39728966594530174}
        percentage_elg_hip = 0.1

Abacus_dir = 'AbacusSummit'

if args.isProduction == 'y':
    args.base_output = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks'
    args.overwrite = False
    if args.new_version is not None:
        Abacus_dir = args.new_version
else:
    if args.base_output == '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks' or args.base_output == '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/':
        args.base_output = os.environ[scratch] 
        print('This is not production, run on user scratch', os.environ[scratch])
    else:
        print('Saving to path', args.base_output)
    if args.new_version is not None:
        Abacus_dir = args.new_version



for real in range(args.realmin, args.realmax):
    if not (args.mockver is None):
        if args.mockver == 'ab_firstgen':
            mockpath = '/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/'
            file_name = 'cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits'
            mockdir = os.path.join(args.base_output, 'FirstGenMocks', 'AbacusSummit')
            
            out_file_name = os.path.join(mockdir, 'forFA{real}.fits'.format(real=real))
            
        
        elif args.mockver == 'ezmocks6':
            mockdir = os.path.join(args.base_output, 'EZMocks_6Gpc')
            mockpath = '/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc'
            out_file_name = os.path.join(mockdir, 'EZMocks_6Gpc_{real}.fits'.format(real=real))

        elif args.mockver == 'ab_secondgen':
            mockpath = args.mockpath
            file_name = 'cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits'
            
            mockdir = os.path.join(args.base_output, 'SecondGenMocks', Abacus_dir)
            #if args.split_snapshot == 'y':
            create_dir(mockdir)
            if not os.path.isfile(os.path.join(mockdir, 'prepare_mock_arguments.txt')):
                with open(os.path.join(mockdir, 'prepare_mock_arguments.txt'), 'w') as f:
                    json.dump(args.__dict__, f, indent=2)

            out_file_name = os.path.join(mockdir, 'forFA{real}.fits'.format(real=real))
            #else:
            #    out_file_name = os.path.join(mockdir, 'forFA{real}.fits'.format(real=real))

        else:
            raise ValueError(args.mockver+' not supported with legacy mockver argument. Use mockpath/mockfilename arguments instead.')
    else:
        mockpath = args.mockpath
        file_name = args.mockfile
        mockdir = args.base_output
        out_file_name = os.path.join(mockdir, 'forFA{0}.fits'.format(real))
        print('generic mock, it needs a mock generation to continue, it will select mockver = ab_secondgen')
        args.mockver = 'ab_secondgen'

    
    print('testing and creating output directory', mockdir)
    create_dir(mockdir)
    print('will write outputs to ', out_file_name)

    
    mockdir = args.base_output


    datat = []
    for type_ in types:
        if args.mockver == 'ab_firstgen' or args.mockver == 'ab_secondgen':
            if args.split_snapshot == 'y':
                datas = []

                for bins in zs[type_]:
                    print(bins)
                    thepath = os.path.join(mockpath, type_, bins, file_name.format(TYPE = type_, Z = bins, PH = "%03d" % real))
                    print('thepath')
                    print(thepath)
                    dat = fitsio.read(thepath, columns=['RA','DEC','Z','Z_COSMO','STATUS'])#f[1].data
                    mask = (dat['Z']>= zs[type_][bins][0])&(dat['Z']< zs[type_][bins][1])
                    datas.append(Table(dat[mask]))
                data = vstack(datas)
                del datas
                del dat
            else:
                thepath = os.path.join(mockpath, type_, zs[type_], file_name.format(TYPE = type_, Z = zs[type_], PH = "%03d" % real))
                data = fitsio.read(thepath, columns=['RA','DEC','Z','Z_COSMO','STATUS'])#f[1].data
        
        elif args.mockver == 'ezmocks6':
            path_ezmock = os.path.join(mockpath, type_, zs[type_])
            if  type_ == "LRG":
                infn1 = os.path.join(path_ezmock, "cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed{real}_NGC.fits".format(real = real))
                infn2 = os.path.join(path_ezmock, "cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed{real}_SGC.fits".format(real = real))
            elif type_ == "ELG":
                infn1 = os.path.join(path_ezmock, "cutsky_ELG_z1.100_EZmock_B6000G1536Z1.1N648012690_b0.345d1.45r40c0.05_seed{real}_NGC.fits".format(real = real))
                infn2 = os.path.join(path_ezmock, "cutsky_ELG_z1.100_EZmock_B6000G1536Z1.1N648012690_b0.345d1.45r40c0.05_seed{real}_SGC.fits".format(real = real))
            elif type_ == "QSO":
                infn1 = os.path.join(path_ezmock, "cutsky_QSO_z1.400_EZmock_B6000G1536Z1.4N27395172_b0.053d1.13r0c0.6_seed{real}_NGC.fits".format(real = real))
                infn2 = os.path.join(path_ezmock, "cutsky_QSO_z1.400_EZmock_B6000G1536Z1.4N27395172_b0.053d1.13r0c0.6_seed{real}_SGC.fits".format(real = real))
            tars1 = Table.read(infn1)
            tars2 = Table.read(infn2)
            tars1["GALCAP"] = "N"
            tars2["GALCAP"] = "S"
            data = vstack([tars1, tars2])


        print(data.dtype.names)
        print(type_, len(data))
        status = data['STATUS'][()]
        idx = np.arange(len(status))

        if args.mockver == 'ab_secondgen':

            mask_main = mask_secondgen(nz=1, foot='Y1')
            idx_main = idx[(status & (mask_main))==mask_main]

            if type_ == 'LRG' or type_ == 'QSO':
                if args.downsampling == 'y':
                    ran_tot = np.random.uniform(size = len(idx_main))
                    idx_main = idx_main[(ran_tot<=downsampling[type_])]
                data = data[idx_main]
                data = Table(data)
                
                data['DESI_TARGET'] = desitar[type_]
                data['PRIORITY_INIT'] = priority[type_]
                data['PRIORITY'] = priority[type_]
                data['NUMOBS_MORE'] = numobs[type_]
                data['NUMOBS_INIT'] = numobs[type_]
                datat.append(data)

            else:

                mask_LOP = mask_secondgen(nz=1, foot='Y1', nz_lop=1)
                idx_LOP = idx[(status & (mask_LOP))==mask_LOP]


                idx_VLO = np.setdiff1d(idx_main, idx_LOP)

                if args.downsampling == 'y':
                    ran_lop = np.random.uniform(size = len(idx_LOP))
                    idx_LOP = idx_LOP[(ran_lop<=downsampling[type_])]
                    ran_vlo = np.random.uniform(size = len(idx_VLO))
                    idx_VLO = idx_VLO[(ran_vlo<=downsampling[type_])]

                data_lop = Table(data[idx_LOP])
                data_vlo = Table(data[idx_VLO])
                
                df_lop=data_lop.to_pandas()
                df_vlo=data_vlo.to_pandas()

                num_HIP_LOP = int(len(df_lop) * percentage_elg_hip)
                df_HIP_LOP = df_lop.sample(n=num_HIP_LOP)
                remaining_LOP = df_lop.drop(df_HIP_LOP.index)
                df_HIP_LOP.reset_index(drop=True, inplace=True)
                remaining_LOP.reset_index(drop=True, inplace=True)

                num_HIP_VLO = int(len(df_vlo) * percentage_elg_hip)
                df_HIP_VLO = df_vlo.sample(n=num_HIP_VLO)
                remaining_VLO = df_vlo.drop(df_HIP_VLO.index)
                df_HIP_VLO.reset_index(drop=True, inplace=True)
                remaining_VLO.reset_index(drop=True, inplace=True)

                remaining_LOP['PRIORITY_INIT'] = 3100
                remaining_LOP['PRIORITY'] = 3100
                remaining_LOP['DESI_TARGET'] = 2**5 + 2**1


                remaining_VLO['PRIORITY_INIT'] = 3000
                remaining_VLO['PRIORITY'] = 3000
                remaining_VLO['DESI_TARGET'] = 2**7 + 2**1

                df_HIP_LOP['PRIORITY_INIT'] = 3200
                df_HIP_LOP['PRIORITY'] = 3200
                df_HIP_LOP['DESI_TARGET'] = 2**6 + 2**1 + 2**5

                df_HIP_VLO['PRIORITY_INIT'] = 3200
                df_HIP_VLO['PRIORITY'] = 3200
                df_HIP_VLO['DESI_TARGET'] = 2**6 + 2**1 + 2**7

                remaining_LOP['NUMOBS_MORE'] = numobs[type_]
                remaining_LOP['NUMOBS_INIT'] = numobs[type_]
                remaining_VLO['NUMOBS_MORE'] = numobs[type_]
                remaining_VLO['NUMOBS_INIT'] = numobs[type_]
                df_HIP_LOP['NUMOBS_MORE'] = numobs[type_]
                df_HIP_LOP['NUMOBS_INIT'] = numobs[type_]
                df_HIP_VLO['NUMOBS_MORE'] = numobs[type_]
                df_HIP_VLO['NUMOBS_INIT'] = numobs[type_]

                datat.append(Table.from_pandas(remaining_LOP))
                datat.append(Table.from_pandas(remaining_VLO))
                datat.append(Table.from_pandas(df_HIP_LOP))
                datat.append(Table.from_pandas(df_HIP_VLO))

        else:
            mask_main = mask_firstgen(main=0, nz=1, Y5=0, sv3=0) #no longer cutting to Y5 footprint because it doesn't actually cover Y1
            if type_ == 'LRG':
                mask_main = mask_firstgen(main=1, nz=1, Y5=0, sv3=0)
            idx_main = idx[(status & (mask_main))==mask_main]
            data = data[idx_main]
            print(len(data))
            data = Table(data)
            data['DESI_TARGET'] = desitar[type_]
            data['PRIORITY_INIT'] = priority[type_]
            data['PRIORITY'] = priority[type_]
            data['NUMOBS_MORE'] = numobs[type_] 
            data['NUMOBS_INIT'] = numobs[type_]

            datat.append(data)
    
    targets = vstack(datat)
    del datat
    if args.mockver != 'ab_secondgen':
        print(len(targets),' in Y5 area')
        selY1 = is_point_in_desi(tiletab,targets['RA'],targets['DEC'])
        targets = targets[selY1]
    print(len(targets),' in Y1 area')

    if args.apply_mask == 'y':
        print('getting nobs and mask bits')
        mask = bitmask.get_nobsandmask(targets)
        maskv = mask.get_nobsandmask()
        maskcols = ['NOBS_G','NOBS_R','NOBS_Z','MASKBITS']
        for col in maskcols:
            targets[col] = maskv[col]
        del maskv
        targets = common.cutphotmask(targets, bits=mainp.imbits)
        


    n=len(targets)
    targets.rename_column('Z_COSMO', 'TRUEZ') 
    targets.rename_column('Z', 'RSDZ') 
    targets['BGS_TARGET'] = np.zeros(n, dtype='i8')
    targets['MWS_TARGET'] = np.zeros(n, dtype='i8')
    targets['SUBPRIORITY'] = np.random.uniform(0, 1, n)
    targets['BRICKNAME'] = np.full(n, '000p0000')    #- required !?!
    targets['OBSCONDITIONS'] = obsconditions.mask(args.prog.upper()) #np.zeros(n, dtype='i8')+int(3) 
    targets['SCND_TARGET'] = np.zeros(n, dtype='i8')+int(0)
    targets['ZWARN'] = np.zeros(n, dtype='i8')+int(0)
    targets['TARGETID'] = np.random.permutation(np.arange(1,n+1))

    targets.write(out_file_name, overwrite = args.overwrite)

    fits.setval(out_file_name, 'EXTNAME', value='TARGETS', ext=1)
    fits.setval(out_file_name, 'OBSCON', value=args.prog.upper(), ext=1)




sys.exit()

