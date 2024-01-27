from astropy.io import fits # Access to FITS (Flexible Image Transport System) files.
from astropy.table import Table, hstack, vstack, Column # A class to represent tables of heterogeneous data.
import fitsio
import numpy as np
import glob
import os
import h5py
import argparse
import sys
import pickle
from desitarget.targetmask import obsconditions
from desimodel.footprint import is_point_in_desi

import LSS.common_tools as common
from LSS.imaging import get_pixel_bitmasknobs as bitmask #get_nobsandmask
from LSS.main.cattools import count_tiles_better
from LSS.globals import main


if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 

parser = argparse.ArgumentParser()
parser.add_argument("--mockver", help="type of mock to use",default=None)
parser.add_argument("--mockpath", help="Location of mock file(s)",default='/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/')
parser.add_argument("--mockfile", help="formattable name of mock file(s). e.g. cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits. TYPE will be replaced with tracer type. PH will be replaced with realization number for simulation of mock.",default='cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits')
#parser.add_argument("--realization", help="number for the realization",default=1,type=int)
parser.add_argument("--tracer", help="if only running on one tracer, specify here. Default is None which will select all tracers for a given program.", default=None, type=str)
parser.add_argument("--realmin", help="number for the realization",default=1,type=int)
parser.add_argument("--realmax", help="number for the realization",default=2,type=int)
parser.add_argument("--prog", help="dark or bright",default='dark')
parser.add_argument("--base_output", help="base directory for output",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/')
parser.add_argument("--prep", help="prepare file for fiberassign?",default='y')
parser.add_argument("--apply_mask", help="apply the same mask as applied to desi targets?",default='y')
parser.add_argument("--par", help="running in parallel?",default='n')


args = parser.parse_args()

tiletab = Table.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/tiles-'+args.prog.upper()+'.fits')

if args.prog == 'dark':
    if args.tracer is None:
        types = ['ELG', 'LRG', 'QSO']
    else:
        types = [args.tracer]
    desitar = {'ELG':34,'LRG':1,'QSO':4}
    priority = {'ELG':3000,'LRG':3200,'QSO':3400}
    mainp = main(tp='QSO',specver='iron')

for real in range(args.realmin,args.realmax):
    if args.mockver == 'ab_firstgen':
        mockpath = '/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/AbacusSummit/CutSky/'
    
        file_name = 'cutsky_{TYPE}_{Z}_AbacusSummit_base_c000_ph{PH}.fits'
        out_file_name = args.base_output+'/FirstGenMocks/AbacusSummit/forFA'+str(real)+'.fits'
        if not os.path.exists(args.base_output+'/FirstGenMocks'):
            os.mkdir(args.base_output+'/FirstGenMocks')
            print('made '+args.base_output+'/FirstGenMocks')
        if not os.path.exists(args.base_output+'/FirstGenMocks/AbacusSummit'):
            os.mkdir(args.base_output+'/FirstGenMocks/AbacusSummit')
            print('made '+args.base_output+'/FirstGenMocks/AbacusSummit')
        mockdir = args.base_output+'/FirstGenMocks/AbacusSummit/'
    
    if args.mockver == 'ezmocks6':
        out_file_name = args.base_output + '/EZMocks_6Gpc/EZMocks_6Gpc_' + str(real) + '.fits'
        if not os.path.exists(args.base_output + '/EZMocks_6Gpc'):
            os.makedirs(args.base_output + '/EZMocks_6Gpc')
            print('made ' + args.base_output + '/EZMocks_6Gpc')
        mockdir = args.base_output + '/EZMocks_6Gpc/'

    elif args.mockver == 'glam':
        if args.mockfile is None:
            file_name = 'lightcone_galaxies_{TYPE}_{PH}.pickle'
        else:
            file_name = args.mockfile
        mockpath = args.mockpath
        #file_name = args.mockfile
        if args.tracer is None:
            out_file_name = args.base_output + '/forFA{0}.fits'.format(real)
        else:
            args.base_output += '/{0}/'.format(args.base_output)
            out_file_name = args.base_output +  '/forFA{0}.fits'.format(real)
    else:
        mockpath = args.mockpath
        file_name = args.mockfile
        out_file_name = args.base_output + '/forFA{0}.fits'.format(real)
    
    print('will write to '+out_file_name)
    if not os.path.exists(args.base_output):
        os.makedirs(args.base_output)
        print('made '+args.base_output)
    
    mockdir = args.base_output
    zs = {'ELG':'z1.100','LRG':'z0.800','QSO':'z1.400'}


    def mask(main=0, nz=0, Y5=0, sv3=0):
        return main * (2**3) + sv3 * (2**2) + Y5 * (2**1) + nz * (2**0)
    if args.prep == 'y':
        datat = []
        for type_ in types:
            if args.mockver == 'ab_firstgen':
                thepath = os.path.join(mockpath, type_, zs[type_], file_name.format(TYPE = type_, Z = zs[type_], PH = "%03d" % real))
                print('thepath')
                print(thepath)
                data = fitsio.read(thepath,columns=['RA','DEC','Z','Z_COSMO','STATUS'])#f[1].data
            elif args.mockver == 'ezmocks6':
                    if  type_ == "LRG":
                        infn1 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed%s_NGC.fits"%real
                        infn2 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed%s_SGC.fits"%real
                    elif type_ == "ELG":
                        infn1 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/ELG/z1.100/cutsky_ELG_z1.100_EZmock_B6000G1536Z1.1N648012690_b0.345d1.45r40c0.05_seed%s_NGC.fits"%real
                        infn2 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/ELG/z1.100/cutsky_ELG_z1.100_EZmock_B6000G1536Z1.1N648012690_b0.345d1.45r40c0.05_seed%s_SGC.fits"%real
                    elif type_ == "QSO":
                        infn1 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/QSO/z1.400/cutsky_QSO_z1.400_EZmock_B6000G1536Z1.4N27395172_b0.053d1.13r0c0.6_seed%s_NGC.fits"%real
                        infn2 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/QSO/z1.400/cutsky_QSO_z1.400_EZmock_B6000G1536Z1.4N27395172_b0.053d1.13r0c0.6_seed%s_SGC.fits"%real
                # infn1 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed1_NGC.fits"
                # infn2 = "/global/cfs/cdirs/desi/cosmosim/FirstGenMocks/EZmock/CutSky_6Gpc/LRG/z0.800/cutsky_LRG_z0.800_EZmock_B6000G1536Z0.8N216424548_b0.385d4r169c0.3_seed1_SGC.fits"
                    tars1 = Table.read(infn1)#fitsio.read(infn1)
                    tars2 = Table.read(infn2)#fitsio.read(infn2)
                    tars1["GALCAP"] = "N"
                    tars2["GALCAP"] = "S"
                    tars = vstack([tars1, tars2])
                    data = tars
                    #tars['TARGETID'] = np.arange(len(tars))
            elif args.mockver.lower() == 'glam':
                #general GLAM filename: lightcone_galaxies_{TRACER}_{str(REAL).zfill(4)}.pickle
                #thepath=os.path.join(mockpath, type_, file_name.format(TYPE = type_, PH = "%04d" % real))
                file_name = 'lightcone_galaxies_{TYPE}_{PH}.pickle'
                thepath=os.path.join(mockpath, file_name.format(TYPE = type_, PH = "%04d" % real))

                temp = pickle.load(open(thepath, 'rb'), fix_imports = True)
                #'xh,yh,zh,vx,vy,vz,mh,rvirh,rsh,TARGETID,Z_COSMO,Z,RA,DEC'
                #gal_id (-> TARGETID), z_cos (-> Z_COSMO), z_obs (-> Z), ra (-> RA), dec (-> DEC)
                data = Table()
                #data['TARGETID'] = temp['gal_id']
                data['Z_COSMO'] = temp['z_cos']
                data['Z'] = temp['z_obs']
                data['RA'] = temp['ra']
                data['DEC'] = temp['dec']
                print('add GALCAP Values here if needed')
                data['STATUS'] = 3*np.ones(len(data), dtype = int)*is_point_in_desi(tiletab, data['RA'], data['DEC'])
                del temp

            else:
                thepath=os.path.join(mockpath, type_, file_name.format(TYPE = type_, PH = "%04d" % real))
                if thepath.endswith('fits') or thepath.endswith('fit') or thepath.endswith('fits.gz') or thepath.endswith('fit.gz'):
                    data = fitsio.read(thepath)
                elif thepath.endswith('pickle'):
                    data = pickle.load(open(thepath, 'rb'), fix_imports = True)
                else:
                    raise NotImplementedError('Only attempting to load generic mock files in fits and pickle formats.')

            #f = fits.open(thepath)
            print(data.dtype.names)
            print(type_,len(data))
            status = data['STATUS'][()]
            idx = np.arange(len(status))

            #def mask(main=0, nz=0, Y5=0, sv3=0):
            #    return main * (2**3) + sv3 * (2**2) + Y5 * (2**1) + nz * (2**0)
            mask_main = mask(main=0, nz=1, Y5=0, sv3=0) #no longer cutting to Y5 footprint because it doesn't actually cover Y1
            if type_ == 'LRG':
                mask_main = mask(main=1, nz=1, Y5=0, sv3=0)
            idx_main = idx[(status & (mask_main))==mask_main]
            data = data[idx_main]
            print(len(data))
            data = Table(data)
            data['DESI_TARGET'] = desitar[type_]
            data['PRIORITY_INIT'] = priority[type_]
            data['PRIORITY'] = priority[type_]
            datat.append(data)
        targets = vstack(datat)
        print(len(targets),' in Y5 area')
        del datat
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
        targets = common.cutphotmask(targets,bits=mainp.imbits)
        

    if args.prep == 'y':
        n=len(targets)
        targets.rename_column('Z_COSMO', 'TRUEZ') 
        targets.rename_column('Z', 'RSDZ') 
        targets['BGS_TARGET'] = np.zeros(n, dtype='i8')
        targets['MWS_TARGET'] = np.zeros(n, dtype='i8')
        targets['SUBPRIORITY'] = np.random.uniform(0, 1, n)
        targets['BRICKNAME'] = np.full(n, '000p0000')    #- required !?!
        targets['OBSCONDITIONS'] = obsconditions.mask(args.prog.upper()) #np.zeros(n, dtype='i8')+int(3) 
        targets['NUMOBS_MORE'] = np.zeros(n, dtype='i8')+int(1) 
        targets['NUMOBS_INIT'] = np.zeros(n, dtype='i8')+int(1)
        targets['SCND_TARGET'] = np.zeros(n, dtype='i8')+int(0)
        targets['ZWARN'] = np.zeros(n, dtype='i8')+int(0)
        targets['TARGETID'] = np.arange(1,n+1)

        targets.write(out_file_name, overwrite = True)

        fits.setval(out_file_name, 'EXTNAME', value='TARGETS', ext=1)
        fits.setval(out_file_name, 'OBSCON', value=args.prog.upper(), ext=1)




sys.exit()

