from astropy.io import fits # Access to FITS (Flexible Image Transport System) files.
from astropy.table import Table, hstack, vstack, Column # A class to represent tables of heterogeneous data.
import fitsio
import numpy as np
import glob
import os
import h5py
import argparse
import sys

from desitarget.targetmask import obsconditions

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
parser.add_argument("--realmin", help="number for the realization",default=1,type=int)
parser.add_argument("--realmax", help="number for the realization",default=2,type=int)
parser.add_argument("--survey", help="points to set of tiles",default='DA02')
parser.add_argument("--prog", help="dark or bright",default='dark')
parser.add_argument("--base_output", help="base directory for output",default='/global/cfs/cdirs/desi/survey/catalogs/main/mocks/')
parser.add_argument("--prep", help="prepare file for fiberassign?",default='y')
parser.add_argument("--runfa", help="run fiberassign",default='n')
parser.add_argument("--par", help="running in parallel?",default='n')

args = parser.parse_args()

if args.prog == 'dark':
    types = ['ELG', 'LRG', 'QSO']
    desitar = {'ELG':34,'LRG':1,'QSO':4}
    priority = {'ELG':3000,'LRG':3200,'QSO':3400}

for real in range(args.realmin,args.realmax):
    if not (args.mockver is None):
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
        else:
            raise ValueError(args.mockver+' not supported with legacy mockver argument. Use mockpath/mockfilename arguments instead.')
    else:
        mockpath = args.mockpath
        file_name = args.mockfile
        out_file_name = args.base_output + '/forFA_Real{0}.fits'.format(real)
    
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
            thepath = os.path.join(mockpath, type_, zs[type_], file_name.format(TYPE = type_, Z = zs[type_], PH = "%03d" % real))
            #f = fits.open(thepath)
            print('thepath')
            print(thepath)
            data = fitsio.read(thepath,columns=['RA','DEC','Z','Z_COSMO','STATUS'])#f[1].data
            print(data.dtype.names)
            print(type_,len(data))
            status = data['STATUS'][()]
            idx = np.arange(len(status))
            mask_main = mask(main=0, nz=1, Y5=1, sv3=0)
            if type_ == 'LRG':
                mask_main = mask(main=1, nz=1, Y5=1, sv3=0)
            idx_main = idx[(status & (mask_main))==mask_main]
            data = data[idx_main]
            print(len(data))
            data = Table(data)
            data['DESI_TARGET'] = desitar[type_]
            data['PRIORITY_INIT'] = priority[type_]
            data['PRIORITY'] = priority[type_]
            datat.append(data)
        targets = vstack(datat)
        print(len(targets))
        del datat

    

    if args.prep == 'y':
        n=len(targets)
        targets.rename_column('Z_COSMO', 'TRUEZ') 
        targets.rename_column('Z', 'RSDZ') 
        targets['BGS_TARGET'] = np.zeros(n, dtype='i8')
        targets['MWS_TARGET'] = np.zeros(n, dtype='i8')
        targets['SUBPRIORITY'] = np.random.uniform(0, 1, n)
        targets['BRICKNAME'] = np.full(n, '000p0000')    #- required !?!
        targets['OBSCONDITIONS'] = obsconditions.mask(program.upper()) #np.zeros(n, dtype='i8')+int(3) 
        targets['NUMOBS_MORE'] = np.zeros(n, dtype='i8')+int(1) 
        targets['NUMOBS_INIT'] = np.zeros(n, dtype='i8')+int(1)
        targets['SCND_TARGET'] = np.zeros(n, dtype='i8')+int(0)
        targets['ZWARN'] = np.zeros(n, dtype='i8')+int(0)
        targets['TARGETID'] = np.arange(1,n+1)

        targets.write(out_file_name, overwrite = True)

        fits.setval(out_file_name, 'EXTNAME', value='TARGETS', ext=1)
        fits.setval(out_file_name, 'OBSCON', value=args.prog.upper(), ext=1)

    if args.runfa == 'y':
        from LSS.mocktools import get_fba_mock
        script_fn = get_fba_mock(mockdir,real,survey=args.survey,prog=args.prog)
        os.system('chmod +x '+script_fn)
        
        os.system(script_fn)



sys.exit()

