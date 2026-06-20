from astropy.table import Table, vstack
import LSS.common_tools as cm
from astropy.io import fits
import numpy as np
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--realization", default=None)
parser.add_argument("--mock_version_forLRG", default='v4.00')
parser.add_argument("--mock_version_forELG", default='v5.0')
parser.add_argument("--mock_version_forQSO", default='v4.00')
parser.add_argument("--mock", default='holi')
parser.add_argument("--LRG_file_name", default='forFA0.fits')
parser.add_argument("--ELG_file_name", default='forFA0_withcontaminants.fits')
parser.add_argument("--QSO_file_name", default='forFA0_withcontaminants.fits')
parser.add_argument("--output_path", default='/pscratch/sd/d/desica/DA2/mocks/holi_v1/')

args = parser.parse_args()

if args.mock == 'holi':
    realizationO = str(args.realization).zfill(4)
    
    elg_path = os.path.join("/global/cfs/cdirs/desi/mocks/cai", args.mock, args.mock_version_forELG, "seed%s" % realizationO, 'ELG', args.ELG_file_name)
    elgs = Table.read(elg_path)

    lrg_path = os.path.join("/global/cfs/cdirs/desi/mocks/cai", args.mock, args.mock_version_forLRG, "seed%s" % realizationO, 'LRG', args.LRG_file_name)
    lrgs = Table.read(lrg_path)

    qso_path = os.path.join("/global/cfs/cdirs/desi/mocks/cai", args.mock, args.mock_version_forQSO, "seed%s" % realizationO, 'QSO', args.QSO_file_name)
    qsos = Table.read(qso_path)

elif args.mock == 'abacus':

    realizationO = str(args.realization).zfill(3)

    elg_path = os.path.join("/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v2.0/AbacusSummit_base_c000_ph%s/CutSky" % realizationO, 'ELG', args.ELG_file_name)
    elgs = Table.read(elg_path)

    lrg_path = os.path.join("/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v2.0/AbacusSummit_base_c000_ph%s/CutSky" % realizationO, 'LRG', args.LRG_file_name)
    lrgs = Table.read(lrg_path)

    qso_path = os.path.join("/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v2.0/AbacusSummit_base_c000_ph%s/CutSky" % realizationO, 'QSO', args.QSO_file_name)
    qsos = Table.read(qso_path)

elif args.mock == 'uchuuref':
    realizationO = str(args.realization).zfill(4)

    elg_path = os.path.join("/global/cfs/cdirs/desi/mocks/cai/Uchuu-SHAM/Y3-v2.0/0000/prep_altmtl", 'ELG', args.ELG_file_name)
    elgs = Table.read(elg_path)

    lrg_path = os.path.join("/global/cfs/cdirs/desi/mocks/cai/Uchuu-SHAM/Y3-v2.0/0000/prep_altmtl", 'LRG', args.LRG_file_name)
    lrgs = Table.read(lrg_path)

    qso_path = os.path.join("/global/cfs/cdirs/desi/mocks/cai/Uchuu-SHAM/Y3-v2.0/0000/prep_altmtl", 'QSO', args.QSO_file_name)
    qsos = Table.read(qso_path)

    #ifil = os.path.join("/global/cfs/cdirs/desi/mocks/cai/Uchuu-SHAM/Y3-v2.0/0000/prep_altmtl", args.tracer, args.file_name)


output_path = os.path.join(args.output_path, 'forFA%s.fits' % args.realization)

if os.path.isfile(output_path):
    print('file', output_path, ' already exist')
    sys.exit()


qso_path = os.path.join(args.output_path, 'qsos')

os.makedirs(qso_path, exist_ok=True)
qsofile = os.path.join(qso_path, 'qso%s.txt' % args.realization)

np.savetxt(qsofile, np.array([qsos['TARGETID'], qsos['RSDZ']]).T, fmt='%d %.3f')
print(f'saving qsos to {qsofile}')

targets = vstack([elgs, lrgs, qsos])

del elgs
del lrgs
del qsos

cm.write_LSS_scratchcp(targets, output_path, extname = 'TARGETS')
fits.setval(output_path, 'OBSCON', value='DARK', ext = 1)

print('DONE', realizationO)

del targets


