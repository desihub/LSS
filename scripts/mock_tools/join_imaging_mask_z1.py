from astropy.table import Table
import LSS.common_tools as cm
from LSS.globals import main
import os
import argparse


mainp = main(tp = 'LRG', specver = 'loa-v1')


parser = argparse.ArgumentParser()
parser.add_argument("--tracer", default = None)
parser.add_argument("--realization", default = None)
parser.add_argument("--mock_version", default = None)
parser.add_argument("--mock", default = 'holi')
parser.add_argument("--file_name", default = 'imforFA0_Y3_noimagingmask_applied.fits')

args = parser.parse_args()

if args.mock == 'holi':
    realization = str(args.realization).zfill(4)
    ifil = os.path.join("/global/cfs/cdirs/desi/mocks/cai", args.mock, args.mock_version, "seed%s" % realization, args.tracer, args.file_name)
elif args.mock == 'abacus':

    realization = str(args.realization).zfill(3)
    ifil = os.path.join("/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v2.0/AbacusSummit_base_c000_ph%s/CutSky" % realization, args.tracer, args.file_name)
    #/global/cfs/projectdirs/desi/mocks/cai/abacus_HF/DR2_v2.0/AbacusSummit_base_c000_ph024/CutSky/QSO/imforFA0_Y3_noimagingmask_applied.fits

elif args.mock == 'uchuuref':
    realization = str(args.realization).zfill(4)
    ifil = os.path.join("/global/cfs/cdirs/desi/mocks/cai/Uchuu-SHAM/Y3-v2.0/0000/prep_altmtl", args.tracer, args.file_name)

output_path = os.path.join(os.path.dirname(ifil), 'forFA0.fits')

print('output file will be', output_path)

if os.path.isfile(output_path):
    print("already exists")
else:
    print("creating file")

    d=Table.read(ifil)                                                                                                                                                   
    targets = cm.cutphotmask(d, bits = mainp.imbits)                                                                                                                                
    cm.write_LSS_scratchcp(targets, output_path, extname='TARGETS')                                                                                              
    print(output_path, 'done')


