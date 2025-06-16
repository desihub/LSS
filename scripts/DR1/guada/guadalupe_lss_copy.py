import os
import glob
import sys
import numpy as np

import errno
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--survey",default='Y1')
parser.add_argument("--pubrelease",default='dr1')
parser.add_argument("--specrel",default='guadalupe/v1.0')
parser.add_argument("--outroot",default=None)
args = parser.parse_args()


def my_ln(infile, publink, isdir=False):
    if isdir:
        if os.path.isdir(infile):
            os.system("rsync -av --ignore-existing --exclude 'recon' " + infile + ' ' + publink)
            #os.system('cp -R ' + infile + ' ' + publink)
            return 'ok'
        else:
            raise Exception(infile + ' does not exist, please revise your script')
    else:
        if os.path.isfile(infile):
            os.system('rsync -av --ignore-existing ' + infile + ' ' + publink)
            ##os.system('ln -s ' + infile + ' ' + publink)
            return 'ok'
        else:
            raise Exception(infile + ' does not exist, please revise your script')


if args.outroot is None:
    args.outroot = os.getenv('SCRATCH')

#outroot should be the equivalent to the DESI_ROOT destination. In guadalupe will be https://data.desi.lbl.gov/public/edr


indir = '/global/cfs/cdirs/desi/survey/catalogs/edav1/da02'
pubdir = os.path.join(args.outroot, 'vac', args.pubrelease, 'lss', args.specrel)


def test_dir(value):
    if not os.path.exists(value):
        try:
            os.makedirs(value, 0o755)
#            print('made %s'%value)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

print(pubdir)
test_dir(pubdir)

#in guadalupe lss, only LSScats will be copied
in_path = os.path.join(indir, 'LSScats')

my_ln(in_path, pubdir, isdir=True)




