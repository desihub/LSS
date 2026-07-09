#!/usr/bin/env python

from astropy.table import Table, vstack
import LSS.common_tools as cm
from astropy.io import fits
# import numpy as np
# import os
import sys
import argparse


def concatenate_tracers(qso_path, elg_path, lrg_path, output_path):
    if False:
        # save only QSO , for debug ?
        # qso_path = os.path.join(args.output_path, 'qsos')
        # os.makedirs(qso_path, exist_ok=True)
        # qsofile = os.path.join(qso_path, 'qso%s.txt' % args.realization)
        # np.savetxt(qsofile, np.array([qsos['TARGETID'], qsos['RSDZ']]).T, fmt='%d %.3f')
        # print(f'saving qsos to {qsofile}')
        pass
    # read Table
    elgs = Table.read(elg_path)    
    lrgs = Table.read(lrg_path)
    qsos = Table.read(qso_path)
    # Concatenate the three tables into one
    targets = vstack([elgs, lrgs, qsos])
    # free 
    del elgs, lrgs, qsos
    # 
    cm.write_LSS_scratchcp(targets, output_path, extname = 'TARGETS')
    fits.setval(output_path, 'OBSCON', value='DARK', ext = 1)
    del targets
    print('DONE: ', output_path)


#
# Main
#

parser = argparse.ArgumentParser()

parser.add_argument("--inputs", nargs="+")
parser.add_argument("--outputs", nargs="+")
args = parser.parse_args() 

print("concatenating tracers to FBA")
print(f"inputs: {args.inputs}")
print(f"outputs: {args.outputs}")

if len(args.inputs) != 3:
    print('Error: must provide exactly 3 input files (QSO, ELG, LRG)')
    sys.exit(1)
    
for f_in in args.inputs:
    if f_in.find('QSO')> 0:
        qso_path = f_in
    elif f_in.find('ELG')> 0:
        elg_path = f_in
    elif f_in.find('LRG')> 0:
        lrg_path = f_in
    else:
        print('Error: input file name must contain QSO, ELG, or LRG')
        sys.exit(1)

concatenate_tracers(qso_path, elg_path, lrg_path, args.outputs[0])