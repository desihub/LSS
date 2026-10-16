#!/usr/bin/env python3

#
# adapt join_imaging_mask.py pour snakemake
#

import os
import argparse

from astropy.table import Table, setdiff
from astropy.utils.diff import report_diff_values


import LSS.common_tools as cm
from LSS.globals import main


print('join_imaging_mask_stdpars.py')

mainp = main(tp="LRG", specver="loa-v1")


def process_one_file(in_f, out_f):
    d = Table.read(in_f)
    #targets2 = cm.cutphotmask(d.copy(), bits=mainp.imbits)
    targets = cm.cutphotmask_inplace(d, bits=mainp.imbits)
    
    #targets2["DEC"][0] = 0.0
    #print(setdiff(targets, targets2))
    cm.write_LSS(targets, out_f, extname="TARGETS")


#
# Main
#
parser = argparse.ArgumentParser()

parser.add_argument("--inputs", nargs="+")
parser.add_argument("--outputs", nargs="+")
args = parser.parse_args()


print("joining imaging mask to mock catalogs")
print(f"inputs: {args.inputs}")
print(f"outputs: {args.outputs}")

for i, o in zip(args.inputs, args.outputs):
    if os.path.exists(o):
        print(f"Output file {o} already exists, skipping")
    else:   
        process_one_file(i, o)
