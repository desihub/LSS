#!/usr/bin/env python3

#
# adapt join_imaging_mask.py pour snakemake
#

print('join_imaging_mask_snake.py')

from astropy.table import Table
import LSS.common_tools as cm
from LSS.globals import main
import argparse

mainp = main(tp="LRG", specver="loa-v1")


def process_one_file(in_f, out_f):
    d = Table.read(in_f)
    targets = cm.cutphotmask(d, bits=mainp.imbits)
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
    process_one_file(i, o)
