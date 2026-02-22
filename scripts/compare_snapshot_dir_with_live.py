#!/usr/bin/python
'''
Useful script for checking the integrity of a directory

example run:
    python compare_snapshot_dir_with_live.py --livedir /global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/daily/LSScats/1/ --snapdate 2021-08-26

Please note snapshots on cfs are only present for a week. Asking for a snapshot older than a week will cause a missing directory error.
'''
import argparse
import os

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--livedir', type=str, required=True, help="Directory you wish to compare with its snapshot")
parser.add_argument('--snapdate', type=str, required=True, help="Date of snapshot to compare against, format YYYY-MM-DD")

uargs = parser.parse_args()
live = set(os.listdir(uargs.livedir))
snapshot = set(os.listdir(uargs.livedir+'/.snapshots/'+uargs.snapdate))

print(f'Files present in the snapshot {uargs.snapdate} but not in the live version: {snapshot-live}')
print(f'Files present in the live version but not the snapshot {uargs.snapdate}: {live-snapshot}')
