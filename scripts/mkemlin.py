#standard python
import sys
import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import healpy as hp
import fitsio
import glob
import argparse
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt
from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desimodel.footprint import is_point_in_desi
import desimodel.footprint as foot
from desitarget import targetmask

#import logging
#logging.getLogger().setLevel(logging.ERROR)


#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
from LSS.globals import main

parser = argparse.ArgumentParser()
parser.add_argument("--prog", help="dark or bright is supported",default='dark')

args = parser.parse_args()
print(args)

specrel = 'daily'
prog = args.prog
progu = prog.upper()


mainp = main(prog)

mt = mainp.mtld
tiles = mainp.tiles

wd = mt['SURVEY'] == 'main'
wd &= mt['ZDONE'] == 'true'
wd &= (mt['FAPRGRM'] == prog | mt['FAPRGRM'] == prog+'1b')
mtd = mt[wd]
print('found '+str(len(mtd))+' '+prog+' time main survey tiles with zdone true for '+specrel+' version of reduced spectra')


tiles4comb = Table()
tiles4comb['TILEID'] = mtd['TILEID']
tiles4comb['ZDATE'] = mtd['ARCHIVEDATE']
tiles4comb['THRUDATE'] = mtd['LASTNIGHT']

tiles.keep_columns(['TILEID','RA','DEC'])
#print(tiles.dtype.names)

tiles4comb = join(tiles4comb,tiles,keys=['TILEID'])

print('check that length of tiles4comb matches '+str(len(tiles4comb)))


outdir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/emtiles/'
guadtiles = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/datcomb_'+prog+'_spec_zdone.fits',columns=['TILEID'])
guadtiles = np.unique(guadtiles['TILEID'])
gtids = np.isin(tiles4comb['TILEID'],guadtiles)
tiles4em = tiles4comb[~gtids]
ndone = 0

def mkEMtile(ii):
    if ii >= len(tiles4em):
        print('out of range!')
    else:    
        tile,zdate,tdate = tiles4em['TILEID'][ii],tiles4em['ZDATE'][ii],tiles4em['THRUDATE'][ii]
        outf = outdir+'emline-'+str(tile)+'.fits'
        if not os.path.isfile(outf):
            tdate = str(tdate)
            ct.combEMdata_daily_old(tile,zdate,tdate,outf=outf)
            print('wrote '+outf)

if __name__ == '__main__':

    from multiprocessing import Pool
    N = 64
    if os.environ['NERSC_HOST'] == 'perlmutter':
        N = 128
        print('using 128 cpus')
    for n in range(0,len(tiles4em),N):
        p = Pool(N)
        inds = []
        for i in range(n,n+N):
            inds.append(i)
        p.map(mkEMtile,inds)
        print(n,len(tiles4em))

