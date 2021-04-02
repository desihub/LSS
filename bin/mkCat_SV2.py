#standard python
import sys
import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import fitsio
import glob
import argparse
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt

sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.SV2.cattools as ct
#except:
#    print('import of LSS.mkCat_singletile.cattools failed')
#    print('are you in LSS/bin?, if not, that is probably why the import failed')   
import LSS.mkCat_singletile.fa4lsscat as fa


mdir = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/sv2/bright/' #location of ledgers
tdir = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.53.0/targets/sv2/resolve/bright/' #location of targets
mtld = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/mtl-done-tiles.ecsv') #log of tiles completed for mtl

#construct a table with the needed tile information
tilel = []
ral = []
decl = []
mtlt = []
fal = []
obsl = []
for tile in mtld['TILEID']:
    fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/0'+str(tile)[:2]+'/fiberassign-0'+str(tile)+'.fits.gz')
    tilel.append(tile)
    ral.append(fht['TILERA'])
    decl.append(fht['TILEDEC'])
    mtlt.append(fht['MTLTIME'])
    fal.append(fht['FA_RUN'])
    obsl.append(fht['OBSCONDITIONS'])
ta = Table()
ta['TILEID'] = tilel
ta['RA'] = ral
ta['DEC'] = decl
ta['MTLTIME'] = mtlt
ta['FA_RUN'] = fal
ta['OBSCONDITIONS'] = obsl

ct.randomtiles_allSV2(ta)