import fitsio
import numpy as np
from astropy.table import Table,join
# system
import os
import subprocess
import sys
import tempfile
import shutil
import re

# time
from time import time
from datetime import datetime, timedelta

# desi
import desitarget
from desitarget.io import read_targets_in_tiles 


def get_fba_mock(mockdir,mocknum,survey='DA02',prog='dark'):
    
    mock_fn = mockdir+'/forFA'+str(mocknum)+'.fits'
    if not os.path.exists(mockdir+'/'+survey):
        os.mkdir(mockdir+'/'+survey)
        print('made '+mockdir+'/'+survey)
    if not os.path.exists(mockdir+'/'+survey+'/fba'+str(mocknum)):
        os.mkdir(mockdir+'/'+survey+'/fba')
        print('made '+mockdir+'/'+survey+'/fba'+str(mocknum))

    tile_fn = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/tiles-'+prog.upper()+'.fits'
    tiles = fitsio.read(tile_fn)
    ts = str(tiles['TILEID'][0]).zfill(6)
    #get info from origin fiberassign file
    fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
    skyf = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/'+'-skies-'+prog.upper()+'.fits'
    outdir = mockdir+'/'+survey+'/fba'+str(mocknum)

    tars = read_targets_in_tiles(mock_fn,tiles)
    tarfn = outdir+'/targs.fits'
    tars.write(tarfn,format='fits',overwrite=True)

    fo = open(outdir+'fa-'+ts+'.sh','w')
    fo.write('#!/bin/bash\n\n')
    fo.write('source /global/common/software/desi/desi_environment.sh master\n')
    fo.write("module swap fiberassign/5.0.0\n")

    fo.write("fba_run")
    fo.write(" --targets "+tarfn)
    fo.write(" --sky "+skyf)
    fo.write(" --footprint "+tile_fn)
    rundate= fht['RUNDATE']
    fo.write(" --rundate "+rundate)
    fo.write(" --fieldrot "+str(fht['FIELDROT']))
    fo.write(" --dir "+outdir)
    fo.write(" --sky_per_petal 40 --standards_per_petal 10")
    fo.write(" --sky_per_slitblock 1")
    fo.write(" --ha "+str(fht['FA_HA']))
    fo.write(" --margin-gfa 0.4 --margin-petal 0.4 --margin-pos 0.05")
    fo.close()    
