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
    #produces script to run to get mock fiberassign files
    mock_fn = mockdir+'/forFA'+str(mocknum)+'.fits'
    if not os.path.exists(mockdir+'/'+survey):
        os.mkdir(mockdir+'/'+survey)
        print('made '+mockdir+'/'+survey)
    if not os.path.exists(mockdir+'/'+survey+'/fba'+str(mocknum)):
        os.mkdir(mockdir+'/'+survey+'/fba'+str(mocknum))
        print('made '+mockdir+'/'+survey+'/fba'+str(mocknum))

    tile_fn = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/tiles-'+prog.upper()+'.fits'
    tiles = Table(fitsio.read(tile_fn,columns=['TILEID','RA','DEC']))
    tiles['OBSCONDITIONS'] = 1
    tiles['IN_DESI'] = 1
    tiles['PROGRAM'] = 'MAIN'
    
    ts = str(tiles['TILEID'][0]).zfill(6)
    #get info from origin fiberassign file
    fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
    skyf = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/skies-'+prog.upper()+'.fits'
    outdir = mockdir+'/'+survey+'/fba'+str(mocknum)
    tile_fn =  outdir+'/tiles.fits'
    tiles.write(tile_fn,overwrite=True)
    tars = read_targets_in_tiles(mock_fn,tiles)
    tarfn = outdir+'/targs.fits'
    Table(tars).write(tarfn,format='fits',overwrite=True)

    fo = open(outdir+'/fa-'+ts+'.sh','w')
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

def get_fba_mock_ran(mockdir,rannum,survey='DA02',prog='dark'):
    #produces script to run to get mock fiberassign files
    from fiberassign.targets import (TargetsAvailable)
    from fiberassign.utils import option_list, GlobalTimers
    from fiberassign.hardware import load_hardware
    from fiberassign.tiles import load_tiles, Tiles
    from fiberassign.targets import (TARGET_TYPE_SCIENCE, TARGET_TYPE_SKY,
                                 TARGET_TYPE_SUPPSKY,
                                 TARGET_TYPE_STANDARD, TARGET_TYPE_SAFE,
                                 Targets, TargetsAvailable,
                                 LocationsAvailable, load_target_file)
    from fiberassign.assign import (Assignment, write_assignment_fits,
                                write_assignment_ascii, merge_results,
                                read_assignment_fits_tile)                                 
 
 
    mock_fn = mockdir+'/ran_forFA'+str(rannum)+'.fits'
    if not os.path.exists(mockdir+'/'+survey):
        os.mkdir(mockdir+'/'+survey)
        print('made '+mockdir+'/'+survey)
    if not os.path.exists(mockdir+'/'+survey+'/fba'+str(mocknum)):
        os.mkdir(mockdir+'/'+survey+'/fba'+str(mocknum))
        print('made '+mockdir+'/'+survey+'/fba'+str(mocknum))

    tile_fn = '/global/cfs/cdirs/desi/survey/catalogs/'+survey+'/LSS/tiles-'+prog.upper()+'.fits'
    tiles = Table(fitsio.read(tile_fn,columns=['TILEID','RA','DEC']))
    tiles['OBSCONDITIONS'] = 1
    tiles['IN_DESI'] = 1
    tiles['PROGRAM'] = 'MAIN'
    
    ts = str(tiles['TILEID'][0]).zfill(6)
    #get info from origin fiberassign file
    fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
    dirout = mockdir+'/'+survey+'/random_fba'+str(rannum)
    tile_fn =  outdir+'/tiles.fits'
    tiles.write(tile_fn,overwrite=True)
    tars = read_targets_in_tiles(mock_fn,tiles)
    tarfn = outdir+'/targs.fits'
    Table(tars).write(tarfn,format='fits',overwrite=True)

    
    from fiberassign.targets import TargetTagalong,create_tagalong
    tagalong = create_tagalong()#TargetTagalong([])
    load_target_file(tgs,tagalong,tarfn)
    print('loaded target file '+tarfn)
    
    hw = load_hardware(rundate=rundate)
    tiles = load_tiles(tiles_file=tile_fn)
    from fiberassign.targets import targets_in_tiles
    tile_targetids, tile_x, tile_y = targets_in_tiles(hw, tgs, tiles,tagalong)
    tgsavail = TargetsAvailable(hw, tiles, tile_targetids, tile_x, tile_y)
    favail = LocationsAvailable(tgsavail)
    asgn = Assignment(tgs, tgsavail, favail,{}) #this is needed for fiberassign 2.4 and higher(?)

    asgn.assign_unused(TARGET_TYPE_SCIENCE)
    write_assignment_fits(tiles,tagalong, asgn, out_dir=dirout, all_targets=True)
    print('wrote assignment files to '+dirout)	
