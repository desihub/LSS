import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import fitsio
import desimodel

##from fiberassign.targets import (TargetsAvailable)
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
import desimodel.io as dmio


def dofa(ranfile,tile,stamp,outdir,id_,mockrea):
    
    namecomb = os.path.join('/global/cscratch1/sd/acarnero/alt_mtls_masterScriptTest_256dirs_rea{MOCKREA}/Univ{UNIV}/fa/SV3'.format(MOCKREA=mockrea, UNIV=id_), str(stamp), 'fa-%s.sh'%str(tile))
    namecomb.format(stamp=stamp, ts=tile)
    run_fba = os.path.join(outdir,'fa-{ts}.sh'.format(ts=tile))
    fout = open(run_fba,'w')
    fin = open(namecomb,'r').readlines()
    for l in fin:
        if l.startswith('source'):
            pass
        elif l.startswith('fba_run'):
            row = l.split(' ')

            for j,r in enumerate(row):
                if r=='--targets':
                    row[j+1] = ranfile
                if r=='--dir':
                    row[j+1] = outdir
            fout.write(" ".join(row))
        else:
            fout.write(l)

    fout.close()
    os.system('bash %s'%run_fba)

def getfatiles(targetf,tilef,dirout='',dt = '2020-03-10T00:00:00',faver='2.3.0'):
    '''
    will write out fiberassignment files for each tile with the FASSIGN, FTARGETS, FAVAIL HDUS
    these are what are required to determine the geometry of what fiberassign thinks could have been observed and also match to actual observations (though FASSIGN is not really necessary)
    targetf is file with all targets to be run through
    tilef lists the tiles to "assign"
    dirout is the directory where this all gets written out !make sure this is unique for every different target!
    '''                                
    tgs = Targets()
    mver = int(faver[:1])
    if mver < 5:
        load_target_file(tgs,targetf)
    else:
        from fiberassign.targets import TargetTagalong,create_tagalong
        tagalong = create_tagalong()#TargetTagalong([])
        load_target_file(tgs,tagalong,targetf)
    print('loaded target file '+targetf)
    
    hw = load_hardware(rundate=dt)
    tiles = load_tiles(tiles_file=tilef)
    #tgsavail = TargetsAvailable(hw, tgs, tiles, tree)
    #favail = LocationsAvailable(tgsavail)
    #del tree
    print('aure',tiles)
    if mver < 3:
        from fiberassign.targets import (TargetTree)
        tree = TargetTree(tgs, 0.01)
     
    if faver == '2.3.0':
        tgsavail = TargetsAvailable(hw, tgs, tiles, tree)
        favail = LocationsAvailable(tgsavail)
        asgn = Assignment(tgs, tgsavail, favail)
    if faver == '2.4.0' or faver == '2.5.0' or faver == '2.5.1':
        tgsavail = TargetsAvailable(hw, tgs, tiles, tree)
        favail = LocationsAvailable(tgsavail)
        asgn = Assignment(tgs, tgsavail, favail,{}) #this is needed for fiberassign 2.4 and higher(?)
    if mver >= 3 and mver < 5:
        from fiberassign.targets import targets_in_tiles
        tile_targetids, tile_x, tile_y = targets_in_tiles(hw, tgs, tiles)
        tgsavail = TargetsAvailable(hw, tiles, tile_targetids, tile_x, tile_y)
        favail = LocationsAvailable(tgsavail)
        asgn = Assignment(tgs, tgsavail, favail,{}) #this is needed for fiberassign 2.4 and higher(?)
    if mver >= 5:
        from fiberassign.targets import targets_in_tiles
        tile_targetids, tile_x, tile_y = targets_in_tiles(hw, tgs, tiles,tagalong)
        tgsavail = TargetsAvailable(hw, tiles, tile_targetids, tile_x, tile_y)
        favail = LocationsAvailable(tgsavail)
        asgn = Assignment(tgs, tgsavail, favail,{}) #this is needed for fiberassign 2.4 and higher(?)

    print('aure',asgn)
    asgn.assign_unused(TARGET_TYPE_SCIENCE)
    if mver < 5:
        write_assignment_fits(tiles, asgn, out_dir=dirout, all_targets=True)
    else:
        write_assignment_fits(tiles,tagalong, asgn, out_dir=dirout, all_targets=True)
    print('wrote assignment files to '+dirout)	


