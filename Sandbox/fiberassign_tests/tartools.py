'''
Some functions to make it easy to generate mtl files, etc., for use in fiberassign tests

'''

import os
import sys
from collections import OrderedDict
import shutil
from datetime import datetime
import glob
import re

import numpy as np
from numpy.lib.recfunctions import append_fields

import matplotlib.pyplot as plt

from astropy.table import Table,join,unique,vstack

from scipy.spatial import KDTree

import fitsio

from desimodel.io import findfile as dm_findfile
from desimodel.io import load_tiles as dm_load_tiles

from desitarget.targetmask import desi_mask, obsconditions

from desitarget.mtl import make_mtl

from fiberassign.targets import (
    Targets,
    TargetTree,
    TargetsAvailable,
    LocationsAvailable,
    load_target_table,
    default_target_masks,
    TARGET_TYPE_SCIENCE, 
    TARGET_TYPE_SKY,
    TARGET_TYPE_SUPPSKY,
    TARGET_TYPE_STANDARD
)

from fiberassign.tiles import (
    load_tiles,
)

from fiberassign.assign import (
    Assignment,
)

from fiberassign.vis import (
    plot_assignment_tile,
)

from fiberassign.qa import qa_targets

from fiberassign.scripts.assign import (
    parse_assign,
    run_assign_full
)

from fiberassign.scripts.merge import (
    parse_merge,
    run_merge
)

# Run the fba_run and fba_merge commandline entrypoints

sky_mask = 0
sky_mask |= desi_mask["SKY"].mask
sky_mask |= desi_mask["SUPP_SKY"].mask
sky_mask |= desi_mask["BAD_SKY"].mask

science_mask = 0
science_mask |= desi_mask["LRG"].mask
science_mask |= desi_mask["ELG"].mask
science_mask |= desi_mask["QSO"].mask

std_mask = 0
std_mask |= desi_mask["STD_FAINT"].mask
std_mask |= desi_mask["STD_WD"].mask
std_mask |= desi_mask["STD_BRIGHT"].mask


def run_assignment(footprint, assign_date = "2020-01-01T00:00:00", indir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/',fullfoot=False,fullsky=''):
    footprint_file = indir+footprint
    science_file = indir + 'mtl_science.fits'
    std_file = indir + 'mtl_std.fits'
    if fullfoot:
        sky_file = fullsky
    else:
        sky_file = indir +'mtl_sky.fits'
    
    outdir=indir+'fiberassign'

    opts = [
        "--rundate", assign_date,
        "--overwrite",
        "--write_all_targets",
        "--footprint", footprint_file,
        "--dir", outdir,
        "--targets", science_file, std_file, sky_file
    ]
    print("  Running raw fiber assignment (fba_run)...")
    ag = parse_assign(opts)
    run_assign_full(ag)
    
    opts = [
        "--skip_raw",
        "--dir", outdir,
        "--targets", science_file, std_file, sky_file
    ]
    print("  Merging input target data (fba_merge_results)...")
    ag = parse_merge(opts)
    run_merge(ag)
    
    return

def update_mtl(obs,oldf='mtl_science_old.fits',science_input='mtl_science.fits',indir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/' ):
    """
    This takes the input MTL and sets the NUMOBS_MORE column based on the
    input dictionary of obs remaining for each target.
    """
    science_input = indir+science_input
    tt = Table.read(science_input)
    tt.write(indir+oldf,format='fits', overwrite=True)
    print("  Loading data from {}".format(science_input), flush=True)
    tdata = None
    with fitsio.FITS(science_input) as fd:
        tdata = fd[1].read()
    
    if "NUMOBS_MORE" not in tdata.dtype.names:
        # create this column based on NUMOBS_INIT
        tdata = append_fields(tdata, "NUMOBS_MORE", tdata["NUMOBS_INIT"])
        
    # Sanity check
    
    if len(obs) != len(tdata):
        msg = "The length of the MTL table ({}) does not match the obs dict ({})".format(
            len(tdata), len(obs)
        )
        raise RuntimeError(msg)
        
    # Now assign the new obs remaining data
    
    print("  Updating observation counts", flush=True)
    tdata["NUMOBS_MORE"][:] = [obs[x] for x in tdata["TARGETID"]]
    
    science_output = science_input
    if os.path.isfile(science_output):
        os.remove(science_output)
        
    print("  Writing updated MTL to {}".format(science_output), flush=True)
    with fitsio.FITS(science_input, "rw") as fd:
        fd.write(tdata)

    del tdata
    return

def mkmtl_assignavail(footprint ,type='ELG',science_input='mtl_science.fits', fba_dir='fiberassign/',indir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/'):#, qso_lyman_rows, qso_tracer_rows):
    #get the unique targetids for assigned and available targets from a set of tiles
    footprint = indir+footprint
    science_input = indir+science_input
    fba_dir = indir+fba_dir
    # Load the footprint
    tile_data = None
    with fitsio.FITS(footprint) as fd:
        tile_data = np.array(fd[1].read())
    

    availids = np.array([])
    assignids = np.array([])
    for tl in tile_data["TILEID"]:
        # For each tile in order of assignment...
        
        # Load assignment and available targets and their properties.
        # NOTE: because we used the --write_all_targets option to fba_run, we get the properties
        # of all available targets in the FTARGETS HDU and have access to those here.
        
        fba_file = os.path.join(fba_dir, "fiberassign-{:06d}.fits".format(tl))
        fassign = None
        ftarget = None
        favail = None
        with fitsio.FITS(fba_file, "r") as fd:
            fassign = fd["FIBERASSIGN"].read()
            ftarget = fd["TARGETS"].read()
            favail = fd["POTENTIAL_ASSIGNMENTS"].read()
        
        # The assigned target IDs
        assign_valid_rows = np.where(fassign["TARGETID"] >= 0)[0]

        assign_tgids = np.sort(fassign["TARGETID"][assign_valid_rows])
        assign_target_rows = np.where(
            np.isin(ftarget["TARGETID"], assign_tgids)
        )[0]

        assign_class_rows = assign_target_rows[
            np.where(
                np.bitwise_and(
                    ftarget["DESI_TARGET"][assign_target_rows],
                    desi_mask[type]
                )
            )[0]
        ]


        assignids = np.concatenate((assignids,ftarget["TARGETID"][assign_class_rows]))

        avail_tgids = np.sort(np.unique(favail["TARGETID"]))
        avail_target_rows = np.where(
            np.isin(ftarget["TARGETID"], avail_tgids)
        )[0]
        avail_class_rows = avail_target_rows[
            np.where(
                np.bitwise_and(
                    ftarget["DESI_TARGET"][avail_target_rows],
                    desi_mask[type]
                )
            )[0]
        ]


        availids = np.concatenate((availids,ftarget["TARGETID"][avail_class_rows]))
        
    assignids = np.unique(assignids)
    w = (assignids > 0) & (assignids*0 == 0)
    assignids = assignids[w]
    print(assignids.dtype.names,assignids.dtype)
    tass = Table()
    tass['TARGETID'] = assignids
    print(len(tass))
    #assignids = Table(assignids,names=['TARGETID'])
    availids = np.unique(availids)
    tav = Table()
    tav['TARGETID'] = availids
    tt = Table.read(science_input)
    #print(len(tt),len(np.unique(tt['TARGETID'])))
    #wass = np.where(
    #    np.isin(tt['TARGETID'],assignids)
    #)[0]
    ttass = join(tass,tt,keys=['TARGETID'],join_type='left')
    ttass = unique(ttass,keys=['TARGETID'])
    print('number of assigned '+type)
    print(len(ttass['TARGETID']),len(assignids),len(np.unique(ttass['TARGETID'])))
    
    ttav = join(tav,tt,keys=['TARGETID'],join_type='left')
    ttav = unique(ttav,keys=['TARGETID'])
    #wave = np.where(np.isin(tt['TARGETID'],availids))[0]
    print('number of available '+type)
    print(len(ttav['TARGETID']),len(availids),len(np.unique(ttav['TARGETID'])))
    plt.plot(ttav['RA'],ttav['DEC'],'k,')
    plt.plot(ttass['RA'],ttass['DEC'],'r,')
    plt.xlim(18,24)
    plt.ylim(7,13)
    plt.show()

def getall_fassign(type,indir,nmonths=70,cadence=28):
    fba_files0 = glob.glob(os.path.join(indir+'0000/',"fba-*.fits"))
    fah = fitsio.read_header(fba_files0[0])
    tile = fah['TILEID']
    fass = fitsio.read(fba_files0[0],ext='FASSIGN')
    
    if type == 'SKY':
        wsk = ((fass['FA_TARGET'] & 2**37) > 0) | ((fass['FA_TARGET'] & 2**36) > 0) | ((fass['FA_TARGET'] & 2**32) > 0)
    fass = fass[wsk]
    tl = np.ones(len(fass),dtype=int)*tile
    #fass = append_fields(fass,'TILE',tl)

    #fass['TILE'] = tile
    for i in range(1,len(fba_files0)):
        fah = fitsio.read_header(fba_files0[i])
        tile = fah['TILEID']
        fai = fitsio.read(fba_files0[i],ext='FASSIGN')
        if type == 'SKY':
            wsk = ((fai['FA_TARGET'] & 2**37) > 0) | ((fai['FA_TARGET'] & 2**36) > 0) | ((fai['FA_TARGET'] & 2**32) > 0)
        fai = fai[wsk]
        #print(len(fai))
        tli = np.ones(len(fai),dtype=int)*tile
        #fai = append_fields(fai,'TILE',tl)
        #fai['TILE'] = tile
        #fass = vstack([fass,fai],metadata_conflicts='silent')
        #print(fass.dtype)
        #print(fai.dtype)
        fass = np.hstack((fass,fai))
        tl = np.hstack((tl,tli))
    fb = np.zeros(len(fass))
    fm = np.ones((len(fass)),dtype=int)*cadence
    #fass['BATCH'] = 0
    #fass['MAXSURVEYMJD'] = cadence
    
    #fass = append_fields(fass, 'BATCH', fb, usemask=False) 
    #fass = append_fields(fass, 'MAXSURVEYMJD', fm, usemask=False) 
    for j in range(1,nmonths):
        print('working on batch '+str(j))
        m = str.zfill(str(j),4)
        fba_filesj = glob.glob(indir+m+"/fba-*.fits")
        for i in range(0,len(fba_filesj)):
            fah = fitsio.read_header(fba_filesj[i])
            tile = fah['TILEID']
            fai = fitsio.read(fba_filesj[i],ext='FASSIGN')
            if type == 'SKY':
                wsk = ((fai['FA_TARGET'] & 2**37) > 0) | ((fai['FA_TARGET'] & 2**36) > 0) | ((fai['FA_TARGET'] & 2**32) > 0)
            fai = fai[wsk]
            #print(len(fai))
            #fai['TILE'] = tile
            tli = np.ones(len(fai),dtype=int)*tile
            tl = np.hstack((tl,tli))
            #fai = append_fields(fai,'TILE',tl)
            fbi = np.ones(len(fai))*j
            fb = np.hstack((fb,fbi))
            #fai = append_fields(fai,'BATCH',fb)
            fmi = np.ones((len(fai)),dtype=int)*cadence*(j+1)
            fm = np.hstack((fm,fmi))
            #fai = append_fields(fai,'MAXSURVEYMJD',fm)

            #fass = vstack([fass,fai],metadata_conflicts='silent')
            fass = np.hstack((fass,fai))

        #fass['BATCH'] = j
        #fass['MAXSURVEYMJD'] = cadence*(j+1) 

        print('after batch '+str(j)+ ' there are '+str(len(fass))+' '+type+' assignments')  
    fass = append_fields(fass,'TILE',tl, usemask=False)
    fass = append_fields(fass,'BATCH',fb, usemask=False)
    fass = append_fields(fass,'MAXSURVEYMJD',fm, usemask=False)
    print(np.unique(fass['BATCH']))
    #fass.write(indir+'all_assigned_'+type+'.fits',format='fits', overwrite=True)
    outf = indir+'all_assigned_'+type+'.fits'
    if os.path.isfile(outf):
        os.remove(outf)

    with fitsio.FITS(outf, "rw") as fd:
        fd.write(fass)

    return True
        
           
        
    

#just get all of the excess sky counts for the fiberassign files in a directory
def sky_counts(indir,nskym=400,nscix = 4500):
    fba_files = glob.glob(os.path.join(indir,"fba-*.fits"))
    next = 0
    ni = 0
    for fl in fba_files:
        fass = fitsio.read(fl,ext='FASSIGN')
        wv = (fass["TARGETID"] >= 0 ) & ((fass['DEVICE_TYPE'] == b'POS') | (fass['DEVICE_TYPE'] == 'POS'))
        fass = fass[wv]
        if len(fass) > 4800:
            wsk = ((fass['FA_TARGET'] & 2**37) > 0) | ((fass['FA_TARGET'] & 2**36) > 0) | ((fass['FA_TARGET'] & 2**32) > 0) | ((fass['FA_TARGET'] & 2**61) > 0)
            wsk &=  ((fass['FA_TARGET'] & 2**2) == 0) & ((fass['FA_TARGET'] & 2**1) == 0) & ((fass['FA_TARGET'] & 2**0) == 0)
            ws = ((fass['FA_TARGET'] & 2**2) > 0) | ((fass['FA_TARGET'] & 2**1) > 0) | ((fass['FA_TARGET'] & 2**0) > 0) | ((fass['FA_TARGET'] & 2**60) > 0) | ((fass['FA_TARGET'] & 2**61) > 0)
            nskyi = len(fass[wsk])
            nscii = len(fass[ws])
            nexti = nskyi-nskym
            nextis = nscix-nscii
            #print(fl,nexti,nextis,len(fass))
            next += nexti
        else:
            ni += 1
            print('not enough assignments',fl,len(fass)) 
    print('total number of extra fibers '+str(next)+ ' across '+str(len(fba_files))+' tiles')
    return next,len(fba_files)-ni

def science_counts(indir):
    fba_files = glob.glob(os.path.join(indir,"fba-*.fits"))
    n0 = 0
    n1 = 0
    n2 = 0
   
    for fl in fba_files:
        fass = fitsio.read(fl,ext='FASSIGN')
        wv = (fass["TARGETID"] >= 0 ) & ((fass['DEVICE_TYPE'] == b'POS') | (fass['DEVICE_TYPE'] == 'POS'))
        fass = fass[wv]
        w0 = ((fass['FA_TARGET'] & 2**0) > 0) 
        w1 = ((fass['FA_TARGET'] & 2**1) > 0) 
        w2 = ((fass['FA_TARGET'] & 2**2) > 0)
        n0 += len(fass[w0])
        n1 += len(fass[w1])
        n2 += len(fass[w2])
    print('total number of science assignments '+str(n0+n1+n2)+ ' across '+str(len(fba_files))+' tiles')
    return n0,n1,n2

def getall_science_counts(indir,nmonths=13,splot=True,title='no pass with 28 day cadence'):
    nt0 = 0
    nt1 = 0
    nt2 = 0
    tl = []
    nl0 = []
    nl1 = []
    nl2 = []
    for i in range(0,nmonths): 
        m = str.zfill(str(i),4)
        n0,n1,n2 = science_counts(indir+m)
        nt0 += n0
        nt1 += n1
        nt2 += n2
        nl0.append(nt0)
        nl1.append(nt1)
        nl2.append(nt2)
        tl.append((i+1)/13.)   
    tl = np.array(tl)
    nl1 = np.array(nl1)
    nl2 = np.array(nl2)
    nl0 = np.array(nl0)
    if splot:
        plt.plot(tl,nl1,'b-',label='ELGS')
        plt.plot(tl,nl0,'r-',label='LRGs')
        plt.plot(tl,nl2,'-',color='purple',label='QSOs')
        plt.plot(tl,tl*10.e6,'k-',label='10e6/year')
        plt.plot(tl,tl*3.e6,'k--',label='3e6/year')
        plt.plot(tl,tl*5.e6,'k:',label='5e6/year')
        #plt.ticklabel_format(style='sci', scilimits=(3,3))
        plt.xlabel('time (years)')
        plt.ylabel('cumulative number of targets')
        #plt.yscale('log')
        plt.title(title)
        plt.legend()
        plt.show()
        plt.plot(tl,nl1/nt1,'b-',label='ELGS')
        plt.plot(tl,nl0/nt0,'r-',label='LRGs')
        plt.plot(tl,nl2/nt2,'-',color='purple',label='QSOs')
        plt.plot(tl,tl*.24,'k--',label='0.24/year')
        plt.ticklabel_format(style='sci', scilimits=(3,3))
        plt.xlabel('time (years)')
        plt.ylabel('fraction completed')
        #plt.yscale('log')
        plt.title(title)
        plt.legend()
        plt.show()

    return True


def getall_sky_counts(indir,nmonths=13,splot=True,title='no pass with 28 day cadence'):
    ne = 0
    nt = 0
    nsl = []
    tl = []
    nel = []
    for i in range(0,nmonths): 
        m = str.zfill(str(i),4)
        nf,nti = sky_counts(indir+m)
        ne += nf
        nel.append(ne)
        nt += nti
        nsl.append(nf/(5000*nti))    
        tl.append((i+1)/13.)   
    print(ne,nt)
    if splot:
        plt.plot(tl,nel,'k-')
        #plt.ticklabel_format(style='sci', scilimits=(3,3))
        plt.xlabel('time (years)')
        plt.ylabel('cumulative number of spare fibers')
        plt.title(title)
        plt.show()
        plt.plot(tl,nsl,'k-')
        plt.xlabel('time (years)')
        plt.ylabel('fraction of spare fibers in previous 28 days')
        plt.title(title)
        plt.show()
    return nsl

# Function to compute the assigned, available, and considered targets for a set of tiles

def assignment_counts(footprint, science_input='mtl_science.fits', fba_dir='fiberassign/',indir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/'):#, qso_lyman_rows, qso_tracer_rows):
    
    footprint = indir+footprint
    science_input = indir+science_input
    fba_dir = indir+fba_dir
    # Load the footprint
    tile_data = None
    with fitsio.FITS(footprint) as fd:
        tile_data = np.array(fd[1].read())
    
    # Load the input science MTL and get the obs remaining for all targets
    mtldata = None
    with fitsio.FITS(science_input) as fd:
        mtldata = fd[1].read()

    obs = None
    if "NUMOBS_MORE" in mtldata.dtype.names:
        obs = {
            x: y for x, y in zip(mtldata["TARGETID"], mtldata["NUMOBS_MORE"])
        }
    else:
        obs = {
            x: y for x, y in zip(mtldata["TARGETID"], mtldata["NUMOBS_INIT"])
        }
    
    qso_rows = (mtldata['DESI_TARGET'] & desi_mask["QSO"].mask) > 0
    qso_lyman_rows = qso_rows & (mtldata['IS_LYA'] == 1)
    qso_tracer_rows = qso_rows &  (mtldata['IS_LYA'] != 1)
    print('all qso, lyman qso, tracer qso')
    print(np.sum(qso_rows),np.sum(qso_lyman_rows),np.sum(qso_tracer_rows))
    qso_lyman = mtldata["TARGETID"][qso_lyman_rows]
    qso_tracer = mtldata["TARGETID"][qso_tracer_rows]

    class_masks = {
        "ELG": desi_mask["ELG"].mask,
        "LRG": desi_mask["LRG"].mask,
        "QSO-lyman": desi_mask["QSO"].mask,
        "QSO-tracer": desi_mask["QSO"].mask,
        "STD": std_mask,
        "SKY": sky_mask
    }
    
    # histogram data to return
    hist_tgassign = dict()
    hist_tgavail = dict()
    hist_tgconsid = dict()
    hist_tgfrac = dict()
    for tgclass, mask in class_masks.items():
        hist_tgassign[tgclass] = list()
        hist_tgavail[tgclass] = list()
        hist_tgconsid[tgclass] = list()
        hist_tgfrac[tgclass] = list()
    
    print("  Accumulating assignment counts for {} tiles...".format(len(tile_data)), flush=True)

    nuelg = np.array([])
    naelg = np.array([])
    for tl in tile_data["TILEID"]:
        # For each tile in order of assignment...
        
        # Load assignment and available targets and their properties.
        # NOTE: because we used the --write_all_targets option to fba_run, we get the properties
        # of all available targets in the FTARGETS HDU and have access to those here.
        
        fba_file = os.path.join(fba_dir, "fiberassign-{:06d}.fits".format(tl))
        fassign = None
        ftarget = None
        favail = None
        with fitsio.FITS(fba_file, "r") as fd:
            fassign = fd["FIBERASSIGN"].read()
            ftarget = fd["TARGETS"].read()
            favail = fd["POTENTIAL_ASSIGNMENTS"].read()
        
        # The assigned target IDs
        assign_valid_rows = np.where(fassign["TARGETID"] >= 0)[0]
        assign_tgids = np.sort(fassign["TARGETID"][assign_valid_rows])
        assign_target_rows = np.where(
            np.isin(ftarget["TARGETID"], assign_tgids)
        )[0]
        
        # The available target IDs
        avail_tgids = np.sort(np.unique(favail["TARGETID"]))
        avail_target_rows = np.where(
            np.isin(ftarget["TARGETID"], avail_tgids)
        )[0]
        
        # For the science classes, we must also look at the obs remaining
        # in order to know which targets were actually considered for assignment
        # (not just reachable).
        
        for tgclass, mask in class_masks.items():
            # The assigned targets in this class
            assign_class_rows = assign_target_rows[
                np.where(
                    np.bitwise_and(
                        ftarget["DESI_TARGET"][assign_target_rows],
                        mask
                    )
                )[0]
            ]

            if tgclass == "QSO-lyman":
                assign_class_rows = assign_class_rows[
                    np.where(
                        np.isin(
                            ftarget["TARGETID"][assign_class_rows], qso_lyman
                        )
                    )[0]
                ]
            elif tgclass == "QSO-tracer":
                assign_class_rows = assign_class_rows[
                    np.where(
                        np.isin(
                            ftarget["TARGETID"][assign_class_rows], qso_tracer
                        )
                    )[0]
                ]
            
            hist_tgassign[tgclass].append(len(assign_class_rows))
                
            # The available targets in this class
            avail_class_rows = avail_target_rows[
                np.where(
                    np.bitwise_and(
                        ftarget["DESI_TARGET"][avail_target_rows],
                        mask
                    )
                )[0]
            ]
            if tgclass == "QSO-lyman":
                avail_class_rows = avail_class_rows[
                    np.where(
                        np.isin(
                            ftarget["TARGETID"][avail_class_rows], qso_lyman
                        )
                    )[0]
                ]
            elif tgclass == "QSO-tracer":
                avail_class_rows = avail_class_rows[
                    np.where(
                        np.isin(
                            ftarget["TARGETID"][avail_class_rows], qso_tracer
                        )
                    )[0]
                ]

            hist_tgavail[tgclass].append(len(avail_class_rows))
            if tgclass == 'ELG':
                nuelg = np.concatenate((nuelg,ftarget["TARGETID"][avail_class_rows]))
                naelg = np.concatenate((naelg,ftarget["TARGETID"][assign_class_rows]))
            
            #print("  target class {}, {} assignments".format(tgclass, len(assign_class_rows)))
            
            if tgclass == "STD" or tgclass == "SKY":
                # the considered targets are the same as the reachable
                hist_tgconsid[tgclass].append(len(avail_class_rows))
                hist_tgfrac[tgclass].append(len(assign_class_rows) / len(avail_class_rows))
            else:
                # compare to obs remaining
                hist_tgconsid[tgclass].append(
                    np.sum(
                        [1 for x in ftarget["TARGETID"][avail_class_rows] if obs[x] > 0]
                    )
                )
                hist_tgfrac[tgclass].append(len(assign_class_rows) / hist_tgconsid[tgclass][-1])

                # Now reduce the obs remaining
                for tgid in ftarget["TARGETID"][assign_class_rows]:
                    if tgclass != 'QSO-lyman':
                        obs[tgid] = 0
                    else:    
                        obs[tgid] -= 1
    
    print('number of unique available ELG targets:')
    print(len(np.unique(nuelg)))
    print('number of unique assigned ELG targets:')
    print(len(np.unique(naelg)))

    # Return our histogram of tile data and also the updated observation counts,
    # which can be used to update the MTL NUMOBS_MORE in a separate function.
    return (obs, hist_tgassign, hist_tgavail, hist_tgconsid, hist_tgfrac)

def mktilefile(obscon=[1,2],target_ra_min=0,target_ra_max=360,target_dec_min=-90,target_dec_max=90,outdir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/'):
    tfn  = os.getenv('DESIMODEL')+'/data/footprint/desi-tiles.fits'
    footprint_data = fitsio.read(tfn)
    tilefile_pass = dict()
    #footprint_file = dm_findfile("footprint/desi-tiles.fits")
    #footprint_data = dm_load_tiles(tilesfile=footprint_file, cache=False)

    tile_radius = 1.65 # degrees
    tile_cut = 2.0 # degrees

    tile_ra_min = target_ra_min #+ tile_cut #this will work near equator...not really bothered if we end up have tiles with no targets on them
    tile_ra_max = target_ra_max #- tile_cut
    tile_dec_min = target_dec_min #+ tile_cut
    tile_dec_max = target_dec_max #- tile_cut

    #obskeep = obsconditions[program]

#     inside = np.where(
#         #np.logical_and(
#             np.logical_and(
#                 np.logical_and(
#                     (footprint_data["RA"] > tile_ra_min), 
#                     (footprint_data["RA"] < tile_ra_max)
#                 ), np.logical_and(
#                     (footprint_data["DEC"] > tile_dec_min), 
#                     (footprint_data["DEC"] < tile_dec_max)
#                 )
#             )
#         #)
#     )[0]
    
    inside = (footprint_data["RA"] > tile_ra_min) & (footprint_data["RA"] < tile_ra_max)
    inside &= (footprint_data["DEC"] > tile_dec_min) & (footprint_data["DEC"] < tile_dec_max)
    inside &= np.isin(footprint_data['OBSCONDITIONS'],obscon)
    inside &= (footprint_data['IN_DESI']==1)

    tiledata = footprint_data[inside]
    
    #tt = Table(tiledata)

    # For each pass, write out a tile file.  Also write out the file for all passes.

    passes = np.unique(tiledata["PASS"])    

    print("Full footprint has {} tiles".format(len(tiledata)))

    tilefile_pass["ALL"] = outdir+'tile_ALL.fits'
    if os.path.isfile(tilefile_pass["ALL"]):
        os.remove(tilefile_pass["ALL"])

    outfd = fitsio.FITS(tilefile_pass["ALL"], "rw")
    outfd.write(None, header=None, extname="PRIMARY")
    outfd.write(tiledata, header=None, extname="TILES")
    outfd.close()


    for ps in passes:
        pstr = "{}".format(ps)
        ps_rows = np.where(tiledata["PASS"] == ps)[0]
        tiledata_pass = tiledata[ps_rows]
        print("Pass {} footprint has {} tiles".format(ps, len(tiledata_pass)))
        tilefile_pass[pstr] = outdir+'tile_'+str(ps)+'.fits'
        if os.path.isfile(tilefile_pass[pstr]):
            os.remove(tilefile_pass[pstr])
        outfd = fitsio.FITS(tilefile_pass[pstr], "rw")
        outfd.write(None, header=None, extname="PRIMARY")
        outfd.write(tiledata_pass, header=None, extname="TILES")
        outfd.close()

    pstr = "2-4"
    ps_rows = np.where(tiledata["PASS"] > 1)[0]
    tiledata_pass = tiledata[ps_rows]
    print("Pass 2-4 footprint has {} tiles".format(len(tiledata_pass)))
    tilefile_pass[pstr] = outdir+'tile_'+pstr+'.fits'
    if os.path.isfile(tilefile_pass[pstr]):
        os.remove(tilefile_pass[pstr])
    outfd = fitsio.FITS(tilefile_pass[pstr], "rw")
    outfd.write(None, header=None, extname="PRIMARY")
    outfd.write(tiledata_pass, header=None, extname="TILES")
    outfd.close()
    
    
def mktarfile(target_ra_min=0,target_ra_max=360,target_dec_min=-90,target_dec_max=90,dr ='dr8',tarver = '0.39.0',outdir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/',prog='dark'):

    # First select the science targets

    input_dir = "/global/cfs/projectdirs/desi/target/catalogs/"+dr+"/"+tarver+"/targets/main/resolve/"+prog

    print('working with target data in '+input_dir+' IS THAT CORRECT???')

    sky_dir = "/global/cfs/projectdirs/desi/target/catalogs/"+dr+"/"+tarver+"/skies"

    print('working with skies data in '+sky_dir+' IS THAT CORRECT???')


    input_files = glob.glob(os.path.join(input_dir, "*.fits"))

    target_data = []

    for file in input_files:
        print("Working on {}".format(os.path.basename(file)), flush=True)
        fd = fitsio.FITS(file, "r")
        fdata = fd[1].read()
        inside = np.where(
            np.logical_and(
                np.logical_and((fdata["RA"] > target_ra_min), (fdata["RA"] < target_ra_max)),
                np.logical_and((fdata["DEC"] > target_dec_min), (fdata["DEC"] < target_dec_max))
            )
        )[0]
        target_data.append(fdata[inside])
        fd.close()

    target_data = np.concatenate(target_data)

    out_file = outdir+"target_science_sample.fits" 
    if os.path.isfile(out_file):
        os.remove(out_file)

    fd = fitsio.FITS(out_file, "rw")
    fd.write(None, header=None, extname="PRIMARY")
    fd.write(target_data, header=None, extname="TARGETS")
    fd.close()
    
    del target_data

    # Now select the sky targets

    print("Working on sky...", flush=True)

    input_files = glob.glob(os.path.join(sky_dir, "*.fits"))

    sky_data = []

    for file in input_files:
        print("Working on {}".format(os.path.basename(file)), flush=True)
        fd = fitsio.FITS(file, "r")
        fdata = fd[1].read()
        inside = np.where(
            np.logical_and(
                np.logical_and((fdata["RA"] > target_ra_min), (fdata["RA"] < target_ra_max)),
                np.logical_and((fdata["DEC"] > target_dec_min), (fdata["DEC"] < target_dec_max))
            )
        )[0]
        sky_data.append(fdata[inside])
        fd.close()

    sky_data = np.concatenate(sky_data)


    out_file = outdir+"target_sky_sample.fits" 
    if os.path.isfile(out_file):
        os.remove(out_file)


    outfd = fitsio.FITS(out_file, "rw")
    outfd.write(None, header=None, extname="PRIMARY")
    outfd.write(sky_data, header=None, extname="TARGETS")
    outfd.close()

    fd.close()

def add_lya(frac=0.2,indir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/'):
    science_file = indir + 'mtl_science.fits'
    ff = fitsio.read(science_file)
    isly = np.zeros(len(ff))
    isly1 = np.ones(len(ff))
    wq = (ff['DESI_TARGET'] & desi_mask["QSO"].mask) > 0
    rv = np.random.rand(len(ff))
    wly = wq & (rv < frac)
    isly[wly] = isly1[wly]
    print('number of targets selected to be lyman alpha:')
    print(np.sum(isly))
    del ff
    tf = Table.read(science_file)
    tf['IS_LYA'] = isly
    tf.write(science_file,format='fits', overwrite=True)
    return True

def splitdarkgray(grayfrac=0.3,indir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/'):
    #from random import random
    science_file = indir + 'mtl_science.fits'
    ft = fitsio.read(science_file)
    rv = np.random.rand(len(ft))
    we = ft['OBSCONDITIONS'] == 3
    #w1 = we & (rv >= grayfrac)
    #w2 = we & (rv < grayfrac)
    #ft[w1]['OBSCONDITIONS'] = 1
    #ft[w2]['OBSCONDITIONS'] = 2
    print('number of targets that were allowed to be observed in either dark or gray')
    print(len(ft[we]))
    for i in range(0,len(ft)):
        if ft[i]['OBSCONDITIONS'] == 3:
            if rv[i] < grayfrac:
                ft[i]['OBSCONDITIONS'] = 2
            else:
                ft[i]['OBSCONDITIONS'] = 1
    we = ft['OBSCONDITIONS'] == 2
    print('number of targets that are now allowed to be observed only in gray')
    print(len(ft[we]))
    # Write MTLs

    if os.path.isfile(science_file):
        os.remove(science_file)
    with fitsio.FITS(science_file, "rw") as fd:
        fd.write(ft)
    

def mkmtl(obscon="DARK|GRAY",target_ra_min=0,target_ra_max=360,target_dec_min=-90,target_dec_max=90,outdir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/',target_sample='/project/projectdirs/desi/users/ajross/dr8tar/target_science_sample.fits'):
    '''
    initially copied from https://github.com/desihub/tutorials/blob/master/FiberAssignAlgorithms_Part2.ipynb
    
    '''
    
    science_file = outdir + 'mtl_science.fits'
    std_file = outdir + 'mtl_std.fits'

    
    # Load the raw science / standard target sample and prune columns

    keep_columns = [
        'TARGETID', 
        'RA', 
        'DEC',
        'RA_IVAR',
        'DEC_IVAR',
        'PMRA',
        'PMDEC',
        'PMRA_IVAR',
        'PMDEC_IVAR',
        'DESI_TARGET', 
        'BGS_TARGET', 
        'MWS_TARGET', 
        'SUBPRIORITY', 
        'BRICKNAME',
        'BRICKID',
        'BRICK_OBJID',
        'PRIORITY_INIT', 
        'NUMOBS_INIT'
    ]

    fd = fitsio.FITS(target_sample)
    fdata = fd[1].read(columns=keep_columns)

    inside = np.where(
        np.logical_and(
            np.logical_and((fdata["RA"] > target_ra_min), (fdata["RA"] < target_ra_max)),
            np.logical_and((fdata["DEC"] > target_dec_min), (fdata["DEC"] < target_dec_max))
        )
    )[0]
    targets_raw = fdata[inside]


    # Get the default target masks for this target file

    (filesurvey, 
     filecol, 
     def_sciencemask, 
     def_stdmask, 
     def_skymask, 
     def_suppskymask,
     def_safemask, 
     def_excludemask) = default_target_masks(targets_raw)

    print("Detected targets for survey '{}', using bitfield column '{}'".format(filesurvey, filecol))

    # Force our science and std masks to a more restrictive set.  Only keep ELG, LRG and QSO targets.
    # Cut any targets with multiple of those set.

#     science_mask = 0
#     science_mask |= desi_mask["LRG"].mask
#     science_mask |= desi_mask["ELG"].mask
#     science_mask |= desi_mask["QSO"].mask
# 
#     std_mask = 0
#     std_mask |= desi_mask["STD_FAINT"].mask
#     std_mask |= desi_mask["STD_WD"].mask
#     std_mask |= desi_mask["STD_BRIGHT"].mask
    
    elg_rows = np.where(
        np.logical_and(
            np.logical_and(
                np.logical_and(
                    np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["ELG"].mask),
                    np.logical_not(
                        np.bitwise_and(targets_raw["DESI_TARGET"], std_mask)
                    )
                ),
                np.logical_not(
                    np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["QSO"].mask)
                )
            ),
            np.logical_not(
                np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["LRG"].mask)
            )
        )
    )[0]

    qso_rows = np.where(
        np.logical_and(
            np.logical_and(
                np.logical_and(
                    np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["QSO"].mask),
                    np.logical_not(
                        np.bitwise_and(targets_raw["DESI_TARGET"], std_mask)
                    )
                ),
                np.logical_not(
                    np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["ELG"].mask)
                )
            ),
            np.logical_not(
                np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["LRG"].mask)
            )
        )
    )[0]

    lrg_rows = np.where(
        np.logical_and(
            np.logical_and(
                np.logical_and(
                    np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["LRG"].mask),
                    np.logical_not(
                        np.bitwise_and(targets_raw["DESI_TARGET"], std_mask)
                    )
                ),
                np.logical_not(
                    np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["QSO"].mask)
                )
            ),
            np.logical_not(
                np.bitwise_and(targets_raw["DESI_TARGET"], desi_mask["ELG"].mask)
            )
        )
    )[0]
    
    n_elg = len(elg_rows)
    n_qso = len(qso_rows)
    n_lrg = len(lrg_rows)

    science_rows = np.concatenate([elg_rows, qso_rows, lrg_rows])

    std_rows = np.where(
        np.logical_and(
            np.bitwise_and(targets_raw["DESI_TARGET"], std_mask),
            np.logical_not(
                np.bitwise_and(targets_raw["DESI_TARGET"], science_mask)
            )
        )
    )[0]

    print(
        "Using {} science and {} standards from input catalog".format(
            len(science_rows),
            len(std_rows)
        )
    )

    # Split out the science and standard targets, although this is actually not necessary for passing
    # to fiberassign.

    science_targets = np.array(targets_raw[science_rows])

    std_targets = np.array(targets_raw[std_rows])

    # Close the input fits file so it doesn't take up extra memory
    del targets_raw
    fd.close()
    del fd

    # We have concatenated the 3 target types in the new table, so now the rows are
    # different:
    elg_rows = np.arange(n_elg, dtype=np.int64)
    qso_rows = np.arange(n_qso, dtype=np.int64) + n_elg
    lrg_rows = np.arange(n_lrg, dtype=np.int64) + n_elg + n_qso

    # Make the MTLs

    
    science_mtl = make_mtl(science_targets, "DARK|GRAY").as_array()
    if len(science_mtl) != len(science_targets):
        print("WARNING:  science MTL has {} rows, input has {}".format(len(science_mtl), len(science_targets)))
    
        

    std_mtl = make_mtl(std_targets, "DARK|GRAY").as_array()
    if len(std_mtl) != len(std_targets):
        print("WARNING:  standards MTL has {} rows, input has {}".format(len(std_mtl), len(std_targets)))

    # Delete the large intermediate arrays
    
    del science_targets
    del std_targets
    
    # Write MTLs

    if os.path.isfile(science_file):
        os.remove(science_file)
    with fitsio.FITS(science_file, "rw") as fd:
        fd.write(science_mtl)

    if os.path.isfile(std_file):
        os.remove(std_file)
    with fitsio.FITS(std_file, "rw") as fd:
        fd.write(std_mtl)    

    print("{} science targets".format(len(science_mtl)))
    print("    {} ELG targets".format(len(elg_rows)))
    print("    {} QSO targets".format(len(qso_rows)))
    print("    {} LRG targets".format(len(lrg_rows)))
    print("{} std targets".format(len(std_mtl)))

    # We'll be loading later science MTLs as we go through the survey, so delete that now.
    # the standards are constant so we'll keep those in memory.

    del science_mtl
    
def mkmtl_sky(target_ra_min=0,target_ra_max=360,target_dec_min=-90,target_dec_max=90,outdir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/',target_sample='/project/projectdirs/desi/users/ajross/dr8tar/target_sky_sample.fits'):
    # Now create the skies file.

    sky_file = outdir + 'mtl_sky.fits'
    
#     keep_columns = [
#         'TARGETID', 
#         'RA', 
#         'DEC', 
#         'DESI_TARGET', 
#         'BGS_TARGET', 
#         'MWS_TARGET', 
#         'SUBPRIORITY', 
#         'BRICKNAME',
#         'BRICKID',
#         'BRICK_OBJID',
#         'FLUX_G',
#         'FLUX_R',
#         'FLUX_Z',
#         'FLUX_IVAR_G',
#         'FLUX_IVAR_R',
#         'FLUX_IVAR_Z',
#         'OBSCONDITIONS'
#     ]


    fd = fitsio.FITS(target_sample)
    fdata = np.array(fd[1].read())#(columns=keep_columns))

    inside = np.where(
        np.logical_and(
            np.logical_and((fdata["RA"] > target_ra_min), (fdata["RA"] < target_ra_max)),
            np.logical_and((fdata["DEC"] > target_dec_min), (fdata["DEC"] < target_dec_max))
        )
    )[0]
    sky_mtl = fdata[inside]


    fd.close()
    del fd

    # Sanity check that these are all sky, supp_sky, or bad_sky

    print("{} input targets in sky file".format(len(sky_mtl)))

    sky_sky_rows = np.where(
        np.bitwise_and(sky_mtl["DESI_TARGET"], desi_mask["SKY"].mask)
    )[0]

    print("  {} SKY targets".format(len(sky_sky_rows)))

    sky_suppsky_rows = np.where(
        np.bitwise_and(sky_mtl["DESI_TARGET"], desi_mask["SUPP_SKY"].mask)
    )[0]

    print("  {} SUPP_SKY targets".format(len(sky_suppsky_rows)))

    sky_badsky_rows = np.where(
        np.bitwise_and(sky_mtl["DESI_TARGET"], desi_mask["BAD_SKY"].mask)
    )[0]

    print("  {} BAD_SKY targets".format(len(sky_badsky_rows)))

#     sky_mask = 0
#     sky_mask |= desi_mask["SKY"].mask
#     sky_mask |= desi_mask["SUPP_SKY"].mask
#     sky_mask |= desi_mask["BAD_SKY"].mask

    sky_unknown_rows = np.where(
        np.logical_not(
            np.bitwise_and(sky_mtl["DESI_TARGET"], sky_mask)
        )
    )[0]

    print("  {} targets are not one of the 3 recognized types".format(len(sky_unknown_rows)))

    if os.path.isfile(sky_file):
        os.remove(sky_file)
    with fitsio.FITS(sky_file, "rw") as fd:
        fd.write(sky_mtl)
        
def get_mtlstats(indir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/',passes=[]):
    science_file = indir + 'mtl_science.fits'
    ff = fitsio.read(science_file)
    types = ['LRG','ELG','QSO']
    print('after 4 passes:')
    for type in types:
        wt = (ff['DESI_TARGET'] & desi_mask[type]) > 0
        ntar = len(ff[wt])
        wtz = wt & (ff['NUMOBS_MORE'] == 0)
        nass = len(ff[wtz])
        print(type + ' total number of targets: '+str(ntar)+' , number with nobs =0 '+str(nass))

    for i in range(0,len(passes)):
        ps = passes[i]
        science_file = indir + 'mtl_science_pass'+str(passes[i+1])+'.fits'
        ff = fitsio.read(science_file)
        types = ['LRG','ELG','QSO']
        print('after '+str(ps)+' passes:')
        for type in types:
            wt = (ff['DESI_TARGET'] & desi_mask[type]) > 0
            ntar = len(ff[wt])
            wtz = wt & (ff['NUMOBS_MORE'] == 0)
            nass = len(ff[wtz])
            print(type + ' total number of targets: '+str(ntar)+' , number with nobs =0 '+str(nass))

def get_graystats(indir='/global/cscratch1/sd/ajross/fiberassigntest/fiducialtargets/temp/'):
    science_file = indir + 'mtl_science.fits'
    ff = fitsio.read(science_file)
    wc = ff['OBSCONDITIONS'] == 2
    print('there were '+str(len(ff[wc]))+' target assigned to be gray time only')
    print(np.unique(ff[wc]['DESI_TARGET']))
    print(np.unique(ff['DESI_TARGET']))
    we = ((ff['DESI_TARGET'] & desi_mask['ELG']) > 0)
    print(str(len(ff[we]))+' are ELG targets')
    wce = wc & ((ff['DESI_TARGET'] & desi_mask['ELG']) > 0)
    print(str(len(ff[wce]))+' are ELG targets')
    wcea = wce & (ff['NUMOBS_MORE'] == 0)
    print(str(len(ff[wcea]))+' were assigned')
        
