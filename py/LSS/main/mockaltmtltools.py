##from desiutil.iers import freeze_iers
##freeze_iers()

import collections.abc
from time import time
import astropy
import astropy.io
import astropy.io.fits as pf
from astropy.table import Table,join

import memory_profiler
from memory_profiler import profile
import gc

#TEMP
#MODULE_PATH = '/global/homes/a/acarnero/.local/lib/python3.10/site-packages/desitarget/__init__.py'
#MODULE_NAME = 'desitarget'
#import importlib
#import sys
#spec = importlib.util.spec_from_file_location(MODULE_NAME, MODULE_PATH)
#module = importlib.util.module_from_spec(spec)
#sys.modules[spec.name] = module
#spec.loader.exec_module(module)
#

import desitarget
from desitarget import io, mtl
from desitarget.cuts import random_fraction_of_trues
from desitarget.mtl import get_mtl_dir, get_mtl_tile_file_name,get_mtl_ledger_format
from desitarget.mtl import get_zcat_dir, get_ztile_file_name, tiles_to_be_processed
from desitarget.mtl import make_zcat,survey_data_model,update_ledger, get_utc_date

from desitarget.targets import initial_priority_numobs, decode_targetid
from desitarget.targetmask import obsconditions, obsmask
from desitarget.targetmask import desi_mask, bgs_mask, mws_mask, zwarn_mask

from desiutil.log import get_logger

import fitsio

import healpy as hp


from LSS.bitweights import pack_bitweights
from LSS.SV3.fatools import get_fba_fromnewmtl
import LSS.SV3.fatools as fatools

import matplotlib.pyplot as plt

import numpy as np
from numpy import random as rand
import numpy.lib.recfunctions as rfn

import os
import pickle
import subprocess
import sys
from time import sleep

import cProfile, pstats
import io as ProfileIO
from pstats import SortKey

import glob


pr = cProfile.Profile()

log = get_logger()

os.environ['DESIMODEL'] = '/global/common/software/desi/cori/desiconda/current/code/desimodel/master'

mtlformatdict = {"PARALLAX": '%16.8f', 'PMRA': '%16.8f', 'PMDEC': '%16.8f'}


zcatdatamodel = np.array([], dtype=[
    ('RA', '>f8'), ('DEC', '>f8'), ('TARGETID', '>i8'),
    ('NUMOBS', '>i4'), ('Z', '>f8'), ('ZWARN', '>i8'), ('ZTILEID', '>i4')
    ])

#mtltilefiledm = np.array([], dtype=[
#    ('TILEID', '>i4'), ('TIMESTAMP', 'U25'),
#    ('VERSION', 'U14'), ('PROGRAM', 'U6'), ('ZDATE', 'U8')
#    ])

mtltilefiledm = np.array([], dtype = [
    ('TILEID', '>i4'), ('TIMESTAMP', '<U25'),
    ('VERSION', '<U14'), ('PROGRAM', '<U6'), 
    ('ZDATE', '>i8'), ('ARCHIVEDATE', '>i8')])

def datesInMonthForYear(yyyy):
    # if divisible by 4
    if (yyyy % 4) == 0:
        # if not divisible by 100, leap year
        if not ((yyyy % 100) == 0):
            monthLengths = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        # if divisible by 100 and 400, leap year
        elif ((yyyy % 400) == 0):
            monthLengths = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        # if divisble by 100 and not 400, no leap year
        else:
            monthLengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    else:
        monthLengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    return monthLengths

def nextDate(date):
    # JL  takes NITE in YYYYMMDD form and increments to the next date
    yyyy, mm, dd = int(str(date)[0:4]), int(str(date)[4:6]), int(str(date)[6:])
    log.info('date = {0}'.format(date))
    monthLengths = datesInMonthForYear(yyyy)
    log.info('monthLengths array is {0}'.format(monthLengths))
    log.info('yyyy, mm, dd = {0}, {1}, {2}'.format(yyyy, mm, dd))
    if dd == monthLengths[mm - 1]:
        if mm == 12:
            mm = '01'
            yyyy = str(yyyy+1)
        else:
            mm = str(mm+1).zfill(2)
        
        dd = '01'
    else:
        dd = str(dd + 1).zfill(2)
    log.info('yyyy, mm, dd = {0}, {1}, {2}'.format(yyyy, mm, dd))
    return ''.join([str(yyyy), str(mm).zfill(2), str(dd).zfill(2)])

def evaluateMask(bits, mask, evalMultipleBits = False):
    if evalMultipleBits:
        return (bits & mask) == mask 
    return (bits & mask) > 0



def flipBit(cat, bit2Flip, cond = None, fieldName = 'DESI_TARGET', mode = 'on'):
    #only works on single bits
    assert( np.abs( np.log2(bit2Flip) - int(np.log2(bit2Flip)) ) < 0.001  )

    if cond is None:
        if mode.lower() == 'on':
            cond = np.invert(evaluateMask(cat[fieldName], bit2Flip))
        elif mode.lower() == 'off':
            cond = evaluateMask(cat[fieldName], bit2Flip)
        #elif mode.lower() == 'both':
        #    cond = np.ones(cat.shape[0], dtype = bool)
        else:
            #raise ValueError('`mode` must be `on` `off` or `both`')
            raise ValueError('`mode` must be `on` or `off`')

    assert( len(cond) ==  len(cat))

    if np.sum(cond) == 0:
        log.warning('This call to flipBit does not flip any bits.')
        return cat


    if mode == 'on':
        cat[fieldName][cond] = cat[fieldName][cond] | bit2Flip
    elif mode == 'off':
        cond = cond & ((cat[fieldName] & bit2Flip) == bit2Flip)
        cat[fieldName][cond] = cat[fieldName][cond] ^ bit2Flip
    #elif mode == 'both':
    #    cat[fieldName][cond] = cat[fieldName][cond] ^ bit2Flip
    else:
        #raise ValueError('`mode` must be `on` `off` or `both`')
        raise ValueError('`mode` must be `on` or `off`')

    return cat

def processTileFile(infile, outfile, startDate, endDate):
    #ztilefile, outputMTLDir + ztilefn, startDate, endDate
    if (startDate is None) and (endDate is None):
        #os.symlink(infile, outfile)
        from shutil import copyfile
        copyfile(infile, outfile)
        return 0
    
        
    if (startDate is None) or (startDate == ''):
        startDate = 0
    else:
        startDate = int(startDate.split('T')[0].replace('-', ''))       
    if (endDate is None) or (endDate == ''):
        endDate = 9999999999
    else:
        endDate = int(endDate.split('T')[0].replace('-', ''))

    origtf = Table.read(infile)

    origtf = origtf[origtf['LASTNIGHT'].astype(int) >= startDate ]
    origtf = origtf[origtf['LASTNIGHT'].astype(int) <= endDate ]


    origtf.write(outfile, overwrite = True, format = 'ascii.ecsv')
    return 0
def uniqueTimestampFATimePairs(tileList, withFlag = False):
    output = []
    for t in tileList: 
        if withFlag:
            datepair = (t['ORIGMTLTIMESTAMP'], t['FAMTLTIME'], t['REPROCFLAG'])
        else:
            datepair = (t['ORIGMTLTIMESTAMP'], t['FAMTLTIME'])

        if datepair in  output:
            continue
        else:
            output.append(datepair)

    return output
def uniqueArchiveDateZDatePairs(tileList, withFlag = False):
    output = []
    for t in tileList: 
        if withFlag:
            datepair = (t['ZDATE'], t['ARCHIVEDATE'], t['REPROCFLAG'])
        else:
            datepair = (t['ZDATE'], t['ARCHIVEDATE'])

        if datepair in  output:
            continue
        else:
            output.append(datepair)

    return output

def findTwin(altFiber, origFiberList, survey = 'sv3', obscon = 'dark'):
    log.critical('this function isn\'t ready yet. Goodbye')
    raise NotImplementedError('Fiber Twin method not implemented yet.')
    if survey == 'sv3':
        if obscon == 'dark':
            altTargBits = altFiber['SV3_DESI_TARGET']
            altTargBitsSec = altFiber['SV3_BGS_TARGET']
            altTargBitsMWS = altFiber['SV3_MWS_TARGET']

            origTargBitList = origFiberList['SV3_DESI_TARGET']
            origTargBitListSec = origFiberList['SV3_BGS_TARGET']
            origTargBitListMWS = origFiberList['SV3_MWS_TARGET']

        elif obscon == 'bright':
            altTargBits = altFiber['SV3_BGS_TARGET']
            altTargBitsSec = altFiber['SV3_DESI_TARGET']
            altTargBitsMWS = altFiber['SV3_MWS_TARGET']

            origTargBitList = origFiberList['SV3_BGS_TARGET']
            origTargBitListSec = origFiberList['SV3_DESI_TARGET']
            origTargBitListMWS = origFiberList['SV3_MWS_TARGET']
        else:
            raise ValueError('Invalid value for \'obscon\': {0}'.format(obscon))
    elif survey == 'main': 
        if obscon == 'dark':
            altTargBits = altFiber['DESI_TARGET']
            altTargBitsSec = altFiber['BGS_TARGET']
            altTargBitsMWS = altFiber['MWS_TARGET']

            origTargBitList = origFiberList['DESI_TARGET']
            origTargBitListSec = origFiberList['BGS_TARGET']
            origTargBitListMWS = origFiberList['MWS_TARGET']

        elif obscon == 'bright':
            altTargBits = altFiber['BGS_TARGET']
            altTargBitsSec = altFiber['DESI_TARGET']
            altTargBitsMWS = altFiber['MWS_TARGET']
            origTargBitList = origFiberList['BGS_TARGET']
            origTargBitListSec = origFiberList['DESI_TARGET']
            origTargBitListMWS = origFiberList['MWS_TARGET']

        else:
            raise ValueError('Invalid value for \'obscon\': {0}'.format(obscon))
    else:
        raise ValueError('Invalid value for \'survey\': {0}'.format(survey))

    altFS = altFiber['FIBERSTATUS']
    origFS = origFiberList['FIBERSTATUS']


    '''
    BGSBits = initialentries['SV3_BGS_TARGET']
    BGSFaintHIP = ((BGSBits & 8) == 8)
    BGSFaintAll = ((BGSBits & 1) == 1) | BGSFaintHIP

    #Set all BGS_FAINT_HIP to BGS_FAINT

    initialentries['SV3_BGS_TARGET'][BGSFaintHIP] = (BGSBits[BGSFaintHIP] & ~8)
    initialentries['PRIORITY'][BGSFaintHIP] = 102000*np.ones(np.sum(BGSFaintHIP))

    NewBGSBits = initialentries['SV3_BGS_TARGET']
    NewBGSFaintHIP = ((BGSBits & 8) == 8)
    NewBGSFaintAll = ((BGSBits & 1) == 1) | NewBGSFaintHIP
    NewBGSPriors = initialentries['PRIORITY']
    #Select 20% of BGS_FAINT to promote using function from 
    BGSFaintNewHIP = random_fraction_of_trues(PromoteFracBGSFaint, BGSFaintAll)
    #Promote them

    initialentries['SV3_BGS_TARGET'][BGSFaintNewHIP] = (BGSBits[BGSFaintNewHIP] | 8)
    initialentries['PRIORITY'][BGSFaintNewHIP] = 102100*np.ones(np.sum(BGSFaintNewHIP)).astype(int)
    '''


def createFAmap(FAReal, FAAlt, TargAlt = None, changeFiberOpt = None, debug = False,
 verbose = False, mock = False, mockTrueZKey = None):
    # Options for 'changeFiberOpt':
    # None: do nothing different to version 1
    # AllTwins: Find a twin fiber with a target of the 
    #   same type and similar Fiber assignment for all
    #   unsimilar target types
    # SomeTwins: Find a twin as above but only for 
    #   assignments where the original fiber was unassigned

    TIDReal = FAReal['TARGETID']
    TIDAlt = FAAlt['TARGETID']
    FibReal = FAReal['FIBER']
    FibAlt = FAAlt['FIBER']

    if not (changeFiberOpt is None):
        raise NotImplementedError('changeFiberOpt is not implemented yet.')
        assert(not(TargAlt is None))
        jTargs = join(FAAlt, TargAlt, keys = "TARGETID")
    
    Real2Alt = {}
    Alt2Real = {}
    if debug:
        inc1 = 0
        inc2 = 0
    negMisMatch = []
    for tr, fr in zip(TIDReal, FibReal):
        taMatch = TIDAlt[FibAlt == fr]
        assert(len(taMatch) == 1)
        if debug:
            try:
                assert(tr == taMatch[0])
            except:
                inc1+=1
        Real2Alt[tr] = taMatch[0]
    
    for ta, fa in zip(TIDAlt, FibAlt):
        trMatch = TIDReal[FibReal == fa]
        try:
            assert(len(trMatch) == 1)
        except:
            if ta < 0:
                negMisMatch.append(ta)
                continue
            else:
                log.info(ta)

                assert(0)
        if debug or verbose:
            try:
                assert(ta == trMatch[0])
            except:
                inc2+=1
        if (changeFiberOpt is None) or (changeFiberOpt == 'SomeTwins') or (ta == trMatch[0]):
            Alt2Real[ta] = trMatch[0]
        elif changeFiberOpt == 'AllTwins':
            #if jTargs['SV3_']
            assert(0)
            pass

    
    if debug or verbose:
        log.info('no matches for negative tas {0}'.format(negMisMatch))
        log.info(inc1)
        log.info(inc2)
    return Alt2Real, Real2Alt



def makeAlternateZCat(zcat, real2AltMap, alt2RealMap, debug = False, verbose = False):
    from collections import Counter
    zcatids = zcat['TARGETID']
    altZCat = Table(zcat)
    if debug:
        failures = 0
        negativeIDs = 0
    for n, i in zip(zcatids, range(len(zcatids))):
        cond = (n == zcatids)
        if debug and (n < 0):
            negativeIDs +=1   
        altid = real2AltMap[n]

        altZCat['TARGETID'][i] = altid
    if debug:
        log.info('negIDs')
        log.info(negativeIDs)
        log.info('failures')
        log.info(failures)
        log.info('testctr')
    d =  Counter(altZCat['TARGETID'])  
    res = [ k for k, v in d.items() if v > 1]
    if debug:
        log.info('res')
        log.info(res)
    if len(res):
        log.info('how many pre dup cuts')
        log.info(zcatids.shape)
        cond2 = np.ones(zcatids.shape, dtype=bool)
        for i in res:
            log.info('test')
            log.info(np.sum(zcatids == i))
            cond2 = cond2 & (altcatids != i)
        log.info("how many post dup cuts")
        log.info(np.sum(cond2))
    else:
        log.info("supposedly, no duplicates")
    return altZCat

def checkMTLChanged(MTLFile1, MTLFile2):
    MTL1 = desitarget.io.read_mtl_ledger(MTLFile1, unique = True)
    MTL2 = desitarget.io.read_mtl_ledger(MTLFile2, unique = True)
    NDiff = 0
    NDiff2 = 0
    NDiff3 = 0
    for tar1 in MTL1:
        tar2 = MTL2[MTL2['TARGETID'] == tar1['TARGETID']]

        if tar1['NUMOBS'] != tar2['NUMOBS']:
            NDiff +=1

        if tar1['TIMESTAMP'] != tar2['TIMESTAMP']:
            NDiff2 +=1
            
        if tar1['SUBPRIORITY'] != tar2['SUBPRIORITY']:
            NDiff3 +=1

    print('Number targets with different NUMOBS')
    print(NDiff)
    print('Number targets with different TIMESTAMP')
    print(NDiff2)
    print('Number targets with different SUBPRIORITY')
    print(NDiff3)

def makeTileTrackerFN(dirName, survey, obscon):
    return dirName + '/{0}survey-{1}obscon-TileTracker.ecsv'.format(survey, obscon.upper())
def makeTileTracker(altmtldir, survey = 'main', obscon = 'DARK', startDate = None,
    endDate = None, overwrite = True):
    """Create action file which orders all actions to do with AMTL in order 
    in which real survey did them.

    Parameters
    ----------
    altmtldir : :class:`str`
        Path to the directory for a single realization of alternate MTL
        ledgers. e.g. /pscratch/u/user/simName/Univ000/
    obscon : :class:`str`, optional, defaults to "dark"
        A string matching ONE obscondition in the desitarget bitmask yaml
        file (i.e. in `desitarget.targetmask.obsconditions`), e.g. "DARK"
        Governs how priorities are set when merging targets.
    survey : :class:`str`, optional, defaults to "main"
        Used to look up the correct ledger, in combination with `obscon`.
        Options are ``'main'`` and ``'svX``' (where X is 1, 2, 3 etc.)
        for the main survey and different iterations of SV, respectively.
    

    Returns
    -------
    
    [Nothing]

    Notes
    -----
    - Writes a tiletracker file to {altmtldir}/{survey.lower()}survey-{obscon.upper()}obscon-TileTracker.ecsv
    """

    TileTrackerFN = makeTileTrackerFN(altmtldir, survey, obscon)

    if (survey.lower() == 'main') or (survey.lower() == 'y1'):
        
        surveyForTSS = 'main'
        if survey.lower() == 'y1':
            TileFN = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/tiles-{0}.fits'.format(obscon.upper())
        elif survey.lower() == 'da2':
            TileFN = '/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/tiles-{0}.fits'.format(obscon.upper())
        else:
            TileFN = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-main.ecsv'
    elif survey.lower() == 'sv3':
        surveyForTSS = 'sv3'
        TileFN = '/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/tiles-{0}.fits'.format(obscon.upper())
    else:
        raise ValueError('only valid values for `survey` are `main` and `sv3.` {0} was provided'.format(survey))

    FABaseDir = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'

    Tiles = Table.read(TileFN)

    TSSFN = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv'

    TSS = Table.read(TSSFN)

    MTLDTFN = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/mtl-done-tiles.ecsv'

    MTLDT = Table.read(MTLDTFN)

    #tiles-specstatus file filtered to only matching obscon and surveySURVEY FAPRGRM
    TSS_Sel = TSS[(TSS['SURVEY'] == surveyForTSS) & (TSS['FAPRGRM'] == obscon.lower())]
    
    TilesSel = np.unique(TSS_Sel['TILEID'])
    
    TileIDs = []
    TypeOfActions = []
    TimesOfActions = []
    doneFlag = []
    archiveDates = []
    
    for tileid in TilesSel:
        print('tileid = {0}'.format(tileid))
        
        ts = str(tileid).zfill(6)
        
        thisTileMTLDT = MTLDT[MTLDT['TILEID'] == tileid]
        
        if len(thisTileMTLDT) > 1:
            thisTileMTLDT.sort('TIMESTAMP')
        elif len(thisTileMTLDT) == 0:
            continue
        else:
            log.info(len(thisTileMTLDT))
            log.info(thisTileMTLDT['ARCHIVEDATE'])
            log.info(thisTileMTLDT['ARCHIVEDATE'][0])
            log.info(type(thisTileMTLDT['ARCHIVEDATE'][0]))
            if thisTileMTLDT['ARCHIVEDATE'][0] > int(endDate):
                continue
        reprocFlag = False
        thisFAFN = FABaseDir + f'/{ts[0:3]}/fiberassign-{ts}.fits'

        thisfhtOrig = fitsio.read_header(thisFAFN)
        thisfadate = thisfhtOrig['MTLTIME']
        thisfadate = desitarget.mtl.add_to_iso_date(thisfadate, 1)
        thisfanite = int(''.join(thisfadate.split('T')[0].split('-')))
        if thisfanite > endDate:
            continue
        
        TileIDs.append(tileid)
        TypeOfActions.append('fa')
        TimesOfActions.append(thisfadate)
        archiveDates.append(thisfanite)
        if thisfanite < startDate:
            doneFlag.append(True)
        else:
            doneFlag.append(False)
        
        for update in thisTileMTLDT:
            
                
            thisupdateTimestamp = update['TIMESTAMP']
            thisupdateNite = int(''.join(thisupdateTimestamp.split('T')[0].split('-')))
            if (thisupdateNite > endDate):
                continue
            
            TileIDs.append(tileid)
            if reprocFlag:
                TypeOfActions.append('reproc')
            else:
                TypeOfActions.append('update')
            TimesOfActions.append(thisupdateTimestamp)
            if (thisupdateNite < startDate):
                doneFlag.append(True)
            else:
                doneFlag.append(False)
            archiveDates.append(update['ARCHIVEDATE'])
            reprocFlag = True
    ActionList = [TileIDs, TypeOfActions, TimesOfActions, doneFlag, archiveDates]
    t = Table(ActionList,
           names=('TILEID', 'ACTIONTYPE', 'ACTIONTIME', 'DONEFLAG', 'ARCHIVEDATE'),
           meta={'Name': 'AltMTLTileTracker', 'StartDate': startDate, 'EndDate': endDate, 'amtldir':altmtldir})
    t.sort(['ACTIONTIME', 'ACTIONTYPE', 'TILEID'])
    
    t.write(TileTrackerFN, format='ascii.ecsv', overwrite = overwrite)




def trimToMTL(notMTL, MTL, debug = False, verbose = False):
    # JL trims a target file, which possesses all of the information in an MTL, down
    # JL to the columns allowed in the MTL data model. 
    allNames = notMTL.dtype.names
    MTLNames = MTL.dtype.names
    for n in allNames:
        if n in MTLNames:
            if debug:
                print('allowed')
                print(n)
            continue
        else:
            if debug:
                print('killed')
                print(n)
            notMTL = rfn.drop_fields(notMTL, n)
    return notMTL


#@profile
def initializeAlternateMTLs(initMTL, outputMTL, nAlt = 2, genSubset = None, seed = 314159, 
    obscon = 'DARK', survey = 'sv3', saveBackup = False, overwrite = False, startDate = None, endDate = None,
    ztilefile = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv', 
    hpnum = None, shuffleBrightPriorities = False, PromoteFracBGSFaint = 0.2, shuffleELGPriorities = False, 
    PromoteFracELG = 0.1, shuffleSubpriorities = True, reproducing = False, usetmp = False, 
    finalDir = None, profile = False, debug = False, verbose = False):
    if profile:
        pr.enable()
    if verbose or debug:
        log.info('starting initializeAltMTLs')

    if (shuffleSubpriorities ^ reproducing):
        pass
    else:
        log.critical('If you are not shuffling subpriorities, you MUST be in debug/reproduction mode.')
        raise ValueError('If you are not shuffling subpriorities, you MUST be in debug/reproduction mode.')

    if ('trunk' in outputMTL.lower()) or  ('ops' in outputMTL.lower()):
        raise ValueError("In order to prevent accidental overwriting of the real MTLs, please remove \'ops\' and \'trunk\' from your MTL output directory")

    if (not usetmp) or (usetmp and (outputMTL.startswith('/dev/shm/') or not(outputMTL.startswith('/tmp/')))):
        pass
    else:
        log.critical('You are trying to write to local tmp directories but \
            your write directory is not in local tmp (/dev/shm/ or /tmp/).')
        log.critical('directory name: {0}'.format(outputMTL))
        raise ValueError('usetmp set to True but output directory not in tmp. Output directory is {0}'.format(outputMTL))

        
    if debug:
        log.info('initMTL')
        log.info(initMTL)
    ztilefn = ztilefile.split('/')[-1]
    fn = initMTL.split('/')[-1]
    log.info('reading initial MTL(s)')
    allentries = Table.read(initMTL) 
    
    meta = allentries.meta
    if verbose or debug:
        log.info('MTL metadata')
        log.info(meta)
        log.info('initial MTL')
        log.info(initMTL)
        log.info('output MTL')
        log.info(outputMTL)
    
    if not ('Univ' in outputMTL):
        log.warning('Code currently relies on using Univ as realization delimiter. \
            Code may function improperly.')
    altmtldir = os.path.dirname(outputMTL).split('Univ')[0]
    origmtldir = os.path.dirname(initMTL).split(survey)[0]
    #zcatdir = os.path.dirname(ztilefile)

    if (startDate is None) or (startDate == ''):

        firstTS = allentries[0]["TIMESTAMP"] 
        initialentries = allentries[allentries["TIMESTAMP"] == firstTS]
        subpriorsInit = initialentries["SUBPRIORITY"]
        startDateShort = 19990101
    else:
        log.debug('startdate')
        log.debug(startDate)
        initialentries = allentries[allentries["TIMESTAMP"] <= startDate]
        subpriorsInit = initialentries["SUBPRIORITY"] 

        origmtltilefn = os.path.join(origmtldir, get_mtl_tile_file_name(secondary=False))
        altmtltilefn = os.path.join(altmtldir, get_mtl_tile_file_name(secondary=False))
        startDateShort = int(startDate.split('T')[0].replace('-', ''))
    if ('T' in endDate) & ('-' in endDate):
        endDateShort = int(endDate.split('T')[0].replace('-', '')) 
    else:
        endDateShort = int(endDate)

    if verbose or debug:
        log.info('generate subset? {0}'.format(genSubset))
    if not genSubset is None:
        if type(genSubset) == int:
            if debug:
                log.info('genSubset Int')
            iterloop = [genSubset]
        elif (type(genSubset) == list) or (type(genSubset) == np.ndarray):
            if debug:
                log.info('genSubset Arraylike')
            iterloop = genSubset
    else:
        if debug:
            log.info('genSubset None')
        iterloop = range(nAlt)
    if verbose or debug:
        log.info('starting iterloop')
    for n in iterloop:
        if verbose or debug:
            log.info('Realization {0:d}'.format(n))
        outputMTLDir = outputMTL.format(n)
        if verbose or debug:
            log.info('outputMTLDir')
            log.info(outputMTLDir)
        outfile = outputMTLDir +'/' + str(survey).lower() + '/' + str(obscon).lower() + '/' + str(fn)
        if verbose or debug:
            log.info('outfile')
            log.info(outfile)
        if os.path.exists(outfile):
            if overwrite: 
                if verbose or debug:
                    log.info('overwrite')
                os.remove(outfile)
            else:
                if verbose or debug:
                    log.info('continuing')
                continue
        if type(hpnum) == str:
            try: 
                hpnum = int(hpnum)
            except:
                log.info('hpnum is string but not integer. Value is {0}'.format(hpnum))
                raise ValueError('hpnum is string but not integer. Value is {0}'.format(hpnum))
            rand.seed(seed + hpnum + n)
        elif isinstance(hpnum, int) or isinstance(hpnum, np.int64):
            rand.seed(seed + hpnum + n)
        elif isinstance(hpnum, float) or isinstance(hpnum, np.float64):
            assert(np.abs(hpnum - int(hpnum)) < 0.01)
            rand.seed(seed + int(hpnum) + n)
        else:
            log.info('hpnum = {0}'.format(hpnum))
            log.info('type(hpnum) = {0}'.format(type(hpnum)))
            assert(0)
            rand.seed(seed + n)
        if verbose or debug:
            log.info('pre creating output dir')
        if not os.path.exists(outputMTLDir):
            os.makedirs(outputMTLDir)
        if not os.path.exists(finalDir.format(n)):
            os.makedirs(finalDir.format(n))
        if not os.path.isfile(finalDir.format(n) + '/' + ztilefn):
            processTileFile(ztilefile, outputMTLDir + ztilefn, startDate, endDate)
            #os.symlink(ztilefile, outputMTLDir + ztilefn)
        thisTileTrackerFN = makeTileTrackerFN(finalDir.format(n), survey, obscon)
        log.info('path to tiletracker = {0}'.format(thisTileTrackerFN))
        if not os.path.isfile(thisTileTrackerFN):
            makeTileTracker(finalDir.format(n), survey = survey, obscon = obscon,overwrite = False,
             startDate = startDateShort, endDate = endDateShort)
            #makeTileTracker(outputMTLDir, survey = survey, obscon = obscon,overwrite = False,
            #startDate = startDateShort, endDate = endDateShort)
        subpriors = initialentries['SUBPRIORITY']

        if (not reproducing) and shuffleSubpriorities:
            newSubpriors = rand.uniform(size = len(subpriors))
        else:
            newSubpriors = np.copy(subpriors)
        try:
            
            assert((np.std(subpriorsInit - newSubpriors) > 0.001) | (len(subpriors) < 2) | ((not shuffleSubpriorities) and reproducing) )
        except:
            log.warning('first shuffle failed')
            log.warning('size of initial subprior array')
            log.warning(len(subpriorsInit))

            newSubpriors = rand.uniform(size = len(subpriors))
            assert((np.std(subpriorsInit - newSubpriors) > 0.001) | (len(subpriors) < 2))

        initialentries['SUBPRIORITY'] = newSubpriors
        

        # add main priority values


        if  (obscon.lower() == 'bright') and (shuffleBrightPriorities):
            if (survey.lower() == 'sv3'):
                BGSHIPBit = 2**3
                BGSBit = 2**0
                BGSPriorityInit = 102000
                BGSHIPPriority = 102100

                BGSBits = initialentries['SV3_BGS_TARGET']
            elif (survey.lower() == 'main'):
                BGSHIPBit = 2**3
                BGSBit = 2**0
                BGSPriorityInit = 2000
                BGSHIPPriority = 2100
                BGSBits = initialentries['BGS_TARGET']
            else:
                raise ValueError('Survey.lower should be `sv3` or `main` but is instead {0:s}'.format(survey.lower()))
            BGSFaintHIP = ((BGSBits & BGSHIPBit) == BGSHIPBit)
            BGSFaintAll = ((BGSBits & BGSBit) == BGSBit) | BGSFaintHIP

            #Set all BGS_FAINT_HIP to BGS_FAINT

            initialentries['SV3_BGS_TARGET'][BGSFaintHIP] = (BGSBits[BGSFaintHIP] & ~BGSHIPBit)
            initialentries['PRIORITY'][BGSFaintHIP] = BGSPriorityInit*np.ones(np.sum(BGSFaintHIP))
            initialentries['TARGET_STATE'][BGSFaintHIP] = np.broadcast_to(np.array(['BGS_FAINT|UNOBS']), BGSFaintHIP.shape)

            #Select 20% of BGS_FAINT to promote using function from desitarget
            BGSFaintNewHIP = random_fraction_of_trues(PromoteFracBGSFaint, BGSFaintAll)
            #Promote them

            initialentries['SV3_BGS_TARGET'][BGSFaintNewHIP] = (BGSBits[BGSFaintNewHIP] | BGSHIPBit)
            initialentries['TARGET_STATE'][BGSFaintNewHIP] = np.broadcast_to(np.array(['BGS_FAINT_HIP|UNOBS']), BGSFaintNewHIP.shape)
            initialentries['PRIORITY'][BGSFaintNewHIP] = BGSHIPPriority*np.ones(np.sum(BGSFaintNewHIP)).astype(int)
            initialentries['PRIORITY_INIT'][BGSFaintNewHIP] = BGSHIPPriority*np.ones(np.sum(BGSFaintNewHIP)).astype(int)

        elif (survey.lower() == 'main') and (obscon.lower() == 'dark') and (shuffleELGPriorities):

            #desi_mask

            #evaluateMask(bit, mask, evalMultipleBits = False):
            #flipBit(cat, bit2Flip, cond = None, fieldName = 'DESI_TARGET', mode = 'on'):

            ELGBits = initialentries['DESI_TARGET']

            #Set up condition arrays to select each type of target class
            LRGs    = evaluateMask(ELGBits, desi_mask['LRG'])
            ELGs    = evaluateMask(ELGBits, desi_mask['ELG'])
            QSOs    = evaluateMask(ELGBits, desi_mask['QSO'])
            ELGHIPs = evaluateMask(ELGBits, desi_mask['ELG_HIP'])
            ELGLOPs = evaluateMask(ELGBits, desi_mask['ELG_LOP'])
            ELGVLOs = evaluateMask(ELGBits, desi_mask['ELG_VLO'])
            log.info('ELGs:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGs)))
            log.info('LRGs:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(LRGs)))
            log.info('QSOs:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(QSOs)))
            log.info('ELGHIPs:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPs)))
            log.info('ELGLOPs:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGLOPs)))
            log.info('ELGVLOs:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGVLOs)))


            log.info('ELGHIPs&ELGLOPs:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPs&ELGLOPs)))
            log.info('ELGHIPs&ELGVLOs:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPs&ELGVLOs)))
            log.info('ELGHIPs&LRGs:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPs&LRGs)))
            log.info('ELGHIPs&QSOs:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPs&QSOs)))
            log.info('ELGLOPs&LRGs:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGLOPs&LRGs)))
            log.info('ELGLOPs&QSOs:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGLOPs&QSOs)))
            log.info('ELGVLOs&LRGs:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGVLOs&LRGs)))
            log.info('ELGVLOs&QSOs:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGVLOs&QSOs)))


            #turn off the ELG_HIP bit
            initialentries = flipBit(initialentries, desi_mask['ELG_HIP'], cond = ELGHIPs, mode = 'off', fieldName = 'DESI_TARGET')

            #reset object priority, priority_init, and numobs_init based on new target bits. 
            outpriority, outnumobs = initial_priority_numobs(initialentries, obscon = 'DARK')
            initialentries['PRIORITY'][ELGHIPs] = outpriority[ELGHIPs]
            initialentries['PRIORITY_INIT'][ELGHIPs] = outpriority[ELGHIPs]
            initialentries['NUMOBS_INIT'][ELGHIPs] = outnumobs[ELGHIPs]


            #JL - reset TARGET_STATES based on new target bits. This step isn't necessary for AMTL function but makes debugging using target states vastly easier. 
            initialentries['TARGET_STATE'][ELGHIPs & ELGVLOs & np.invert(LRGs) & np.invert(QSOs)] = np.broadcast_to(np.array(['ELG_VLO|UNOBS']), np.sum(ELGHIPs & ELGVLOs & np.invert(LRGs) & np.invert(QSOs) ) )

            initialentries['TARGET_STATE'][ELGHIPs & ELGLOPs & np.invert(LRGs) & np.invert(QSOs)] = np.broadcast_to(np.array(['ELG_LOP|UNOBS']), np.sum(ELGHIPs & ELGLOPs & np.invert(LRGs) & np.invert(QSOs) ) )

            initialentries['TARGET_STATE'][ELGHIPs & LRGs] = np.broadcast_to(np.array(['LRG|UNOBS']), np.sum(ELGHIPs & LRGs) )


            #For Debug. New Target bit flags after demoting all ELG_HIPs
            ELGBitsMid = initialentries['DESI_TARGET']
            LRGsMid    = evaluateMask(ELGBitsMid, desi_mask['LRG'])
            ELGsMid    = evaluateMask(ELGBitsMid, desi_mask['ELG'])
            QSOsMid    = evaluateMask(ELGBitsMid, desi_mask['QSO'])
            ELGHIPsMid = evaluateMask(ELGBitsMid, desi_mask['ELG_HIP'])
            ELGLOPsMid = evaluateMask(ELGBitsMid, desi_mask['ELG_LOP'])
            ELGVLOsMid = evaluateMask(ELGBitsMid, desi_mask['ELG_VLO'])
            log.info('ELGsMid:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGsMid)))
            log.info('LRGsMid:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(LRGsMid)))
            log.info('QSOsMid:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(QSOsMid)))
            log.info('ELGHIPsMid:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPsMid)))
            log.info('ELGLOPsMid:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGLOPsMid)))
            log.info('ELGVLOsMid:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGVLOsMid)))


            log.info('ELGHIPsMid&ELGLOPsMid:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPsMid&ELGLOPsMid)))
            log.info('ELGHIPsMid&ELGVLOsMid:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPsMid&ELGVLOsMid)))
            log.info('ELGHIPsMid&LRGsMid:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPsMid&LRGsMid)))
            log.info('ELGHIPsMid&QSOsMid:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPsMid&QSOsMid)))
            log.info('ELGLOPsMid&LRGsMid:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGLOPsMid&LRGsMid)))
            log.info('ELGLOPsMid&QSOsMid:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGLOPsMid&QSOsMid)))
            log.info('ELGVLOsMid&LRGsMid:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGVLOsMid&LRGsMid)))
            log.info('ELGVLOsMid&QSOsMid:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGVLOsMid&QSOsMid)))



            #Determine which 10% of ELGLOP and ELGVLO will be promoted to ELGHIP. These are done separately.

            #ELGNewHIP = random_fraction_of_trues(PromoteFracELG, ELGLOPs)

            #ELGNewHIP = ELGNewHIP | random_fraction_of_trues(PromoteFracELG, ELGVLOs)

            chosenLOP = rand.random(len(ELGLOPs)) < 0.1
            ELGNewHIP_FromLOP = ELGLOPs & chosenLOP 

            chosenVLO = rand.random(len(ELGVLOs)) < 0.1
            ELGNewHIP_FromVLO = ELGVLOs & chosenVLO

            ELGNewHIP = ELGNewHIP_FromLOP | ELGNewHIP_FromVLO

            log.info('ELGNewHIP_FromVLO:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGNewHIP_FromVLO)))
            log.info('ELGNewHIP_FromLOP:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGNewHIP_FromLOP)))
            log.info('ELGNewHIP:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGNewHIP)))

            #promote the just-determined 10% of ELG_LOP/ELG_VLO
            initialentries = flipBit(initialentries, desi_mask['ELG_HIP'], cond = ELGNewHIP, mode = 'on', fieldName = 'DESI_TARGET')


            #For Debug. New Target bit flags after promoting 10% of ELGs to HIP
            ELGBitsFinal = initialentries['DESI_TARGET']
            LRGsFinal    = evaluateMask(ELGBitsFinal, desi_mask['LRG'])
            ELGsFinal    = evaluateMask(ELGBitsFinal, desi_mask['ELG'])
            QSOsFinal    = evaluateMask(ELGBitsFinal, desi_mask['QSO'])
            ELGHIPsFinal = evaluateMask(ELGBitsFinal, desi_mask['ELG_HIP'])
            ELGLOPsFinal = evaluateMask(ELGBitsFinal, desi_mask['ELG_LOP'])
            ELGVLOsFinal = evaluateMask(ELGBitsFinal, desi_mask['ELG_VLO'])
            log.info('ELGsFinal:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGsFinal)))
            log.info('LRGsFinal:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(LRGsFinal)))
            log.info('QSOsFinal:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(QSOsFinal)))
            log.info('ELGHIPsFinal:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPsFinal)))
            log.info('ELGLOPsFinal:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGLOPsFinal)))
            log.info('ELGVLOsFinal:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGVLOsFinal)))


            log.info('ELGHIPsFinal&ELGLOPsFinal:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPsFinal&ELGLOPsFinal)))
            log.info('ELGHIPsFinal&ELGVLOsFinal:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPsFinal&ELGVLOsFinal)))
            log.info('ELGHIPsFinal&LRGsFinal:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPsFinal&LRGsFinal)))
            log.info('ELGHIPsFinal&QSOsFinal:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGHIPsFinal&QSOsFinal)))
            log.info('ELGLOPsFinal&LRGsFinal:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGLOPsFinal&LRGsFinal)))
            log.info('ELGLOPsFinal&QSOsFinal:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGLOPsFinal&QSOsFinal)))
            log.info('ELGVLOsFinal&LRGsFinal:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGVLOsFinal&LRGsFinal)))
            log.info('ELGVLOsFinal&QSOsFinal:HPNUM:{0}:Sum:{1}'.format(hpnum, np.sum(ELGVLOsFinal&QSOsFinal)))



            #reset object priority, priority_init, and numobs_init based on new target bits. 
            outpriority, outnumobs = initial_priority_numobs(initialentries, obscon = 'DARK')
            initialentries['PRIORITY'][ELGNewHIP] = outpriority[ELGNewHIP]
            initialentries['PRIORITY_INIT'][ELGNewHIP] = outpriority[ELGNewHIP]
            initialentries['NUMOBS_INIT'][ELGNewHIP] = outnumobs[ELGNewHIP]

            #JL - reset TARGET_STATES based on new target bits. This step isn't necessary for AMTL function but makes debugging using target states vastly easier. 
            initialentries['TARGET_STATE'][ELGNewHIP & np.invert(QSOs)] = np.broadcast_to(np.array(['ELG_HIP|UNOBS']), np.sum(ELGNewHIP & np.invert(QSOs)  ) )

        retval = desitarget.io.write_mtl(outputMTLDir, initialentries, survey=survey, obscon=obscon, extra=meta, nsidefile=meta['FILENSID'], hpxlist = [meta['FILEHPX']])
        if debug or verbose:
            log.info('(nowrite = False) ntargs, fn = {0}'.format(retval))
        log.info('wrote MTLs to {0}'.format(outputMTLDir))
        if saveBackup and (not usetmp):
            if not os.path.exists(str(outputMTLDir) +'/' + str(survey).lower() + '/' +str(obscon).lower() + '/orig/'):
                os.makedirs(str(outputMTLDir) +'/' + str(survey).lower() + '/' +str(obscon).lower() + '/orig/')
            
            
            if not os.path.exists(str(outputMTLDir) +'/' + str(survey).lower() + '/' +str(obscon).lower() + '/orig/' + str(fn)):
                from shutil import copyfile
                copyfile(str(outputMTLDir) +'/' + str(survey).lower() + '/' + str(obscon).lower() + '/' + str(fn), str(outputMTLDir) +'/' + str(survey).lower() + '/' +str(obscon).lower() + '/orig/' + str(fn))
        if usetmp:
            from shutil import copyfile

            if not os.path.exists(str(finalDir.format(n)) +'/' + str(survey).lower() + '/' +str(obscon).lower() ):
                os.makedirs(str(finalDir.format(n)) +'/' + str(survey).lower() + '/' +str(obscon).lower() )
            if saveBackup and (not os.path.exists(str(finalDir.format(n)) +'/' + str(survey).lower() + '/' +str(obscon).lower() + '/orig/')):
                os.makedirs(str(finalDir.format(n)) +'/' + str(survey).lower() + '/' +str(obscon).lower() + '/orig/')
            if debug:
                log.info('tempdir contents before copying')
                log.info(glob.glob(outputMTLDir + '/*' ))
                log.info(glob.glob(outputMTLDir + '/main/dark/*' ))
            copyfile(str(outputMTLDir) +'/' + str(survey).lower() + '/' + str(obscon).lower() + '/' + str(fn), str(finalDir.format(n)) +'/' + str(survey).lower() + '/' +str(obscon).lower() + '/' + str(fn))
            if debug:
                log.info('tempdir contents after copying')
                log.info(glob.glob(outputMTLDir + '/*' ))
                log.info(glob.glob(outputMTLDir + '/main/dark/*' ))

            if saveBackup and not os.path.exists(str(outputMTLDir) +'/' + str(survey).lower() + '/' +str(obscon).lower() + '/orig/' + str(fn)):
                #JL Potentially move the saveBackup copying to an afterburner
                #JL to speed up afterburner process. Copy all at once
                copyfile(str(outputMTLDir) +'/' + str(survey).lower() + '/' + str(obscon).lower() + '/' + str(fn), str(finalDir.format(n)) +'/' + str(survey).lower() + '/' +str(obscon).lower() + '/orig/' + str(fn))
                
            os.remove(str(outputMTLDir) +'/' + str(survey).lower() + '/' + str(obscon).lower() + '/' + str(fn))
            if debug:
                log.info('tempdir contents after removing')
                log.info(glob.glob(outputMTLDir + '/*' ))
    if usetmp:
        
        if verbose or debug:
            log.info('cleaning up tmpdir')
            log.info(glob.glob(outputMTLDir + '*' ))
        f2c = glob.glob(outputMTLDir + '*' )
        if verbose or debug:
            log.info('finaldir')
            log.info(finalDir.format(n))
        for tempfn in f2c:
            if '.' in str(os.path.split(tempfn)[1]):
                if verbose or debug:
                    log.info('copying tempfn: {0}'.format(tempfn))
                copyfile(tempfn , str(finalDir.format(n)) +'/' + os.path.basename(tempfn) )

        if verbose or debug:
            log.info('tempdir contents after copying')
            log.info(glob.glob(outputMTLDir + '*' ))

    if profile:
        pr.disable()
        s = ProfileIO.StringIO()
        sortby = SortKey.CUMULATIVE
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        if usetmp:

            ps.dump_stats(str(finalDir.format(n)) +'/' + str(survey).lower() + '/' + str(obscon).lower() + '/' + str(fn) + '.prof')
        else:
            ps.dump_stats(str(outputMTLDir) +'/' + str(survey).lower() + '/' + str(obscon).lower() + '/' + str(fn) + '.prof')
        print(s.getvalue())
        
        

def quickRestartFxn(ndirs = 1, altmtlbasedir = None, survey = 'sv3', obscon = 'dark', multiproc =False, nproc = None, verbose = False, debug = False):
    if verbose or debug:
        log.info('quick restart running')
    from shutil import copyfile, move
    from glob import glob as ls
    if multiproc:
        iterloop = range(nproc, nproc+1)
    else:
        iterloop = range(ndirs)
    for nRestart in iterloop:
        if verbose or debug:
            log.info(nRestart)
        altmtldirRestart = altmtlbasedir + '/Univ{0:03d}/'.format(nRestart)
        if os.path.exists(altmtldirRestart + 'mtl-done-tiles.ecsv'):
            move(altmtldirRestart + 'mtl-done-tiles.ecsv',altmtldirRestart + 'mtl-done-tiles.ecsv.old')
        restartMTLs = ls(altmtldirRestart +'/' + survey + '/' + obscon + '/' + '/orig/*')
        for fn in restartMTLs:
            copyfile(fn, altmtldirRestart +'/' + survey + '/' + obscon + '/' + fn.split('/')[-1])

def do_fiberassignment(altmtldir, FATiles, survey = 'sv3', obscon = 'dark', 
    verbose = False, debug = False, getosubp = False, redoFA = False, mock = False, reproducing = False):
    #FATiles = tiles_to_be_processed_alt(altmtldir, obscon = obscon, survey = survey, today = today, mode = 'fa')
    if len(FATiles):
        try:
            log.info('FATiles[0] = {0}'.format(FATiles[0]))
            if isinstance(FATiles[0], (collections.abc.Sequence, np.ndarray)):
                pass 
            else:
                FATiles = [FATiles]
        except:
            log.info('cannot access element 0 of FATiles')
    log.info('FATiles = {0}'.format(FATiles))


    OrigFAs = []
    AltFAs = []
    AltFAs2 = []
    TSs = []
    fadates = []


    #if len(FATiles):
    #    log.info('len FATiles = {0}'.format(len(FATiles)))
    #    pass 
    #else:
    #    return OrigFAs, AltFAs, AltFAs2, TSs, fadates, FATiles
    for t in FATiles:
        log.info('t = {0}'.format(t))
        #JL This loop takes each of the original fiberassignments for each of the tiles on $date
        #JL and opens them to obtain information for the alternative fiber assignments.
        #JL Then it runs the alternative fiber assignments, stores the results in an array (AltFAs)
        #JL while also storing the original fiber assignment files in a different array (OrigFA)

        ts = str(t['TILEID']).zfill(6)
        #JL Full path to the original fiber assignment from the real survey
        FAOrigName = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
        fhtOrig = fitsio.read_header(FAOrigName)
        fadate = fhtOrig['RUNDATE']
        # e.g. DESIROOT/target/catalogs/dr9/1.0.0/targets/main/resolve/dark
        targver = fhtOrig['TARG'].split('/targets')[0].split('/')[-1]
        assert(not ('/' in targver))
        log.info('fadate = {0}'.format(fadate))
        #JL stripping out the time of fiber assignment to leave only the date
        #JL THIS SHOULD ONLY BE USED IN DIRECTORY NAMES. THE ACTUAL RUNDATE VALUE SHOULD INCLUDE A TIME
        fadate = ''.join(fadate.split('T')[0].split('-'))
        log.info('fadate stripped = {0}'.format(fadate))
        fbadirbase = altmtldir + '/fa/' + survey.upper() +  '/' + fadate + '/'
        
        log.info('fbadirbase = {0}'.format(fbadirbase))
        log.info('ts = {0}'.format(ts))
        ##log.info('t[reprocflag] (should be false if here)= {0}'.format(t['REPROCFLAG']))
        ##assert(not bool(t['REPROCFLAG']))
        #if str(ts) == str(3414).zfill(6):
        #    raise ValueError('Not only do I create the backup here but I also need to fix the reproc flag')
        
        if getosubp:
            #JL When we are trying to reproduce a prior survey and/or debug, create a separate
            #JL directory in fbadirbase + /orig/ to store the reproduced FA files. 
            FAAltName = fbadirbase + '/orig/fba-' + ts+ '.fits'
            #FAMapName = fbadirbase + '/orig/famap-' + ts + '.pickle'
            fbadir = fbadirbase + '/orig/'
        else:

            #JL For normal "alternate" operations, store the fiber assignments
            #JL in the fbadirbase directory. 

            FAAltName = fbadirbase + '/fba-' + ts+ '.fits'
            #FAMapName = fbadirbase + '/famap-' + ts + '.pickle'
            fbadir = fbadirbase
        if verbose or debug:
            log.info('FAOrigName = {0}'.format(FAOrigName))

            log.info('FAAltName = {0}'.format(FAAltName))

        #JL Sometimes fiberassign leaves around temp files if a run is aborted. 
        #JL This command removes those temp files to prevent endless crashes. 
        if os.path.exists(FAAltName + '.tmp'):
            os.remove(FAAltName + '.tmp')
        #JL If the alternate fiberassignment was already performed, don't repeat it
        #JL Unless the 'redoFA' flag is set to true
        if verbose or debug:
            log.info('redoFA = {0}'.format(redoFA))
            log.info('FAAltName = {0}'.format(FAAltName))

        if  redoFA or (not os.path.exists(FAAltName)):
            if verbose and os.path.exists(FAAltName):
                log.info('repeating fiberassignment')
            elif verbose:
                log.info('fiberassignment not found, running fiberassignment')
            if verbose:
                log.info(ts)
                log.info(altmtldir + survey.lower())
                log.info(fbadir)
                log.info(getosubp)
                log.info(redoFA)
            if getosubp and verbose:
                log.info('checking contents of fiberassign directory before calling get_fba_from_newmtl')
                log.info(glob.glob(fbadir + '/*' ))
            #get_fba_fromnewmtl(ts,mtldir=altmtldir + survey.lower() + '/',outdir=fbadirbase, getosubp = getosubp, overwriteFA = redoFA, verbose = verbose, mock = mock, targver = targver)#, targets = targets)
            get_fba_fromnewmtl(ts,mtldir=altmtldir + survey.lower() + '/',outdir=fbadirbase, getosubp = getosubp, overwriteFA = redoFA, verbose = verbose, mock = mock, targver = targver, reproducing = reproducing)#, targets = targets)
            command_run = (['bash', fbadir + 'fa-' + ts + '.sh']) 
            if verbose:
                log.info('fa command_run')
                log.info(command_run)
            result = subprocess.run(command_run, capture_output = True)
        else: 
            log.info('not repeating fiberassignment')
        log.info('adding fiberassignments to arrays')
        OrigFAs.append(pf.open(FAOrigName)[1].data)
        AltFAs.append(pf.open(FAAltName)[1].data)
        AltFAs2.append(pf.open(FAAltName)[2].data)
        TSs.append(ts)
        fadates.append(fadate)
        
    return OrigFAs, AltFAs, AltFAs2, TSs, fadates, FATiles

def make_fibermaps(altmtldir, OrigFAs, AltFAs, AltFAs2, TSs, fadates, tiles, survey = 'sv3', obscon = 'dark', changeFiberOpt = None, verbose = False, debug = False, getosubp = False, redoFA = False):
    A2RMap = {}
    R2AMap = {}
    if verbose:
        log.info('beginning loop through FA files')
    for ofa, afa, afa2, ts, fadate, t in zip(OrigFAs, AltFAs, AltFAs2, TSs, fadates, tiles):
        log.info('ts = {0}'.format(ts))
        if changeFiberOpt is None:
            A2RMap, R2AMap = createFAmap(ofa, afa, changeFiberOpt = changeFiberOpt)
        else:
            raise NotImplementedError('changeFiberOpt has not yet been implemented')

            #FAOrigName = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
            
            A2RMap, R2AMap = createFAmap(ofa, afa, TargAlt = afa2, changeFiberOpt = changeFiberOpt)

        fbadirbase = altmtldir + '/fa/' + survey.upper() +  '/' + fadate + '/'
        if getosubp:
            FAAltName = altmtldir + '/fa/' + survey.upper() +  '/' + fadate + '/orig/fba-' + ts+ '.fits'
            FAMapName = fbadirbase + '/orig/famap-' + ts + '.pickle'
            fbadir = altmtldir + '/fa/' + survey.upper() +  '/' + fadate + '/orig/'
        else:

            FAAltName = altmtldir + '/fa/' + survey.upper() +  '/' + fadate + '/fba-' + ts+ '.fits'
            FAMapName = fbadirbase + '/famap-' + ts + '.pickle'
            fbadir = fbadirbase


        if debug:
            log.info('ts = {0}'.format(ts))
            log.info('FAMapName = {0}'.format(FAMapName))
        
        
        if redoFA or (not (os.path.isfile(FAMapName))):
            if verbose:
                log.info('dumping out fiber map to pickle file')
            with open(FAMapName, 'wb') as handle:
                pickle.dump((A2RMap, R2AMap), handle, protocol=pickle.HIGHEST_PROTOCOL)
        #thisUTCDate = get_utc_date(survey=survey)
        if verbose:
            log.info('---')
            log.info('unique keys in R2AMap = {0:d}'.format(np.unique(R2AMap.keys()).shape[0]))
            log.info('---')

            log.info('---')
            log.info('unique keys in A2RMap = {0:d}'.format(np.unique(A2RMap.keys()).shape[0]))
            log.info('---')
        #retval = write_amtl_tile_tracker(altmtldir, [t], obscon = obscon, survey = survey, mode = 'fa')
        retval = write_amtl_tile_tracker(altmtldir, [t], obscon = obscon, survey = survey)
        log.info('write_amtl_tile_tracker retval = {0}'.format(retval))

    return A2RMap, R2AMap
def update_alt_ledger(altmtldir,althpdirname, altmtltilefn,  actions, survey = 'sv3', obscon = 'dark', today = None, 
    getosubp = False, zcatdir = None, mock = False, numobs_from_ledger = True, targets = None, verbose = False, debug = False):
    if verbose or debug:
        log.info('today = {0}'.format(today))
        log.info('obscon = {0}'.format(obscon))
        log.info('survey = {0}'.format(survey))
    #UpdateTiles = tiles_to_be_processed_alt(altmtldir, obscon = obscon, survey = survey, today = today, mode = 'update')
    #log.info('updatetiles = {0}'.format(UpdateTiles))
    # ADM grab the zcat directory (in case we're relying on $ZCAT_DIR).
    zcatdir = get_zcat_dir(zcatdir)
    # ADM And contruct the associated ZTILE filename.
    ztilefn = os.path.join(zcatdir, get_ztile_file_name())
    #if len(UpdateTiles):
    #    pass 
    #else:
    #    return althpdirname, altmtltilefn, ztilefn, None
    #isinstance(FATiles[0], (collections.abc.Sequence, np.ndarray))
    if not isinstance(actions['TILEID'], (collections.abc.Sequence, np.ndarray)):
        actions = [actions]
    log.info('actions = {0}'.format(actions))
    for t in actions:
        log.info('t = {0}'.format(t))
        if t['ACTIONTYPE'].lower() == 'reproc':
            raise ValueError('Reprocessing should be handled elsewhere.')
            #raise ValueError('Make sure backup is made and reprocessing logic is correct before beginning reprocessing.')
        ts = str(t['TILEID']).zfill(6)

        FAOrigName = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz'
        fhtOrig = fitsio.read_header(FAOrigName)
        fadate = fhtOrig['RUNDATE']
        fadate = ''.join(fadate.split('T')[0].split('-'))
        fbadirbase = altmtldir + '/fa/' + survey.upper() +  '/' + fadate + '/'
        log.info('t = {0}'.format(t))
        log.info('fbadirbase = {0}'.format(fbadirbase))
        log.info('ts = {0}'.format(ts))

        if getosubp:
            FAMapName = fbadirbase + '/orig/famap-' + ts + '.pickle'
        else:
            FAMapName = fbadirbase + '/famap-' + ts + '.pickle'

        log.info('FAMapName = {0}'.format(FAMapName))
        with open(FAMapName,'rb') as fl:
            (A2RMap, R2AMap) = pickle.load(fl,fix_imports = True)

        # ADM create the catalog of updated redshifts.
        log.info('making zcats')
        log.info('zcatdir = {0}'.format(zcatdir))
        log.info('t = {0}'.format(t))
        zcat = make_zcat(zcatdir, [t], obscon, survey)

        altZCat = makeAlternateZCat(zcat, R2AMap, A2RMap, debug = debug, verbose = verbose)
        # ADM insist that for an MTL loop with real observations, the zcat
        # ADM must conform to the data model. In particular, it must include
        # ADM ZTILEID, and other columns addes for the Main Survey. These
        # ADM columns may not be needed for non-ledger simulations.
        # ADM Note that the data model differs with survey type.
        #print('AURE')
        #print('zcat columns', zcat.columns)
        #print('altZCat columns', altZCat.columns)

        zcatdm = survey_data_model(zcatdatamodel, survey=survey)
        if zcat.dtype.descr != zcatdm.dtype.descr:
            msg = "altZCat data model must be {} not {}!".format(
            ##TEMPmsg = "zcat data model must be {} not {}!".format(
                zcatdm.dtype.descr, altZCat.dtype.descr)
            log.critical(msg)
        #    log.critical('Changing format')
            #print('BEFORE')
            #print(altZCat)
            #correct_dtype = np.dtype([('RA', '>f8'), ('DEC', '>f8'), ('TARGETID', '>i8'), ('NUMOBS', '>i4'),
            #                          ('Z', '>f8'), ('ZWARN', '>i8'), ('ZTILEID', '>i4'), ('Z_QN', '>f8'), 
            #                          ('IS_QSO_QN', '>i2'), ('DELTACHI2', '>f8')])
            #corrected_data = Table()
            #for colname in altZCat.colnames:
            #    corrected_data[colname] = altZCat[colname].astype(correct_dtype[colname])
            #print('AFTER')
            #print(corrected_data)
            #altZCat = corrected_data
        #if altZCat.dtype.descr != zcatdm.dtype.descr:
        #    msg = "altZCat data model must be {} not {}!".format(
        #        zcatdm.dtype.descr, altZCat.dtype.descr)
        #    log.critical(msg)

            raise ValueError(msg)
        #else:
        #    print('WE PASS IT!!')
        # ADM useful to know how many targets were updated.
        _, _, _, _, sky, _ = decode_targetid(zcat["TARGETID"])
        ntargs, nsky = np.sum(sky == 0), np.sum(sky)
        msg = "Update state for {} targets".format(ntargs)
        msg += " (the zcats also contain {} skies with +ve TARGETIDs)".format(nsky)
        log.info(msg)
        didUpdateHappen = False
        # ADM update the appropriate ledger.

        if targets is None:
            raise ValueError('If processing mocks, you MUST specify a target file')
        log.info('update loc a')
        update_ledger(althpdirname, altZCat, obscon=obscon.upper(),
                  numobs_from_ledger=numobs_from_ledger, tabform='ascii.ecsv')#, targets = targets)
        didUpdateHappen = True
        assert(didUpdateHappen)
        if verbose or debug:
            log.info('if main, should sleep 1 second')
        #thisUTCDate = get_utc_date(survey=survey)
        if survey == "main":
            sleep(1)
            if verbose or debug:
                log.info('has slept one second')
            #t["ALTARCHIVEDATE"] = thisUTCDate
        if verbose or debug:
            log.info('now writing to amtl_tile_tracker')
        #io.write_mtl_tile_file(altmtltilefn,dateTiles)
        #write_amtl_tile_tracker(altmtldir, dateTiles, thisUTCDate, obscon = obscon, survey = survey)
        log.info('changes are being registered')
        log.info('altmtldir = {0}'.format(altmtldir))
        log.info('t = {0}'.format(t))
        #log.info('thisUTCDate = {0}'.format(thisUTCDate))
        log.info('today = {0}'.format(today))
        #retval = write_amtl_tile_tracker(altmtldir, [t], obscon = obscon, survey = survey, mode = 'update')
        retval = write_amtl_tile_tracker(altmtldir, [t], obscon = obscon, survey = survey)
        log.info('write_amtl_tile_tracker retval = {0}'.format(retval))
        if verbose or debug:
            log.info('has written to amtl_tile_tracker')

    return althpdirname, altmtltilefn, ztilefn, actions
#@profile
def loop_alt_ledger(obscon, survey='sv3', zcatdir=None, mtldir=None,
                altmtlbasedir=None, ndirs = 3, numobs_from_ledger=True, 
                secondary=False, singletile = None, singleDate = None, debugOrig = False, 
                    getosubp = False, quickRestart = False, redoFA = False,
                    multiproc = False, nproc = None, testDoubleDate = False, 
                    changeFiberOpt = None, targets = None, mock = False,
                    debug = False, verbose = False, reproducing = False):
    """Execute full MTL loop, including reading files, updating ledgers.

    Parameters
    ----------
    obscon : :class:`str`
        A string matching ONE obscondition in the desitarget bitmask yaml
        file (i.e. in `desitarget.targetmask.obsconditions`), e.g. "DARK"
        Governs how priorities are set when merging targets.
    survey : :class:`str`, optional, defaults to "main"
        Used to look up the correct ledger, in combination with `obscon`.
        Options are ``'main'`` and ``'svX``' (where X is 1, 2, 3 etc.)
        for the main survey and different iterations of SV, respectively.
    zcatdir : :class:`str`, optional, defaults to ``None``
        Full path to the "daily" directory that hosts redshift catalogs.
        If this is ``None``, look up the redshift catalog directory from
        the $ZCAT_DIR environment variable.
    mtldir : :class:`str`, optional, defaults to ``None``
        Full path to the directory that hosts the MTL ledgers and the MTL
        tile file. If ``None``, then look up the MTL directory from the
        $MTL_DIR environment variable.
    altmtlbasedir : :class:`str`, optional, defaults to ``None``
        Formattable path to a directory that hosts alternate MTL ledgers  
        If ``None``, then look up the MTL directory from the
        $ALT_MTL_DIR environment variable. This will fail since that variable
        is not currently set in the desi code setup.
    ndirs : :class:`int`, optional, defaults to ``3``
        Number of alternate MTLs to process within altmtlbasedir
    numobs_from_ledger : :class:`bool`, optional, defaults to ``True``
        If ``True`` then inherit the number of observations so far from
        the ledger rather than expecting it to have a reasonable value
        in the `zcat.`
    secondary : :class:`bool`, optional, defaults to ``False``
        If ``True`` then process secondary targets instead of primaries
        for passed `survey` and `obscon`.
    quickRestart : :class:`bool`, optional, defaults to ``False``
        If ``True`` then copy original alternate MTLs from 
        altmtlbasedir/Univ*/survey/obscon/orig and  
    redoFA : :class:`bool`, optional, defaults to ``False``
        If ``True`` then automatically redo fiberassignment regardless of
        existence of fiberassign file in alternate fiberassign directory
    multiproc : :class:`bool`, optional, defaults to ``False``
        If ``True`` then run a single MTL update in a directory specified by
        nproc.
    nproc : :class:`int`, optional, defaults to None
        If multiproc is ``True`` this must be specified. Integer determines 
        directory of alternate MTLs to update.

    Returns
    -------
    :class:`str`
        The directory containing the ledger that was updated.
    :class:`str`
        The name of the MTL tile file that was updated.
    :class:`str`
        The name of the ZTILE file that was used to link TILEIDs to
        observing conditions and to determine if tiles were "done".
    :class:`~numpy.array`
        Information for the tiles that were processed.

    Notes
    -----
    - Assumes all of the relevant ledgers have already been made by,
      e.g., :func:`~LSS.SV3.altmtltools.initializeAlternateMTLs()`.
    """


    if mock:
        if targets is None:
            raise ValueError('If processing mocks, you MUST specify a target file')
    if debug:
        log.info('getosubp value: {0}'.format(getosubp))
    if ('trunk' in altmtlbasedir.lower()) or  ('ops' in altmtlbasedir.lower()):
        raise ValueError("In order to prevent accidental overwriting of the real MTLs, please remove \'ops\' and \'trunk\' from your MTL output directory")
    assert((singleDate is None) or (type(singleDate) == bool))
    if multiproc:
        import multiprocessing as mp
        import logging

        logger=mp.log_to_stderr(logging.DEBUG)

    ### JL - Start of directory/loop variable construction ###


    # ADM first grab all of the relevant files.
    # ADM grab the MTL directory (in case we're relying on $MTL_DIR).
    ##mtldir = get_mtl_dir(mtldir)
    # ADM construct the full path to the mtl tile file.
    ##mtltilefn = os.path.join(mtldir, get_mtl_tile_file_name(secondary=secondary))
    # ADM construct the relevant sub-directory for this survey and
    # ADM set of observing conditions..
    form = get_mtl_ledger_format()
    resolve = True
    msg = "running on {} ledger with obscon={} and survey={}"
    if secondary:
        log.info(msg.format("SECONDARY", obscon, survey))
        resolve = None
    else:
        log.info(msg.format("PRIMARY", obscon, survey))
    
    
    
    if altmtlbasedir is None:
        log.critical('This will automatically find the alt mtl dir in the future but fails now. Bye.')
        assert(0)
    if debugOrig:
        iterloop = range(1)
    elif multiproc:
        iterloop = range(nproc, nproc+1)
    else:
        iterloop = range(ndirs)
    ### JL - End of directory/loop variable construction ###



    if quickRestart:
        raise NotImplementedError('There is no way the quick restart will work properly post refactor.')
        quickRestartFxn(ndirs = ndirs, altmtlbasedir = altmtlbasedir, survey = survey, obscon = obscon, multiproc = multiproc, nproc = nproc)

    ### JL - this loop is through all realizations serially or (usually) one realization parallelized
    for n in iterloop:
        if debugOrig:
            altmtldir = altmtlbasedir
        else:
            altmtldir = os.path.join(altmtlbasedir.format(mock_number=n), 'Univ000/')
            #altmtldir = altmtlbasedir + '/Univ{0:03d}/'.format(n)
        altmtltilefn = os.path.join(altmtldir, get_mtl_tile_file_name(secondary=secondary))
 
        althpdirname = desitarget.io.find_target_files(altmtldir, flavor="mtl", resolve=resolve,
                                     survey=survey, obscon=obscon, ender=form)
        
        altMTLTileTrackerFN = makeTileTrackerFN(altmtldir, survey = survey, obscon = obscon)
        altMTLTileTracker = Table.read(altMTLTileTrackerFN)
        #today = altMTLTileTracker.meta['Today']
        #endDate = altMTLTileTracker.meta['EndDate']

        actionList = altMTLTileTracker[np.invert(altMTLTileTracker['DONEFLAG'])]

        actionList.sort(['ACTIONTIME'])
        #if not (singletile is None):
        #   tiles = tiles[tiles['TILEID'] == singletile]
        
        #if testDoubleDate:
        #    raise NotImplementedError('this block needs to be moved for new organization of tiletracker.')
        #    log.info('Testing Rosette with Doubled Date only')
        #    cond1 = ((tiles['TILEID'] >= 298) & (tiles['TILEID'] <= 324))
        #    cond2 = ((tiles['TILEID'] >= 475) & (tiles['TILEID'] <= 477))
        #    log.info(tiles[tiles['TILEID' ] == 314])
        #    log.info(tiles[tiles['TILEID' ] == 315])
        #    tiles = tiles[cond1 | cond2 ]
        

        #for ots,famtlt,reprocFlag in datepairs:
        #while int(today) <= int(endDate):
        for action in actionList:

            if action['ACTIONTYPE'] == 'fa':

                OrigFAs, AltFAs, AltFAs2, TSs, fadates, tiles = do_fiberassignment(altmtldir, [action], survey = survey, obscon = obscon ,verbose = verbose, debug = debug, getosubp = getosubp, redoFA = redoFA, mock = mock, reproducing = reproducing)
                assert(len(OrigFAs))
                A2RMap, R2AMap = make_fibermaps(altmtldir, OrigFAs, AltFAs, AltFAs2, TSs, fadates, tiles, changeFiberOpt = changeFiberOpt, verbose = verbose, debug = debug, survey = survey , obscon = obscon, getosubp = getosubp, redoFA = redoFA )
            elif action['ACTIONTYPE'] == 'update':
                althpdirname, altmtltilefn, ztilefn, tiles = update_alt_ledger(altmtldir,althpdirname, altmtltilefn, action, survey = survey, obscon = obscon ,getosubp = getosubp, zcatdir = zcatdir, mock = mock, numobs_from_ledger = numobs_from_ledger, targets = targets, verbose = verbose, debug = debug)
            elif action['ACTIONTYPE'] == 'reproc':
                #returns timedict

                #raise NotImplementedError('make backup here before reprocessing. Then resume Debugging.')
                retval = reprocess_alt_ledger(altmtldir, action, obscon=obscon, survey = survey)
                if debug or verbose:
                    log.info(f'retval = {retval}')
                
            else:
                raise ValueError('actiontype must be `fa`, `update`, or `reproc`.')

            gc.collect()
            #retval = write_amtl_tile_tracker(altmtldir, None, None, today, obscon = obscon, survey = survey, mode = 'endofday')
            #log.info('write_amtl_tile_tracker retval = {0}'.format(retval))

            #today = nextDate(today)
            #log.info('----------')
            #log.info('----------')
            #log.info('----------')
            #log.info('moving to next day: {0}'.format(today))
            #log.info('----------')
            #log.info('----------')
            #log.info('----------')

            
        return althpdirname, altmtltilefn, altMTLTileTrackerFN, actionList

def plotMTLProb(mtlBaseDir, ndirs = 10, hplist = None, obscon = 'dark', survey = 'sv3', outFileName = None, outFileType = '.png', jupyter = False, debug = False, verbose = False):
    """Plots probability that targets were observed among {ndirs} alternate realizations
    of SV3. Uses default matplotlib colorbar to plot between 1-{ndirs} observations.

    Parameters
    ----------
    mtlBaseDir : :class:`str`
        The home directory of your alternate MTLs. Should not contain obscon
        or survey. String should be formattable (i.e. '/path/to/dirs/Univ{0:03d}')
    ndirs   : :class:`int`
        The number of alternate realizations to plot.
    survey : :class:`str`, optional, defaults to "sv3"
        Used to look up the correct ledger, in combination with `obscon`.
        Options are ``'main'`` and ``'svX``' (where X is 1, 2, 3 etc.)
        for the main survey and different iterations of SV, respectively.
    obscon : :class:`str`, optional, defaults to "dark"
        Used to look up the correct ledger, in combination with `survey`.
        Options are ``'dark'`` and ``'bright``' 
    hplist : :class:`arraylike`, optional, defaults to None
        List of healpixels to plot. If None, defaults to plotting all available
        healpixels
    outFileName : :class:`str`, optional, defaults to None
        If desired, save file to this location. This will 
        usually be desired, but was made optional for use in
        ipython notebooks.
    outFileType : :class:`str`, optional, defaults to '.png'
        If desired, save file with name "outFileName" with
        type/suffix outFileType. This will usually be desired, 
        but was made optional for use in ipython notebooks.


    
    Returns
    -------
    Nothing
    
    """
    ObsFlagList = np.array([])
    for i in range(ndirs):
        mtldir = mtlBaseDir.format(i) + '/' + survey + '/' + obscon
        MTL = np.sort(desitarget.io.read_mtl_in_hp(mtldir, 32, hplist, unique=True, isodate=None, returnfn=False, initial=False, leq=False), order = 'TARGETID')
        try:
            ObsFlagList = np.column_stack((ObsFlagList,MTL['NUMOBS'] > 0.5))
        except:
            log.info('This message should appear once, only for the first realization.')
            ObsFlagList = MTL['NUMOBS'] > 0.5
    if verbose or debug:
        log.info(ObsFlagList.shape)
    ObsArr = np.sum(ObsFlagList, axis = 1)


    #MTLList[i] = rfn.append_fields(MTLList[i], 'OBSFLAG', MTLList[i]['NUMOBS'] > 0, dtypes=np.dtype(bool))
    
    hist, bins = np.histogram(ObsArr, bins = np.arange(ndirs)+ 0.1)

    plt.figure()
    plt.plot(bins[1:]- 0.01, hist)
    plt.xlabel('Number of Realizations in which a target was observed')
    plt.ylabel('Number of targets')
    #plt.yscale('log')
    if len(hplist )> 100:
        plt.ylim(0, 8000)
    elif (len(hplist) > 4) & (obscon == 'dark'):
        plt.ylim(0, 1500)
    elif (len(hplist) > 4) & (obscon == 'bright'):
        plt.ylim(0, 500)
    if not (outFileName is None):
        plt.savefig(outFileName + '_vsNtarget' + outFileType)
    if not jupyter:
        plt.close()
    plt.figure()
    plt.scatter(MTL['RA'][ObsArr > 0], MTL['DEC'][ObsArr > 0], c = ObsArr[ObsArr > 0], s = 0.1)
    plt.xlabel('RA')
    plt.ylabel('DEC')
    cbar = plt.colorbar()
    cbar.set_label('Number of Realizations in which target was observed')
    if not (outFileName is None):
        plt.savefig(outFileName + '_vsRADEC' + outFileType )
    if not jupyter:
        plt.close()

#@profile
def makeBitweights(mtlBaseDir, ndirs = 64, hplist = None, obscon = 'dark', survey = 'sv3', debug = False, obsprob = False, splitByReal = False, verbose = False):
    """Takes a set of {ndirs} realizations of DESI/SV3 and converts their MTLs into bitweights
    and an optional PROBOBS, the probability that the target was observed over the realizations

    Parameters
    ----------
    mtlBaseDir : :class:`str`
        The home directory of your alternate MTLs. Should not contain obscon
        or survey. String should be formattable (i.e. '/path/to/dirs/Univ{0:03d}')
    ndirs   : :class:`int`
        The number of alternate realizations to process. 
    survey : :class:`str`, optional, defaults to "sv3"
        Used to look up the correct ledger, in combination with `obscon`.
        Options are ``'main'`` and ``'svX``' (where X is 1, 2, 3 etc.)
        for the main survey and different iterations of SV, respectively.
    obscon : :class:`str`, optional, defaults to "dark"
        Used to look up the correct ledger, in combination with `survey`.
        Options are ``'dark'`` and ``'bright``' 
    hplist : :class:`arraylike`, optional, defaults to None
        List of healpixels to plot. If None, defaults to plotting all available
        healpixels
    debug : :class:`bool`, optional, defaults to False
        If True, prints extra information showing input observation information
        and output bitweight information for the first few targets as well as
        the first few targets that were observed in at least one realization
    obsprob: class:`bool`, optional, defaults to False
        If True, returns TARGETID, BITWEIGHT, and OBSPROB. Else returns TARGETID 
        and BITWEIGHT only
    splitByReal: class:`bool`, optional, defaults to False
        If True, run for only a single realization but for all healpixels in hplist.

    Returns
    -------
    :class:`~numpy.array`
        Array of Target IDs 
    :class:`~numpy.array`
        Array of bitweights for those target ids
    :class:`~numpy.array`, optional if obsprob is True
        Array of probabilities a target gets observed over {ndirs} realizations
        
    """
    
    TIDs = None
    if splitByReal:

        from mpi4py import MPI
        if debug or verbose:
            log.info('mtlbasedir')
            log.info(mtlBaseDir)
            log.info(mtlBaseDir.format(0))
        ntar = desitarget.io.read_mtl_in_hp(mtlBaseDir.format(0) + '/' + survey + '/' + obscon, 32, hplist, unique=True, isodate=None, returnfn=False, initial=False, leq=False).shape[0]
        
        comm = MPI.COMM_WORLD
        mpi_procs = comm.size
        mpi_rank = comm.rank
        if debug or verbose:
            log.info('running on {0:d} cores'.format(mpi_procs))
        n_realization = ndirs
        realizations = np.arange(ndirs, dtype=np.int32)
        my_realizations = np.array_split(realizations, mpi_procs)[mpi_rank]
        MyObsFlagList = np.empty((my_realizations.shape[0], ntar), dtype = bool)
        #MTL = np.sort(desitarget.io.read_mtl_in_hp(mtldir, 32, hplist, unique=True, isodate=None, returnfn=False, initial=False, leq=False), order = 'TARGETID')
        for i, r in enumerate(my_realizations):
            mtldir = mtlBaseDir.format(i) + '/' + survey + '/' + obscon
            MTL = np.sort(desitarget.io.read_mtl_in_hp(mtldir, 32, hplist, unique=True, isodate=None, returnfn=False, initial=False, leq=False), order = 'TARGETID') 
            if TIDs is None:
                TIDs = MTL['TARGETID']
            else:
                assert(np.array_equal(TIDs, MTL['TARGETID']))

            MyObsFlagList[i][:] = MTL['NUMOBS'] > 0.5
        
        ObsFlagList = None
        bitweights = None
        obsprobs = None
        #gather_weights = None
        if mpi_rank == 0:
            #gather_weights = np.empty(len(bitweights), dtype=bool)
            ObsFlagList = np.empty ((ndirs, ntar), dtype = bool)
        comm.Gather(MyObsFlagList, ObsFlagList, root=0)
        if mpi_rank == 0:    
            if debug or verbose:
                print(ObsFlagList.shape)
            ObsArr = np.sum(ObsFlagList, axis = 0)
            obsprobs = ObsArr/ndirs
            if debug or verbose:
                print(np.min(ObsArr))
                print(np.max(ObsArr))
                print("ObsFlagList shape here: {0}".format(ObsFlagList.shape))
            bitweights = pack_bitweights(ObsFlagList.T)
            if debug or verbose:
                print('bitweights shape here: {0}'.format(bitweights.shape))
                print('TIDs shape here: {0}'.format(TIDs.shape))
            assert(not (TIDs is None))
        if obsprob:
            return TIDs, bitweights, obsprobs
        else:
            return TIDs, bitweights
            
    else:
        ObsFlagList = np.empty(ndirs)
        for i in range(ndirs):
            mtldir = mtlBaseDir.format(i) + '/' + survey + '/' + obscon
            MTL = np.sort(desitarget.io.read_mtl_in_hp(mtldir, 32, hplist, unique=True, isodate=None, returnfn=False, initial=False, leq=False), order = 'TARGETID')
            if TIDs is None:
                TIDs = MTL['TARGETID']
            else:
                assert(np.array_equal(TIDs, MTL['TARGETID']))
            try:
                ObsFlagList = np.column_stack((ObsFlagList,MTL['NUMOBS'] > 0.5))
            except:
                log.info('hplist[0] = {0:d}'.format(hplist[0]))
                log.info('This message should only appear once for the first realization.')
                ObsFlagList = MTL['NUMOBS'] > 0.5
        if debug or verbose:
            log.info(ObsFlagList.shape)
        ObsArr = np.sum(ObsFlagList, axis = 1)
        if debug or verbose:
            log.info(np.min(ObsArr))
            log.info(np.max(ObsArr))
        bitweights = pack_bitweights(ObsFlagList)

        assert(not (TIDs is None))
        if obsprob:
            
            obsprobs = ObsArr/ndirs

            return TIDs, bitweights, obsprobs
        else:
            return TIDs, bitweights





def writeBitweights(mtlBaseDir, ndirs = None, hplist = None, debug = False, outdir = None, obscon = "dark", survey = 'sv3', overwrite = False, allFiles = False, splitByReal = False, splitNChunks = None, verbose = False):
    """Takes a set of {ndirs} realizations of DESI/SV3 and converts their MTLs into bitweights
    and an optional PROBOBS, the probability that the target was observed over the realizations.
    Then writes them to (a) file(s)

    Parameters
    ----------
    mtlBaseDir : :class:`str`
        The home directory of your alternate MTLs. Should not contain obscon
        or survey. String should be formattable (i.e. '/path/to/dirs/Univ{0:03d}')
    ndirs   : :class:`int`
        The number of alternate realizations to process. 
    survey : :class:`str`, optional, defaults to "sv3"
        Used to look up the correct ledger, in combination with `obscon`.
        Options are ``'main'`` and ``'svX``' (where X is 1, 2, 3 etc.)
        for the main survey and different iterations of SV, respectively.
    obscon : :class:`str`, optional, defaults to "dark"
        Used to look up the correct ledger, in combination with `survey`.
        Options are ``'dark'`` and ``'bright``' 
    hplist : :class:`arraylike`, optional, defaults to None
        List of healpixels to plot. If None, defaults to plotting all available
        healpixels
    debug : :class:`bool`, optional, defaults to False
        If True, prints extra information showing input observation information
        and output bitweight information for the first few targets as well as
        the first few targets that were observed in at least one realization
    obsprob: class:`bool`, optional, defaults to False
        If True, returns TARGETID, BITWEIGHT, and OBSPROB. Else returns TARGETID 
        and BITWEIGHT only
    outdir : :class:`str`, optional, defaults to None
        The base directory in which to create the BitweightFiles output directory. 
        If None, defaults to one level above mtlBaseDir
    overwrite: class:`bool`, optional, defaults to False
        If True, will clobber already existing bitweight files.
    ***Fix this Option to autograb all tiles***
    allfiles: class:`bool`, optional, defaults to False
        If True, do not generate a bitweight file for each healpixel, but generate
        one "allTiles" file for the combination of healpixels
    splitByReal: class:`bool`, optional, defaults to False
        If True, run for only a single realization but for all healpixels in hplist
    

    Returns
    -------
    :class:`~numpy.array`
        Array of Target IDs 
    :class:`~numpy.array`
        Array of bitweights for those target ids
    :class:`~numpy.array`, optional if obsprob is True
        Array of probabilities a target gets observed over {ndirs} realizations
        
    """
    if outdir is None:
        log.info('No outdir provided')
        outdir = mtlBaseDir.split('/')[:-1]
        log.info('autogenerated outdir')
        log.info(outdir)
    if splitByReal:
        from mpi4py import MPI        
        comm = MPI.COMM_WORLD
        mpi_procs = comm.size
        mpi_rank = comm.rank
        if mpi_rank == 0:
            if not os.path.exists(outdir + '/BitweightFiles/' + survey + '/' + obscon):
                os.makedirs(outdir + '/BitweightFiles/' + survey + '/' + obscon)
    elif not os.path.exists(outdir + '/BitweightFiles/' + survey + '/' + obscon):
        os.makedirs(outdir + '/BitweightFiles/' + survey + '/' + obscon)
    if type(hplist) == int:
        hplist = [hplist]
    if allFiles:
        hpstring = 'AllTiles'
    else:
        hpstring = 'hp-'

        for hp in hplist:
            hpstring += str(hp)
    fn = outdir + '/BitweightFiles/' + survey + '/' + obscon + '/{0}bw-{1}-'.format(survey.lower(), obscon.lower()) + hpstring + '.fits'
    #    outdir + '/BitweightFiles/' + survey + '/' + obscon + '/{0}bw-{1}-'.format(survey.lower(), obscon.lower()) + hpstring + '.fits'
    if (not overwrite) and os.path.exists(fn):
        print('overwrite')
        print(overwrite)
        print('fn')
        print(fn)
        return None
    
    if not (splitNChunks is None):
        if debug or verbose:
            log.info('makeBitweights1')
            log.info("splitting into {0} chunks".format(splitNChunks))
        splits = np.array_split(hplist, int(splitNChunks))


        for i, split in enumerate(splits):
            if debug or verbose:
                log.info('split {0}'.format(i))
                log.info(split)
            if i == 0:
                TIDs, bitweights, obsprobs = makeBitweights(mtlBaseDir, ndirs = ndirs, hplist = split, debug = False, obsprob = True, obscon = obscon, survey = survey, splitByReal = splitByReal)
            else:
                TIDsTemp, bitweightsTemp, obsprobsTemp = makeBitweights(mtlBaseDir, ndirs = ndirs, hplist = split, debug = False, obsprob = True, obscon = obscon, survey = survey, splitByReal = splitByReal)
                
                if mpi_rank == 0:
                    if debug or verbose:
                        log.info('----')
                        log.info('mpi_rank: {0}'.format(mpi_rank))
                        log.info("TIDs shape: {0}".format(TIDs.shape))
                        log.info("bitweights shape: {0}".format(bitweights.shape))
                        log.info("obsprobs shape: {0}".format(obsprobs.shape))
                        log.info('----')
                        log.info('mpi_rank: {0}'.format(mpi_rank))
                        log.info("TIDsTemp shape: {0}".format(TIDsTemp.shape))
                        log.info("bitweightsTemp shape: {0}".format(bitweightsTemp.shape))
                        log.info("obsprobsTemp shape: {0}".format(obsprobsTemp.shape))
                    TIDs = np.hstack((TIDs, TIDsTemp))
                    bitweights = np.vstack((bitweights, bitweightsTemp))
                    obsprobs = np.hstack((obsprobs, obsprobsTemp))
    else:
        if debug or verbose:
            log.info('makeBitweights2')
        TIDs, bitweights, obsprobs = makeBitweights(mtlBaseDir, ndirs = ndirs, hplist = hplist, debug = False, obsprob = True, obscon = obscon, survey = survey, splitByReal = splitByReal)
    if splitByReal:
        if debug or verbose:
            log.info('----')
            log.info('mpi_rank: {0}'.format(mpi_rank))
        if mpi_rank == 0:
            if debug or verbose:
                log.info("TIDs shape: {0}".format(TIDs.shape))
                log.info("bitweights shape: {0}".format(bitweights.shape))
                log.info("obsprobs shape: {0}".format(obsprobs.shape))
            data = Table({'TARGETID': TIDs, 'BITWEIGHTS': bitweights, 'PROB_OBS': obsprobs},
                      names=['TARGETID', 'BITWEIGHTS', 'PROB_OBS'])
            
            data.write(fn, overwrite = overwrite)
    else:
        if debug or verbose:
            log.info("TIDs shape: {0}".format(TIDs.shape))
            log.info("bitweights shape: {0}".format(bitweights.shape))
            log.info("obsprobs shape: {0}".format(obsprobs.shape))
        data = Table({'TARGETID': TIDs, 'BITWEIGHTS': bitweights, 'PROB_OBS': obsprobs},
              names=['TARGETID', 'BITWEIGHTS', 'PROB_OBS'])
    
        data.write(fn, overwrite = overwrite)
    
def reprocess_alt_ledger(altmtldir, action, obscon="dark", survey = 'main', zcatdir = None):
    """
    Reprocess HEALPixel-split ledgers for targets with new redshifts.

    Parameters
    ----------
    hpdirname : :class:`str`
        Full path to a directory containing an MTL ledger that has been
        partitioned by HEALPixel (i.e. as made by `make_ledger`).
    zcat : :class:`~astropy.table.Table`, optional
        Redshift catalog table with columns ``TARGETID``, ``NUMOBS``,
        ``Z``, ``ZWARN``, ``ZTILEID``, and ``msaddcols`` at the top of
        the code for the Main Survey.
    obscon : :class:`str`, optional, defaults to "DARK"
        A string matching ONE obscondition in the desitarget bitmask yaml
        file (i.e. in `desitarget.targetmask.obsconditions`), e.g. "DARK"
        Governs how priorities are set using "obsconditions". Basically a
        check on whether the files in `hpdirname` are as expected.

    Returns
    -------
    :class:`dict`
        A dictionary where the keys are the integer TILEIDs and the values
        are the TIMESTAMP at which that tile was reprocessed.

    """
    tileid = action['TILEID']
    ts = str(tileid).zfill(6)
    FABaseDir = '/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'
    FAFN = FABaseDir + f'/{ts[0:3]}/fiberassign-{ts}.fits'

    fhtOrig = fitsio.read_header(FAFN)
    fadate = fhtOrig['RUNDATE']
    fanite = int(''.join(fadate.split('T')[0].split('-')))

    hpdirname = altmtldir + f'/{survey.lower()}/{obscon.lower()}/'

    fbadirbase = altmtldir + '/fa/' + survey.upper() +  '/' + str(fanite) + '/'

    #if getosubp:
    #    FAMapName = fbadirbase + '/orig/famap-' + ts + '.pickle'
    #else:
    FAMapName = fbadirbase + '/famap-' + ts + '.pickle'
    with open(FAMapName,'rb') as fl:
        (A2RMap, R2AMap) = pickle.load(fl,fix_imports = True)

    #zcat = make_zcat(zcatdir, dateTiles, obscon, survey)
    zcatdir = get_zcat_dir(zcatdir)
    zcat = make_zcat(zcatdir, [action], obscon, survey, allow_overlaps = True)
    log.info('ts = {0}'.format(ts))
    altZCat = makeAlternateZCat(zcat, R2AMap, A2RMap)



    #if getosubp:
    #    FAMapName = fbadirbase + '/orig/famap-' + ts + '.pickle'
    #else:
    #    FAMapName = fbadirbase + '/famap-' + ts + '.pickle'
    #with open(FAMapName,'rb') as fl:
    #    (A2RMapTemp, R2AMapTemp) = pickle.load(fl,fix_imports = True)
    t0 = time()
    log.info("Reprocessing based on altZCat with {} entries...t={:.1f}s"
             .format(len(altZCat), time()-t0))

    # ADM the output dictionary.
    timedict = {}

    # ADM bits that correspond to a "bad" observation in the zwarn_mask.
    Mxbad = "BAD_SPECQA|BAD_PETALQA|NODATA"

    # ADM find the general format for the ledger files in `hpdirname`.
    # ADM also returning the obsconditions.
    fileform, oc = desitarget.io.find_mtl_file_format_from_header(hpdirname, returnoc=True)
    # ADM also find the format for any associated override ledgers.
    overrideff = desitarget.io.find_mtl_file_format_from_header(hpdirname,
                                                     forceoverride=True)

    # ADM check the obscondition is as expected.
    if obscon != oc:
        msg = "File is type {} but requested behavior is {}".format(oc, obscon)
        log.critical(msg)
        raise RuntimeError(msg)

    # ADM check the altZCat has unique TARGETID/TILEID combinations.
    tiletarg = [str(tt["ZTILEID"]) + "-" + str(tt["TARGETID"]) for tt in zcat]
    if len(set(tiletarg)) != len(tiletarg):
        msg = "Passed altZCat does NOT have unique TARGETID/TILEID combinations!!!"
        log.critical(msg)
        raise RuntimeError(msg)

    # ADM record the set of tiles that are being reprocessed.
    reproctiles = set(altZCat["ZTILEID"])

    # ADM read ALL targets from the relevant ledgers.
    log.info("Reading (all instances of) targets for {} tiles...t={:.1f}s"
             .format(len(reproctiles), time()-t0))
    nside = desitarget.mtl._get_mtl_nside()
    theta, phi = np.radians(90-altZCat["DEC"]), np.radians(altZCat["RA"])
    pixnum = hp.ang2pix(nside, theta, phi, nest=True)
    pixnum = list(set(pixnum))
    targets = desitarget.io.read_mtl_in_hp(hpdirname, nside, pixnum, unique=False, tabform='ascii.ecsv')

    # ADM remove OVERRIDE entries, which should never need reprocessing.
    targets, _ = desitarget.mtl.remove_overrides(targets)

    # ADM sort by TIMESTAMP to ensure tiles are listed chronologically.
    targets = targets[np.argsort(targets["TIMESTAMP"])]

    # ADM for speed, we only need to work with targets with a altZCat entry.
    ntargs = len(targets)
    nuniq = len(set(targets["TARGETID"]))
    log.info("Read {} targets with {} unique TARGETIDs...t={:.1f}s"
             .format(ntargs, nuniq, time()-t0))
    log.info("Limiting targets to {} (unique) TARGETIDs in the altZCat...t={:.1f}s"
             .format(len(set(altZCat["TARGETID"])), time()-t0))
    s = set(altZCat["TARGETID"])
    ii = np.array([tid in s for tid in targets["TARGETID"]])
    targets = targets[ii]
    nuniq = len(set(targets["TARGETID"]))
    log.info("Retained {}/{} targets with {} unique TARGETIDs...t={:.1f}s"
             .format(len(targets), ntargs, nuniq, time()-t0))

    # ADM split off the updated target states from the unobserved states.
    _, ii = np.unique(targets["TARGETID"], return_index=True)
    unobs = targets[sorted(ii)]
    # ADM this should remove both original UNOBS states and any resets
    # ADM to UNOBS due to reprocessing data that turned out to be bad.
    targets = targets[targets["ZTILEID"] != -1]
    # ADM every target should have been unobserved at some point.
    if len(set(targets["TARGETID"]) - set(unobs["TARGETID"])) != 0:
        msg = "Some targets don't have a corresponding UNOBS state!!!"
        log.critical(msg)
        raise RuntimeError(msg)
    # ADM each target should have only one UNOBS state.
    if len(set(unobs["TARGETID"])) != len(unobs["TARGETID"]):
        msg = "Passed ledgers have multiple UNOBS states!!!"
        log.critical(msg)
        raise RuntimeError(msg)

    log.info("{} ({}) targets are in the unobserved (observed) state...t={:.1f}s"
             .format(len(unobs), len(targets), time()-t0))

    # ADM store first-time-through tile order to reproduce processing.
    # ADM ONLY WORKS because we sorted by TIMESTAMP, above!
    _, ii = np.unique(targets["ZTILEID"], return_index=True)
    # ADM remember to sort ii so that the first tiles appear first.
    orderedtiles = targets["ZTILEID"][sorted(ii)]

    # ADM assemble a altZCat for all previous and reprocessed observations.
    altZCatfromtargs = np.zeros(len(targets), dtype=zcat.dtype)
    for col in altZCat.dtype.names:
        altZCatfromtargs[col] = targets[col]
    # ADM note that we'll retain the TIMESTAMPed order of the old ledger
    # ADM entries and new redshifts will (deliberately) be listed last.
    allaltZCat = np.concatenate([altZCatfromtargs, altZCat])
    log.info("Assembled a altZCat of {} total observations...t={:.1f}s"
             .format(len(allaltZCat), time()-t0))

    # ADM determine the FINAL observation for each TILED-TARGETID combo.
    # ADM must flip first as np.unique finds the FIRST unique entries.
    allaltZCat = np.flip(allaltZCat)
    # ADM create a unique hash of TILEID and TARGETID.
    tiletarg = [str(tt["ZTILEID"]) + "-" + str(tt["TARGETID"]) for tt in allaltZCat]
    # ADM find the final unique combination of TILEID and TARGETID.
    _, ii = np.unique(tiletarg, return_index=True)
    # ADM make sure to retain exact reverse-ordering.
    ii = sorted(ii)
    # ADM condition on indexes-of-uniqueness and flip back.
    allaltZCat = np.flip(allaltZCat[ii])
    log.info("Found {} final TARGETID/TILEID combinations...t={:.1f}s"
             .format(len(allaltZCat), time()-t0))

    # ADM mock up a dictionary of timestamps in advance. This is faster
    # ADM as no delays need to be built into the code.
    now = get_utc_date(survey="main")
    timestamps = {t: desitarget.mtl.add_to_iso_date(now, s) for s, t in enumerate(orderedtiles)}

    # ADM make_mtl() expects altZCats to be in Table form.
    allaltZCat = Table(allaltZCat)
    # ADM a merged target list to track and record the final states.
    mtl = Table(unobs)
    # ADM to hold the final list of updates per-tile.
    donemtl = []

    # ADM loop through the tiles in order and update the MTL state.
    for tileid in orderedtiles:
        # ADM the timestamp for this tile.
        timestamp = timestamps[tileid]

        # ADM restrict to the observations on this tile.
        altZCatmini = allaltZCat[allaltZCat["ZTILEID"] == tileid]
        # ADM check there are only unique TARGETIDs on each tile!
        if len(set(altZCatmini["TARGETID"])) != len(altZCatmini):
            msg = "There are duplicate TARGETIDs on tile {}".format(tileid)
            log.critical(msg)
            raise RuntimeError(msg)

        # ADM update NUMOBS in the altZCat using previous MTL totals.
        mii, zii = desitarget.geomask.match(mtl["TARGETID"], altZCatmini["TARGETID"])
        altZCatmini["NUMOBS"][zii] = mtl["NUMOBS"][mii] + 1

        # ADM restrict to just objects in the altZCat that match an UNOBS
        # ADM target (i,e that match something in the MTL).
        log.info("Processing {}/{} observations from altZCat on tile {}...t={:.1f}s"
                 .format(len(zii), len(altZCatmini), tileid, time()-t0))
        log.info("(i.e. removed secondaries-if-running-primaries or vice versa)")
        altZCatmini = altZCatmini[zii]

        # ADM ------
        # ADM NOTE: We could use trimtozcat=False without matching, and
        # ADM just continually update the overall mtl list. But, make_mtl
        # ADM doesn't track NUMOBS just NUMOBS_MORE, so we need to add
        # ADM complexity somewhere, hence trimtozcat=True/matching-back.
        # ADM ------
        # ADM push the observations on this tile through MTL.
        zmtl = desitarget.mtl.make_mtl(mtl, oc, zcat=altZCatmini, trimtozcat=True, trimcols=True)

        # ADM match back to overall merged target list to update states.
        mii, zii = desitarget.geomask.match(mtl["TARGETID"], zmtl["TARGETID"])
        # ADM update the overall merged target list.
        for col in mtl.dtype.names:
            if col.upper() == 'RA':
                continue
            elif col.upper() == 'DEC':
                continue
            mtl[col][mii] = zmtl[col][zii]
        # ADM also update the TIMESTAMP for changes on this tile.
        mtl["TIMESTAMP"][mii] = timestamp

        # ADM trimtozcat=True discards BAD observations. Retain these.
        tidmiss = list(set(altZCatmini["TARGETID"]) - set(zmtl["TARGETID"]))
        tii = desitarget.geomask.match_to(altZCatmini["TARGETID"], tidmiss)
        zbadmiss = altZCatmini[tii]
        # ADM check all of the missing observations are, indeed, bad.
        if np.any(zbadmiss["ZWARN"] & zwarn_mask.mask(Mxbad) == 0):
            msg = "Some objects skipped by make_mtl() on tile {} are not BAD!!!"
            msg = msg.format(tileid)
            log.critical(msg)
            raise RuntimeError(msg)
        log.info("Adding back {} bad observations from altZCat...t={:.1f}s"
                 .format(len(zbadmiss), time()-t0))

        # ADM update redshift information in MTL for bad observations.
        mii, zii = desitarget.geomask.match(mtl["TARGETID"], zbadmiss["TARGETID"])
        # ADM update the overall merged target list.
        # ADM Never update NUMOBS or NUMOBS_MORE using bad observations.
        for col in set(zbadmiss.dtype.names) - set(["NUMOBS", "NUMOBS_MORE", "RA", "DEC"]):
            if col.upper() == 'RA':
                continue
            elif col.upper() == 'DEC':
                continue
            mtl[col][mii] = zbadmiss[col][zii]
        # ADM also update the TIMESTAMP for changes on this tile.
        mtl["TIMESTAMP"][mii] = timestamp

        # ADM record the information to add to the output ledgers...
        donemtl.append(mtl[mtl["ZTILEID"] == tileid])

        # ADM if this tile was actually reprocessed (rather than being a
        # ADM later overlapping tile) record the TIMESTAMP...
        if tileid in reproctiles:
            timedict[tileid] = timestamp

    # ADM collect the results.
    mtl = Table(np.concatenate(donemtl))

    # ADM re-collect everything on pixels for writing to ledgers.
    nside = desitarget.mtl._get_mtl_nside()
    theta, phi = np.radians(90-mtl["DEC"]), np.radians(mtl["RA"])
    pixnum = hp.ang2pix(nside, theta, phi, nest=True)

    # ADM loop through the pixels and update the ledger, depending
    # ADM on whether we're working with .fits or .ecsv files.
    ender = get_mtl_ledger_format()
    for pix in set(pixnum):
        # ADM grab the targets in the pixel.
        ii = pixnum == pix
        mtlpix = mtl[ii]

        # ADM the correct filenames for this pixel number.
        fn = fileform.format(pix)
        overfn = overrideff.format(pix)

        # ADM if an override ledger exists, update it and recover its
        # ADM relevant MTL entries.
        if os.path.exists(overfn):
            overmtl = process_overrides(overfn)
            # ADM add any override entries TO THE END OF THE LEDGER.
            mtlpix = vstack([mtlpix, overmtl])

        # ADM if we're working with .ecsv, simply append to the ledger.
        if ender == 'ecsv':
            f = open(fn, "a")
            astropy.io.ascii.write(mtlpix, f, format='no_header', formats=mtlformatdict)
            f.close()
        # ADM otherwise, for FITS, we'll have to read in the whole file.
        else:
            ledger, hd = fitsio.read(fn, extname="MTL", header=True)
            done = np.concatenate([ledger, mtlpix.as_array()])
            fitsio.write(fn+'.tmp', done, extname='MTL', header=hd, clobber=True)
            os.rename(fn+'.tmp', fn)
    retval = write_amtl_tile_tracker(altmtldir, [action], obscon = obscon, survey = survey)
    return timedict    

 
def write_amtl_tile_tracker(dirname, tiles, obscon = 'dark', survey = 'main'):
    """Write AMTL Processing times into TileTrackers

    Parameters
    ----------
    dirname : :class:`str`
        The path to the AMTL directory.
    tiles : :class`astropy.table or numpy.recarray`
        The tiles which were processed in this AMTL loop iteration
    timestamp : :class:`str`
        the time at which the AMTL updates were performed
    obscon : :class:`str`
        The observing conditions of the tiles that were processed. "dark" or "bright"
    survey : :class:`str`
        The survey of the tiles that were processed. "main" or "sv3"

    Returns
    -------
    :class:`int`
        The number of targets that were written to file.
    :class:`str`
        The name of the file to which targets were written.
    """
    #if len(tiles) == 1:
    #    tiles = [tiles]
    TileTrackerFN =  makeTileTrackerFN(dirname, survey, obscon)
    log.info(TileTrackerFN)
    if os.path.isfile(TileTrackerFN):
        TileTracker = Table.read(TileTrackerFN, format = 'ascii.ecsv')

    #if mode.lower() == 'update':
    #    dateKey = 'ALTARCHIVEDATE'
    #elif mode.lower() == 'fa':
    #    dateKey = 'ALTFADATE'
    #elif mode.lower() == 'endofday':
    #    TileTracker.meta['Today'] = today
    #    TileTracker.write(TileTrackerFN, format = 'ascii.ecsv', overwrite = True)
    #    return 'only wrote today in metadata'
    for t in tiles:
        log.info('t = {0}'.format(t))
        tileid = t['TILEID']
        #reprocFlag = t['REPROCFLAG']
        actionType = t['ACTIONTYPE']
        cond = (TileTracker['TILEID'] == tileid) & (TileTracker['ACTIONTYPE'] == actionType)
        log.info('for tile {0}, number of matching tiles = {1}'.format(tileid, np.sum(cond)))
        #debugTrap = np.copy(TileTracker[dateKey])
        TileTracker['DONEFLAG'][cond] = True
    
    assert(not (np.all(np.invert(TileTracker['DONEFLAG']))))

    #if mode == 'update':
    #    todaysTiles = TileTracker[TileTracker['ORIGMTLDATE'] == today]
    #    #if np.sum(todaysTiles['ALTARCHIVEDATE'] == None) == 0:
            
    TileTracker.write(TileTrackerFN, format = 'ascii.ecsv', overwrite = True)
    return 'done'
