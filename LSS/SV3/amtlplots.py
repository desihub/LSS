import fitsio
from astropy.table import Table, join, vstack
import matplotlib.pyplot as plt
import numpy as np
import sys
import importlib.util
import scipy.stats as stats
import desitarget.io
import matplotlib as mpl
from desitarget.geomask import is_in_hp, nside2nside, pixarea2nside
from desimodel.footprint import is_point_in_desi, tiles2pix
from desitarget.mtl import add_to_iso_date


def tile2timestamp(MTLDir, TileID):
    doneTilesFile = Table.read(MTLDir + 'mtl-done-tiles.ecsv')
    cond = doneTilesFile['TILEID'] == TileID
    return doneTilesFile[cond]['TIMESTAMP']

def warpCoordsAroundCenter(RAs, Decs, RACent = None):
    if RACent is None:
        print('WARNING: Taking RA Center to be mean of all RA Coordinates.')
        RACent = np.mean(RAs)
    
    RACorr = RACent + (RACent - RAs)/np.cos(np.deg2rad(Decs))
    return RACorr, Decs

def TB2Colors(mtl):
    TB = mtl['DESI_TARGET']
    numobs = mtl['NUMOBS']
    ttypes = ['ELG', 'QSO', 'LRG']
    tbits = [2,4,1]
    tcolors = ['g', 'b', 'r']
    ColorArray = np.zeros(TB.shape, dtype = str)
    for ttype, tbit, tc in zip(ttypes, tbits, tcolors):
        for tb, no, ii in zip(TB, numobs, range(TB.shape[0])):
            isTtype = (tb & tbit) == tbit
            if isTtype:
                #print('assigning a color')
                np.put(ColorArray, ii, tc)
                
    return ColorArray


def completenessPlot2D(BaseDir, HPList, obscon = 'dark', survey = 'main', isodate = None, spacing = 5, CompQty = 'PROB_OBS',
                     RALimits = (-180, +180), nbinsRA = 360, nbinsDec = 180, DecLimits = (-90, +90), radius = 2, RACent = None, DecCent = None, useBitweightFile = False,
                    figsize = (12,8), rosetteNum = 'Not Specified', datadir = '/global/cfs/cdirs/desi/public/edr/vac/edr/lss/v2.0/LSScats/full/', dataFN = 'all', verbose = True):
    if not RACent is None:
        assert(RALimits[0] < RACent)
        assert(RALimits[1] > RACent)
        
    if not DecCent is None:
        assert(DecLimits[0] < DecCent)
        assert(DecLimits[1] > DecCent)
    #MTL = desitarget.io.read_mtl_in_hp(BaseDir + 'Univ000/{0}/{1}/'.format(survey, obscon), 32, GoodHPList, unique=True, isodate=isodate, returnfn=False, initial=False, leq=False)
    if (survey == 'sv3') and (BaseDir == '/global/cfs/cdirs/desi/public/edr/vac/edr/lss/v2.0/altmtl/'):
        if useBitweightFile:
            raise ValueError('Cannot use a separate bitweight file directly if using the EDR products.')
        TIDsInHP = desitarget.io.read_targets_in_hp(BaseDir + 'Univ000/{0}/'.format(obscon), 32, HPList, columns=['TARGETID'], mtl = True)
    elif useBitweightFile:
        TIDsInHP = desitarget.io.read_targets_in_hp(BaseDir + 'Univ000/{0}/{1}/'.format(survey, obscon), 32, HPList, columns=['TARGETID'], mtl = True)
        try:
            BWFile = Table.read(BaseDir + '/BitweightFiles/{0}/{1}/{0}bw-{1}-allTiles.fits'.format(survey, obscon))
        except:
            concatenateBWFiles(BaseDir, hplist, obscon = obscon, skipFailures=True, overwrite = True)
            BWFile = Table.read(BaseDir + '/BitweightFiles/{0}/{1}/{0}bw-{1}-allTiles.fits'.format(survey, obscon))
    else:
        TIDsInHP = desitarget.io.read_targets_in_hp(BaseDir + 'Univ000/{0}/{1}/'.format(survey, obscon), 32, HPList, columns=['TARGETID'], mtl = True)
        
    
    TIDsInHP = Table([TIDsInHP['TARGETID']], names = ['TARGETID'])
    if (obscon == 'dark') and (dataFN == 'all'):
        dataFile = vstack((Table.read(datadir + 'LRG_full.dat.fits'), Table.read(datadir + 'ELGnotqso_full.dat.fits'), Table.read(datadir + 'QSO_full.dat.fits')))
    elif (obscon == 'bright') and (dataFN == 'all'):
        dataFile = Table.read(datadir + 'BGS_ANY_full.dat.fits')
    elif type(dataFN) == type('a'):
        dataFile = Table.read(datadir + dataFN)
    else: 
        dataFile = None
        
    if not (dataFile is None):
        combo = join(TIDsInHP, dataFile, keys = ['TARGETID'])
    else:
        pass
    
    if useBitweightFile:
        combo = join(combo, BWFile, keys = ['TARGETID'])
        
    if (survey == 'sv3'):
        NRosettes, counts = np.unique(combo['ROSETTE_NUMBER'], return_counts = True)
        if len(NRosettes) == 1:
            ind = np.argmin(combo['ROSETTE_R'])
            rc = combo['RA'][ind]
            dc = combo['DEC'][ind]
            RALimits = ((rc - radius/np.cos(np.radians(dc))), (rc + radius/np.cos(np.radians(dc))))
            DecLimits = ((dc - radius), (dc + radius))
        elif (len(NRosettes) > 1) & (len(NRosettes) < 5):
            indRN = np.argmax(counts)
            combo = combo[combo['ROSETTE_NUMBER'] == NRosettes[indRN]]
            indRad = np.argmin(combo['ROSETTE_R'])
            rc = combo['RA'][indRad]
            dc = combo['DEC'][indRad]
            RALimits = ((rc - radius/np.cos(np.radians(dc))), (rc + radius/np.cos(np.radians(dc))))
            DecLimits = ((dc - radius), (dc + radius))
            
            
    
            
    if RALimits is None:
        binsRA = nbinsRA
    else:
        binsRA = np.linspace(RALimits[0], RALimits[-1], nbinsRA)
        
    if DecLimits is None:
        binsDec = nbinsDec
    else:
        binsDec = np.linspace(DecLimits[0], DecLimits[-1], nbinsDec)

    if CompQty == 'diff':
        hist2D, binsRA, binsDec, binnum = stats.binned_statistic_2d(combo['RA'], combo['DEC'], combo['PROB_OBS'] - combo['FRACZ_TILELOCID'], statistic = 'mean', bins = (binsRA, binsDec))
        #hist2DFZ, binsRA, binsDec, binnum = stats.binned_statistic_2d(combo['RA'], combo['DEC'], combo['FRACZ_TILELOCID'], statistic = 'mean', bins = (binsRA, binsDec))
        #hist2D = hist2DPO - hist2DFZ
    else:
        hist2D, binsRA, binsDec, binnum = stats.binned_statistic_2d(combo['RA'], combo['DEC'], combo[CompQty], statistic = 'mean', bins = (binsRA, binsDec))

    plt.figure(figsize = figsize)
    if dataFN == 'all':
        plt.title('rosetteNum = {0}; all {1} tracers'.format(rosetteNum, obscon))
    else:
        plt.title('rosetteNum = {0}; {1}'.format(rosetteNum, dataFN.split('.')[0]))

    
    
    
    if CompQty == 'diff':
        cmap = mpl.cm.RdBu
        norm = mpl.colors.Normalize(vmin=-0.2, vmax=+0.2)
    else:
        cmap = mpl.cm.viridis
        norm = mpl.colors.Normalize(vmin=0, vmax=1)
    
    mappable = mpl.cm.ScalarMappable(cmap = cmap, norm = norm)

    cb = plt.colorbar(mappable = mappable)
    if CompQty == 'diff':
        cb.set_label(label='PROB_OBS - FRACZ_TILELOCID',size = 20, labelpad = 10)
    else:
        cb.set_label(label='{0}'.format(CompQty),size = 20, labelpad = 10)



    plt.imshow(hist2D.T, origin = 'lower', cmap = cmap, norm = norm)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.xticks(ticks = np.arange(nbinsRA)[0::spacing], labels = np.around(binsRA[0::spacing], decimals = 2))
    plt.yticks(ticks = np.arange(nbinsDec)[0::spacing], labels = np.around(binsDec[0::spacing], decimals = 2))
    plt.tight_layout()
    if verbose:
        print('RALimits = {0}'.format(RALimits))
        print('DecLimits = {0}'.format(DecLimits))
        print('nanminRA, nanmaxRA, nanminDec, nanmaxDec')
        print(np.nanmin(combo['RA']))
        print(np.nanmax(combo['RA']))
        print(np.nanmin(combo['DEC']))
        print(np.nanmax(combo['DEC']))