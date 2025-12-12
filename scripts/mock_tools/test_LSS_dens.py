import sys,os
import numpy as np
import fitsio
import LSS.common_tools as common


def comp_numbers(tracer,mockdir,dataver='loa-v1/LSScats/v2/',survey='DA2',rootdir='/global/cfs/cdirs/desi/survey/catalogs/'):
    datadir = rootdir+survey+'/LSS/'+dataver+'/'
    dat = fitsio.read(datadir+tracer+'_full_noveto.dat.fits',columns=['ZWARN'])
    mock = fitsio.read(mockdir+tracer+'_full_noveto.dat.fits',columns=['ZWARN'])
    ndat = len(dat[dat['ZWARN']!=999999])
    nmock = len(mock[mock['ZWARN']!=999999])
    print(tracer+' numbers in full_noveto, ndat, nmock, ratio')
    print(ndat,nmock,ndat/nmock)

    dat = fitsio.read(datadir+tracer+'_full_HPmapcut.dat.fits',columns=['ZWARN'])
    mock = fitsio.read(mockdir+tracer+'_full_HPmapcut.dat.fits',columns=['ZWARN'])
    ndat = len(dat[dat['ZWARN']!=999999])
    nmock = len(mock[mock['ZWARN']!=999999])
    print(tracer+' number in area in full_HPmapcut, ndat, nmock, ratio')
    print(len(dat),len(mock),len(dat)/len(mock))

    print(tracer+' number assigned in full_HPmapcut, ndat, nmock, ratio')
    print(ndat,nmock,ndat/nmock)
    dat = fitsio.read(datadir+'nonKP/'+tracer+'_clustering.dat.fits',columns=['Z'])
    mock = common.read_hdf5_blosc(mockdir+tracer+'_clustering.dat.h5',columns=['Z'])
    if tracer == 'QSO':
        zmin = 0.8
        zmax = 3.5
    if tracer == 'LRG':
        zmin = 0.4
        zmax = 1.1
    if tracer[:3] == 'ELG':
        zmin = 0.8
        zmax = 1.6
    seld = dat['Z'] > zmin
    seld &= dat['Z'] < zmax
    #print(len(dat))
    dat = dat[seld]
    #print(len(dat))
    selm = mock['Z'] > zmin
    selm &= mock['Z'] < zmax
    mock = mock[selm]
    
    print('numbers in clustering, ndat, nmock, ratio')
    print(len(dat),len(mock),len(dat)/len(mock))



mockdir = sys.argv[1]

comp_numbers('QSO',mockdir)
comp_numbers('ELG_LOPnotqso',mockdir)
comp_numbers('LRG',mockdir)
