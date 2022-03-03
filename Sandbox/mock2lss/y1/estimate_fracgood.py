import glob, os, sys
import fitsio as F
import numpy as np

path='/global/cscratch1/sd/acarnero/codes/LSS/Sandbox/mock2lss/fiberassigment'
fba_files = glob.glob(path+"/mocks*.fits")
allavail_tar=np.concatenate([F.read(f,ext="FAVAIL") for f in fba_files])
allassign_tar=np.concatenate([F.read(f,ext="FASSIGN",columns=["FIBER","TARGETID","LOCATION","FIBERSTATUS"]) for f in fba_files])
targ=F.read("/global/cscratch1/sd/acarnero/codes/LSS/Sandbox/mock2lss/mockTargets_000_FirstGen_CutSky_alltracers_sv3bits.fits")

mask_elg = np.where(targ['SV3_DESI_TARGET'] == 2)
mask_lrg = np.where(targ['SV3_DESI_TARGET'] == 1)
mask_qso = np.where(targ['SV3_DESI_TARGET'] == 4)
targ_elg = targ[mask_elg]
targ_lrg = targ[mask_lrg]
targ_qso = targ[mask_qso]

isas_elg=np.zeros(len(targ_elg),dtype=bool)
isas_lrg=np.zeros(len(targ_lrg),dtype=bool)
isas_qso=np.zeros(len(targ_qso),dtype=bool)

'''idas_elg=np.isin(targ_elg["TARGETID"], allassign_tar["TARGETID"])
idas_lrg=np.isin(targ_lrg["TARGETID"], allassign_tar["TARGETID"])
idas_qso=np.isin(targ_qso["TARGETID"], allassign_tar["TARGETID"])'''

idas_elg=np.isin(targ_elg["TARGETID"], allavail_tar["TARGETID"])
idas_lrg=np.isin(targ_lrg["TARGETID"], allavail_tar["TARGETID"])
idas_qso=np.isin(targ_qso["TARGETID"], allavail_tar["TARGETID"])

isas_elg[idas_elg] = True
isas_lrg[idas_lrg] = True
isas_qso[idas_qso] = True
assign_elg=targ_elg[idas_elg]
assign_lrg=targ_lrg[idas_lrg]
assign_qso=targ_qso[idas_qso]
print('Fraction ELG', len(assign_elg)/len(targ_elg))
print('Fraction LRG', len(assign_lrg)/len(targ_lrg))
print('Fraction QSO', len(assign_qso)/len(targ_qso))
