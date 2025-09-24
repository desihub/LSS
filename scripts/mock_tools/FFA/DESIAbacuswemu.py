from matplotlib import pyplot as plt
import matplotlib.cm as cm
from astropy.table import Table
from astropy.time import Time
from astropy.io import (fits, ascii)
import numpy as np
import os
import healpy as hp
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--mocknumber')
parser.add_argument('--tracer')
parser.add_argument('--basedir')
parser.add_argument('--outdir', default = "/global/cfs/cdirs/desi/survey/catalogs/main/mocks/FAemu_preliminary/sikandar/Updated_Code_CFC/fof_v1.0/in/")
parser.add_argument('--overwrite', default = 'n', help = 'do you want to overwrite the final output file?')
args = parser.parse_args()


mock = args.mocknumber
galcap = 'B'
target = args.tracer
print('mock =', mock)
print('galcap =', galcap)
print('target =', target)


RA = []
dec = []
targetid = []
truez = []
rsdz = []
gcap = []
nz = []
raw_nz = []
ntile = []

basedir = args.basedir
namin_gal = basedir + "/mock%s/pota-DARK_joined_%s.fits"%(mock, target)
hdul = fits.open(namin_gal)
hdul.info()
hdul[1].header
data_gal = hdul[1].data



RA = np.append(RA, data_gal['RA'])
dec = np.append(dec, data_gal['DEC'])
targetid = np.append(targetid, data_gal['TARGETID'])
truez = np.append(truez,data_gal['TRUEZ'])
rsdz = np.append(rsdz, data_gal['RSDZ'])
#gcap = np.append(gcap, data_gal['GALCAP'])
#nz = np.append(nz, data_gal['NZ'])
#raw_nz = np.append(raw_nz, data_gal['RAW_NZ'])
ntile = np.append(ntile, data_gal['NTILE'])

slope_gcap0, ra_gcap0 = -2.3, 97
slope_gcap1, ra_gcap1 = 2.3, 287

gcap = np.array(['S' for i in range(0, len(RA))])

gcap[np.where((dec > slope_gcap0 * (RA - ra_gcap0))
            & (dec > slope_gcap1 * (RA - ra_gcap1)))] = 'N'



if target == 'ELG':
    target0 = 34
if target == 'LRG':
    target0 = 1
if target == 'QSO':
    target0 = 4
print(target)
print(target0)
#id_tarcut = np.where(data_ref["DESI_target"] == target0)

DESI_target = target0 * np.ones(len(RA), dtype=int)
rchmsk = np.ones(len(RA), dtype=bool)                                            
selmsk = np.ones(len(RA), dtype=bool) 
#ntile = -np.ones(len(RA), dtype=int) 

# print(DESI_target)
# print(rchmsk)
# print(selmsk)
# print(ntile)

outdir = args.outdir

if not os.path.exists(outdir):
    os.makedirs(outdir)
    print('made ' + outdir)

namout_txt = outdir + target +'_Abacus2mock_forFAemu_m'+ mock +'.dat'

#RA_cut = [20.0, 40.0]
#dec_cut = [-10.0, 10.0]
RA_cut = [-10.0, 370.0]
dec_cut = [-100.0, 100.0]

id_out = np.where((RA > RA_cut[0]) & (RA < RA_cut[1])
                  & (dec > dec_cut[0]) & (dec < dec_cut[1])
                 )

print(np.shape(id_out))

if args.overwrite == 'y':
    overwriteflag = True
else:
    overwriteflag = False

ascii.write([RA[id_out],
             dec[id_out],
             truez[id_out],
             rsdz[id_out],
             DESI_target[id_out],
             rchmsk[id_out],
             selmsk[id_out],
             ntile[id_out],
             gcap[id_out],
             targetid[id_out]
            ],
            namout_txt,
            format='no_header',
            overwrite=overwriteflag,
            formats={'col0': '%8.5f',
                     'col1': '%8.5f',
                     'col2': '%8.5f',
                     'col3': '%8.5f',
                     'col7': '%8.5f',
                     'col9': '%12i'
                    }
           )

print('wrote file:', namout_txt)