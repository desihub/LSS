'''
Calculate NTILE and create a joined table
'''

import numpy as np
import os
from astropy.table import Table, join, unique
import argparse
import LSS.main.cattools as ct
import LSS.common_tools as common
from LSS.globals import main
from LSS.main.cattools import count_tiles_better
import fitsio

parser = argparse.ArgumentParser()
parser.add_argument("--indir")
parser.add_argument("--outdir", default = None)
parser.add_argument("--tracer")
parser.add_argument("--tileloc_add", default = 'n', help='set to y if your input catalog does not have TILELOCID column yet')
parser.add_argument("--skip_spec", default = 'n', help='set to y if spectrograph information already accounted for; combdark_wdupspec_zdone.fits already exists')
args = parser.parse_args()

if args.outdir == None:
    out_dir = args.indir
else:
    out_dir = args.outdir

if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    print('made ' + out_dir)


intablefn = args.indir + "pota-DARK" +  ".fits"
intable = Table.read(intablefn)
intable = intable[intable["COLLISION"]==False]



specfo = "/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/datcomb_dark_spec_zdone.fits"
print('loading specf file '+specfo)
specf = Table(fitsio.read(specfo))

#mainp = main(type,'daily')
mainp = main(type,'iron','Y1')

mt = mainp.mtld
tiles = mainp.tiles
imbits = mainp.imbits #mask bits applied to targeting
ebits = mainp.ebits #extra mask bits we think should be applied


tsnrcut = mainp.tsnrcut
dchi2 = mainp.dchi2
tnsrcol = mainp.tsnrcol        
zmin = mainp.zmin
zmax = mainp.zmax
badfib = mainp.badfib

wd = mt['SURVEY'] == 'main'
wd &= mt['ZDONE'] == 'true'
wd &= mt['FAPRGRM'] == 'dark'
wd &=mt['ZDATE'] < 20220900

mtld = mt[wd]
sel = np.isin(specf['TILEID'],mtld['TILEID'])
specf = specf[sel]
specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']           
print('loaded specf file '+specfo)
specfc = common.cut_specdat(specf,badfib=mainp.badfib)
gtl = np.unique(specfc['TILELOCID'])
# #common.write_LSS(gtl, out_dir + "/gtl.fits")

if args.skip_spec != 'y':
    print('combining with spec info')
    mask_coll = True
    #kc = ['LOCATION','FIBER','TILEID','TILELOCID','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG','PRIORITY']
    kc = ['LOCATION','FIBER','TILEID','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG']

    ct.combran_wdupspec(rann=0,tp='dark',lspecdir='',specf=specf,infile =intablefn,keepcols=kc,mask_coll=mask_coll, alt_out = out_dir, mock_priority_mask = 'n')
else: 
    print('skipping combran_wdupspec')


# COUNTING

intablefn = out_dir + '/combdark_wdupspec_zdone.fits' 
intable = Table.read(intablefn)


if args.tileloc_add == 'y':
    intable['TILELOCID'] = 10000*intable['TILEID'] + intable['LOCATION']
    intable.write(out_dir + '/combdark_wdupspec_zdone.fits' , overwrite = True)

intable = intable[intable["COLLISION"]==False]

if args.tracer == "ELG":
    intable = intable[intable["DESI_TARGET"]==34]
elif args.tracer == "LRG":
    intable = intable[intable["DESI_TARGET"]==1]
elif args.tracer == "QSO":
    intable = intable[intable["DESI_TARGET"]==4]

    
intable.write(out_dir + '/combdark_wdupspec_zdone_' + args.tracer + '.fits', overwrite = True)
    
print('counting tiles')

tc = ct.count_tiles_better('mock',pd = args.tracer, gtl=gtl, indir = out_dir)
print("going to write")
common.write_LSS(tc, out_dir + "/pota-DARK_withntile_" + args.tracer +".fits")
print("done counting")

## CUT INTABLE TO GTL AS WELL!!
wg = np.isin(intable['TILELOCID'],gtl)
intable = intable[wg]
#intable.write(out_dir + '/goodtilelocid_' + args.tracer + '.fits', overwrite = True)


print("joining full")
fulljj = join(intable, tc, keys = "TARGETID")
fulljj.write(out_dir + '/pota-DARK_joined_unreduced_%s.fits'%args.tracer)

print("joining unique")
reduced = intable[np.unique(intable["TARGETID"], return_index=True)[1]]
jj = join(reduced, tc, keys = "TARGETID")
#jj.write(args.indir + 'jointest.fits', overwrite = True)


#jju = unique(jj, keys = "TARGETID")

#jj = jj['TARGETID','RA','DEC', 'RSDZ', 'TRUEZ', 'GALCAP', 'NZ', 'RAW_NZ', 'COLLISION','NTILE','TILES','TILELOCIDS']
common.write_LSS(jj, out_dir + "/pota-DARK_joined_%s.fits"%args.tracer)
print("done")