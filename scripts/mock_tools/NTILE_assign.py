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
from datetime import datetime
startTime = datetime.now()

parser = argparse.ArgumentParser()
parser.add_argument("--indir")
parser.add_argument("--outdir", default = None)
parser.add_argument("--tracer")
parser.add_argument("--tileloc_add", default = 'n', help='set to y if your input catalog does not have TILELOCID column yet')
parser.add_argument("--skip_spec", default = 'n', help='set to y if spectrograph information already accounted for; comb{program}_wdupspec_zdone.fits already exists')
parser.add_argument("--prog", default = 'DARK')
args = parser.parse_args()

if args.outdir == None:
    out_dir = args.indir
else:
    out_dir = args.outdir

if not os.path.exists(out_dir):
    os.makedirs(out_dir)
    print('made ' + out_dir)

if args.prog == 'DARK':
    prog_small = 'dark'
elif args.prog == 'BRIGHT':
    prog_small = 'bright'
print("Doing " + args.prog + " program")

intablefn = args.indir + "pota-" + args.prog +  ".fits"
intable = Table.read(intablefn)
intable = intable[intable["COLLISION"]==False]



specfo = "/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/datcomb_" + prog_small + "_spec_zdone.fits"
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
wd &= mt['FAPRGRM'] == prog_small
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
    t0 = datetime.now()
    print('combining with spec info')
    mask_coll = True
    #kc = ['LOCATION','FIBER','TILEID','TILELOCID','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG','PRIORITY']
    kc = ['LOCATION','FIBER','TILEID','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG']

    ct.combran_wdupspec(rann=0,tp=prog_small,lspecdir='',specf=specf,infile =intablefn,keepcols=kc,mask_coll=mask_coll, alt_out = out_dir, mock_priority_mask = 'n')
    t1 = datetime.now()
    print("combran_wdupspec time:")
    print(t1 - t0)
else: 
    print('skipping combran_wdupspec')


# COUNTING

intablefn = out_dir + '/comb'+ prog_small +'_wdupspec_zdone.fits' 
intable = Table.read(intablefn)


if args.tileloc_add == 'y':
    intable['TILELOCID'] = 10000*intable['TILEID'] + intable['LOCATION']
    intable.write(out_dir + '/comb' + prog_small + '_wdupspec_zdone.fits' , overwrite = True)

intable = intable[intable["COLLISION"]==False]

if args.tracer == "ELG":
    intable = intable[intable["DESI_TARGET"]==34]
elif args.tracer == "LRG":
    intable = intable[intable["DESI_TARGET"]==1]
elif args.tracer == "QSO":
    intable = intable[intable["DESI_TARGET"]==4]
elif args.tracer == "BGS":
    intable = intable[intable["DESI_TARGET"]==2**60]

    
intable.write(out_dir + '/comb' + prog_small + '_wdupspec_zdone_' + args.tracer + '.fits', overwrite = True)
    
print('counting tiles')

t2 = datetime.now()
tc = ct.count_tiles_better('mock',pd = args.tracer, gtl=gtl, indir = out_dir, prog_ = 'bright')
t3 = datetime.now()
print("count_tiles_better time:")
print(t3 - t2)
print("going to write")
common.write_LSS(tc, out_dir + "/pota-" + args.prog + "_withntile_" + args.tracer +".fits")
print("done counting")

## CUT INTABLE TO GTL AS WELL!!
wg = np.isin(intable['TILELOCID'],gtl)
intable = intable[wg]
#intable.write(out_dir + '/goodtilelocid_' + args.tracer + '.fits', overwrite = True)


print("joining full")
t4 = datetime.now()
fulljj = join(intable, tc, keys = "TARGETID")
t5 = datetime.now()
print("unreduced joining time:")
print(t5 - t4)
fulljj.write(out_dir + '/pota-' + args.prog + '_joined_unreduced_%s.fits'%args.tracer)
t6 = datetime.now()
print("unreduced writing time:")
print(t6 - t5)

print("joining unique")
t7 = datetime.now()
reduced = intable[np.unique(intable["TARGETID"], return_index=True)[1]]
t8 = datetime.now()
jj = join(reduced, tc, keys = "TARGETID")
t9 = datetime.now()
print("reduced joining time:")
print(t9 - t8)
#jj.write(args.indir + 'jointest.fits', overwrite = True)


#jju = unique(jj, keys = "TARGETID")

#jj = jj['TARGETID','RA','DEC', 'RSDZ', 'TRUEZ', 'GALCAP', 'NZ', 'RAW_NZ', 'COLLISION','NTILE','TILES','TILELOCIDS']
common.write_LSS(jj, out_dir + "/pota-"+ args.prog + "_joined_%s.fits"%args.tracer)
t10 = datetime.now()
print("reduced writing time:")
print(t10 - t9)
print("done")
print("script time:")
print(datetime.now() - startTime)