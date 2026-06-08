import fitsio
from astropy.table import Table
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import sys,os

#set LSSCODE, e.g. via export LSSCODE=$HOME on the command line if that is where you cloned the repo 
#sys.path.append(os.environ['LSSCODE']+'/LSS/py')
from LSS import common_tools as common
from LSS.globals import main


import argparse

parser = argparse.ArgumentParser()
#parser.add_argument("--indir", help="base directory for catalogs", default='/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v1.5/')
#parser.add_argument("--outdir", help="output directory for the plots", default=None)
parser.add_argument("--survey", help="e.g., Y1, DA2, main", default='DA2')
parser.add_argument("--specver", help="e.g., iron, loa-v1, daily", default='loa-v1')
parser.add_argument("--version", help="e.g., test, v1.1", default='v1.1')
parser.add_argument("--tracers", help="which tracers to plot completeness for; use 'all' to make all tracer types", default=None)
parser.add_argument("--mk_ntile_dark", help="whether to make the dark time ntile plot", default='n')
parser.add_argument("--mk_ntile_bright", help="whether to make the bright time ntile plot", default='n')
parser.add_argument("--fontsize", help="base fontsize for labels", default=None)
args = parser.parse_args()


indir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/'+args.specver+'/LSScats/'+args.version+'/'
outdir = '/global/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/analysis/'+args.specver+'/LSScats/'+args.version+'/'
#if args.outdir is None:
    #outdir = args.indir.replace('dvs_ro','global')+'KPplots/'
#    outdir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/KPplots/'
#else:
#    outdir = args.outdir

if not os.path.exists(outdir):
    os.makedirs(outdir)
    print('made '+outdir)

if args.tracers == 'all':
    list_tracer = ['BGS_BRIGHT','ELG_LOPnotqso','LRG','QSO']
elif args.tracers is not None:
    list_tracer = [args.tracers]




nran = 0

for tr in list_tracer:
    mainp = main(tr,args.specver,survey=args.survey)
    if args.survey == 'Y1':
        data_full = fitsio.read(indir+'{}_0_full_noveto.ran.fits'.format(tr))
    else:
        prog = 'dark'
        if 'BGS' in tr:
            prog = 'bright'
        data_full = fitsio.read(indir+'{}_0_full_noveto.ran.fits'.format(prog))
    maxp = 3200
    if 'BGS' in tr:
        maxp = 2100
    if 'QSO' in tr:
        maxp = 3400
    area_tot = len(data_full)/2500
    sel_gh = data_full['GOODHARDLOC']
    area_bh = len(data_full[~sel_gh])/2500
    sel_pri = data_full['PRIORITY'] <= maxp#data_full['GOODPRI']
    area_bp = len(data_full[~sel_pri])/2500
    if mainp.reccircmasks is not None:
        for maskfn in mainp.reccircmasks:
            mask = common.maskcircandrec(data_full,maskfn)
            data_full = data_full[~mask]
    ebits = mainp.ebits
    if ebits == 'lrg_mask':
        sel = data_full['lrg_mask'] == 0
        data_full = data_full[sel]
    else:
        data_full = common.cutphotmask(data_full,ebits)
    area_im = area_tot - len(data_full)/2500
    print(tr,'total area:'+str(area_tot),'\n area and fraction in hardware mask:'+str(area_bh),str(area_bh/area_tot),'\n area and fraction in priority mask:'+str(area_bp),str(area_bp/area_tot),'\n area and fraction in imaging mask:'+str(area_im),str(area_im/area_tot))



