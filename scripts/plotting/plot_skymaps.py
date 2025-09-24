import fitsio
from astropy.table import Table
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import sys,os
from desi_y1_plotting.kp3 import KP3StylePaper
style = KP3StylePaper()

#set LSSCODE, e.g. via export LSSCODE=$HOME on the command line if that is where you cloned the repo 
#sys.path.append(os.environ['LSSCODE']+'/LSS/py')
from LSS import common_tools as common


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

@mpl.rc_context(style._rcparams)
def plot_comp_sindec_panel(tracer,vmin=0,vmax=1,column='COMP_TILE',size_fac=1,addedges='y',fontsize=None):
    #
    datfile = indir+ tracer + '_full_HPmapcut.dat.fits'
    vmin = 0.
    vmax = 1.
    title  = tracer[:3] + ' fiber assignment completeness'
    
    #read in data and convert coordinates
    if fontsize is not None:
        mpl.rcParams.update({'font.size': int(fontsize)})
    cols = ['RA','DEC',column,'ZWARN','DELTACHI2','Z_not4clus']
    if tracer[:3] == 'ELG':
        cols.append('o2c')

    dt = Table(fitsio.read(datfile,columns=cols))
    sel_gz = common.goodz_infull(tracer[:3],dt)
    sel_obs = dt['ZWARN'] != 999999
    dt = dt[sel_obs & sel_gz]
    ras = dt['RA']
    ras[ras > 300] -= 360
    sindec = np.sin(dt['DEC']*np.pi/180.)
    yr = (np.max(sindec)-np.min(sindec))*1.05
    xr = 360*1.1/90#(np.max(ra)-np.min(ra))*1.1/90
    if size_fac == 2:
        xfac = 2.*size_fac
        yfac = 2.3*size_fac
    if size_fac == 1:
        xfac = 2.*size_fac
        yfac = 2.5*size_fac #if size_fac is lower, label text takes up more relative room
       
    fig = plt.figure(figsize=(xr*xfac, yr*yfac))
    ps = 0.1
    ax = fig.add_subplot(111)
    ax.set_aspect(90) #this is the correct aspect to use given sin(dec) of 0,1 corresponds to 0,90 degrees

    mp = plt.scatter(ras,sindec,c=dt[column],edgecolor='none',vmin=vmin,vmax=vmax,s=ps)

    #read in footprint outline
    if addedges == 'y':
        out = Table.read("/dvs_ro/cfs/cdirs/desi/users/raichoor/footprints/desi-boundaries.ecsv")
        if tracer[:3] == 'BGS':
            selprog = out['PROGRAM'] == 'BRIGHT'
        else:    
            selprog = out['PROGRAM'] == 'DARK'
        out = out[selprog]    
        caps= ['NGC','SGC']
        ras = out['RA']
        ras[ras > 300] -= 360
        sindec = np.sin(out['DEC']*np.pi/180.)

        for cap in caps:
            sel = out['CAP'] == cap
            #plt.fill(ras[sel],sindec[sel],color='lightgray',hatch="x",edgecolor='w',zorder=0)
            plt.plot(ras[sel],sindec[sel],color='k',zorder=1)
    
    
    plt.colorbar(mp,fraction=0.05)#,shrink=2/2.3) 
    plt.xlabel('R.A. [deg]')#,fontsize=14)
    decvals = np.array([-30,-15,0,15,30,60,90])
    sindecvals = np.sin(np.radians(decvals))
    plt.yticks(sindecvals,decvals)
    plt.ylabel('Dec. [deg]')#,fontsize=14)
    plt.xlim(300,-60.)
    plt.ylim(sindecvals[0],sindecvals[-1])
    plt.title(title)#,fontsize=14)
    #plt.grid()
    plt.tight_layout()
    plt.savefig(outdir+tracer+'_comptile.png', bbox_inches='tight')    

    return

@mpl.rc_context(style._rcparams)
def plot_ntile_sindec_panel(tracer,column='NTILE',size_fac=1,title='',addedges='y',fontsize=None):
    if fontsize is not None:
        mpl.rcParams.update({'font.size': int(fontsize)})

    #
    #datfile = args.indir+ tracer + '_0_full_noveto.ran.fits'
    datfile = indir+ tracer + '_0_full_HPmapcut.ran.fits'
    
    #read in data and convert coordinates
    cols = ['RA','DEC',column]

    dt = Table(fitsio.read(datfile,columns=cols))
    ras = dt['RA']
    ras[ras > 300] -= 360
    sindec = np.sin(dt['DEC']*np.pi/180.)
    yr = (np.max(sindec)-np.min(sindec))*1.05
    xr = 360*1.1/90#(np.max(ra)-np.min(ra))*1.1/90
    if size_fac == 2:
        xfac = 2.*size_fac
        yfac = 2.3*size_fac
    if size_fac == 1:
        xfac = 2.*size_fac
        yfac = 2.5*size_fac #if size_fac is lower, label text takes up more relative room
       
    fig = plt.figure(figsize=(xr*xfac, yr*yfac))
    ps = 0.1
    ax = fig.add_subplot(111)
    ax.set_aspect(90) #this is the correct aspect to use given sin(dec) of 0,1 corresponds to 0,90 degrees
    cmap = plt.get_cmap('jet', 7)
    vmax = 7.5
    ticks=np.arange(1,8 )
    if tracer[:3] == 'BGS':
        vmax = 4.5
        cmap = plt.get_cmap('RdBu_r', 4)
        ticks=np.arange(1,5)
    mp = plt.scatter(ras,sindec,c=dt[column],edgecolor='none',cmap=cmap,vmin=.5,vmax=vmax,s=ps)
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(mp,fraction=0.05,ticks=ticks)#,shrink=2/2.3) 
    plt.xlabel('R.A. [deg]')#,fontsize=14)
    decvals = np.array([-30,-15,0,15,30,60,90])
    sindecvals = np.sin(np.radians(decvals))
    plt.yticks(sindecvals,decvals)
    plt.ylabel('Dec. [deg]')#,fontsize=14)
    plt.xlim(300,-60.)
    plt.ylim(sindecvals[0],sindecvals[-1])
    plt.title(title)#,fontsize=14)
    #plt.grid()
    #read in footprint outline
    if addedges == 'y':
        out = Table.read("/dvs_ro/cfs/cdirs/desi/users/raichoor/footprints/desi-boundaries.ecsv")
        if tracer[:3] == 'BGS':
            selprog = out['PROGRAM'] == 'BRIGHT'
        else:    
            selprog = out['PROGRAM'] == 'DARK'
        out = out[selprog]    
        caps= ['NGC','SGC']
        ras = out['RA']
        ras[ras > 300] -= 360
        sindec = np.sin(out['DEC']*np.pi/180.)

        for cap in caps:
            sel = out['CAP'] == cap
            #plt.fill(ras[sel],sindec[sel],color='lightgray',hatch="x",edgecolor='w',zorder=0)
            plt.plot(ras[sel],sindec[sel],color='k',zorder=1)

    plt.tight_layout()
    plt.savefig(outdir+tracer+'_ntile.png', bbox_inches='tight',dpi=500)    

    return



if args.tracers is not None:
    for tracer in list_tracer:
        # Read in completeness map data
        plot_comp_sindec_panel(tracer,fontsize=args.fontsize)

if args.mk_ntile_dark == 'y':
    plot_ntile_sindec_panel('QSO',title='number of overlapping dark time tiles',fontsize=args.fontsize)

if args.mk_ntile_bright == 'y':
    plot_ntile_sindec_panel('BGS_BRIGHT',title='number of overlapping bright time tiles',fontsize=args.fontsize)


