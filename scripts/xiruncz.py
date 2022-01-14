#make sure to type these two commands:
#export OMP_NUM_THREADS=64
#module load gsl
#python xiruncz.py --type ELG_HIP
import subprocess
import sys
import argparse
import os
#sys.path.append('../py')
#import LSS.mkCat_singletile.xitools as xt
#import LSS.SV3.xitools as xt


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is desi catalog directory",default='/global/cfs/cdirs/desi/survey/catalogs')
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='everest')
parser.add_argument("--survey",help="e.g., SV3 or main",default='SV3')
parser.add_argument("--nran",help="number of random files to combine together (1-18 available)",default=10)

args = parser.parse_args()

type = args.type
basedir = args.basedir
version = args.version
specrel = args.verspec
survey = args.survey
nran = int(args.nran)

if survey == 'SV3':
    import LSS.SV3.xitools as xt
if survey == 'main':
    import LSS.main.xitools as xt

lssdir = basedir+'/'+survey+'/LSS/'+specrel+'/LSScats/'
#dirout = svdir+'LSScats/'+version+'/'

zmask = ['']
minn = 0

subt = None
if type == 'LRGAlltiles' or type == 'LRGAlltiles_main':
    zl = [0.32,0.6,0.8,1.05,1.3]
    #minn = 2
    #zmin=0.32
    #zmax=1.05

if type == 'LRG':
    zl = [0.4,0.6,0.8,1.1]
#    minn = 5
    #zmin=0.32
    #zmax=1.05


if type == 'LRG_OPT':
    subt = type
    zmin=0.6
    zmax=1.
    type = 'LRG'

if type == 'LRG_IR':
    subt = type
    zmin=0.6
    zmax=1.
    type = 'LRG'


if type[:3] == 'ELG':# or type == 'ELG_HIP':
    #minn = 5
    zl = [0.8,1.1,1.5]
    #zmask = ['','_zmask']
    
    #zmin = 0.8
    #zmax = 1.6

#if type == 'ELG_HIP':
#    zmin = 0.8
#    zmax = 1.6
if type == 'ELG_HIP16':
    minn = 5
    zl = [1,1.6]
    type = 'ELG_HIP'

if type == 'ELG16':
    minn = 5
    zl = [1,1.6]
    type = 'ELG'


if type == 'ELGlz':
    zmin = 0.6
    zmax = 0.8
    type = 'ELG'

if type == 'ELGmz':
    zmin = 0.8
    zmax = 1.1
    type = 'ELG'

if type == 'ELGhz':
    zmin = 1.1
    zmax = 1.6
    type = 'ELG'    

if type == 'ELGmhz':
    zmin = 0.6
    zmax = 1.497
    type = 'ELG'    

if type == 'ELGhz497':
    zmin = 1.1
    zmax = 1.497
    type = 'ELG'

if type == 'QSO':
    zl = [0.8,1.1,1.5,2.1]
    #zmin = 1.
    #zmax = 2.1

if type == 'QSOhiz':
    zmin = 1.6
    zmax = 2.1
    type = 'QSO'

if type == 'QSOlya':
    #zmin = 2.1
    #zmax = 3.5
    zl = [2.1,3.5]
    type = 'QSO'


if type == 'QSO_RF_4PASS':
    subt = type
    zmin = 1.6
    zmax = 2.1
    type = 'QSO'

if type == 'ELG_FDR_GFIB':
    subt = type
    zmin = 1.1
    zmax = 1.6
    type = 'ELG'

if type[:3] == 'BGS':
    #minn = 2
    zl = [0.1,0.3,0.5]
    #zmin = 0.1
    #zmax = 0.5 
    
if type == 'BGS_hiz':
    zmin = 0.3
    zmax = 0.5
    type = 'BGS_ANY'        

ranwt1=False

regl = ['_N','_S']

if survey == 'main':
    regl = ['_DN','_DS','_N','_S']

for i in range(0,len(zl)):
    if i == len(zl)-1:
        zmin=zl[0]
        zmax=zl[-1]
    else:
        zmin = zl[i]
        zmax = zl[i+1]
    print(zmin,zmax)
    for zma in zmask:
        for reg in regl:
            xt.prep4czxi(type,zmin,zmax,nran=nran,indir=lssdir,ver=version,minn=minn,reg=zma+reg,outdir=os.environ['CSCRATCH']+'/cz/',ranwt1=ranwt1,subt=subt)
            subprocess.run(['chmod','+x','czpc.sh'])
            subprocess.run('./czpc.sh')
            fa = ''
            if ranwt1:
                fa = 'ranwt1'
            if subt is not None:
                fa += subt    
            xt.calcxi_dataCZ(type,zmin,zmax,minn=minn,reg=zma+reg,ver=version,fa=fa)


        xt.prep4czxi(type,zmin,zmax,nran=nran,indir=lssdir,ver=version,minn=minn,reg=zma,outdir=os.environ['CSCRATCH']+'/cz/',ranwt1=ranwt1,subt=subt)
        subprocess.run(['chmod','+x','czpc.sh'])
        subprocess.run('./czpc.sh')
        fa = ''
        if ranwt1:
            fa = 'ranwt1'
        if subt is not None:
            fa += subt    
        xt.calcxi_dataCZ(type,zmin,zmax,minn=minn,ver=version,fa=fa,reg=zma)

