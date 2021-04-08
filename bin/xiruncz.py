#make sure to type these two commands:
#export OMP_NUM_THREADS=64
#module load gsl
import subprocess
import sys
import argparse
import os
sys.path.append('../py')
#import LSS.mkCat_singletile.xitools as xt
import LSS.SV3.xitools as xt


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is desi catalog directory",default='/global/cfs/cdirs/desi/survey/catalogs')
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')

args = parser.parse_args()

type = args.type
basedir = args.basedir
version = args.version


lssdir = basedir+'/SV3/LSS/LSScats/'
#dirout = svdir+'LSScats/'+version+'/'

subt = None
if type == 'LRG':
    zmin=0.32
    zmax=1.05

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


if type == 'ELG':
    zmin = 0.8
    zmax = 1.6

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
    zmin = 1.
    zmax = 2.1

if type == 'QSOhiz':
    zmin = 1.6
    zmax = 2.1
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
   

if type == 'BGS_ANY':
    zmin = 0.1
    zmax = 0.5 
    
if type == 'BGS_hiz':
    zmin = 0.3
    zmax = 0.5
    type = 'BGS_ANY'        

ranwt1=False
xt.prep4czxi(type,zmin,zmax,nran=10,indir=lssdir,ver=version,outdir=os.environ['CSCRATCH']+'/cz/',ranwt1=ranwt1,subt=subt)
subprocess.run(['chmod','+x','czpc.sh'])
subprocess.run('./czpc.sh')
fa = ''
if ranwt1:
    fa = 'ranwt1'
if subt is not None:
    fa += subt    
xt.calcxi_dataCZ(type,zmin,zmax,ver=version,fa=fa)
