#make sure to type these two commands:
#export OMP_NUM_THREADS=64
#module load gsl
import subprocess
import sys
import argparse
import os
sys.path.append('../py')
import LSS.mkCat_singletile.xitools as xt


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')

args = parser.parse_args()

type = args.type
basedir = args.basedir
version = args.version


lssdir = basedir+'/SV1/LSS/LSScats/'
#dirout = svdir+'LSScats/'+version+'/'


if type == 'LRG':
    zmin=0.6
    zmax=1.

if type == 'ELG':
    zmin = 0.6
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

xt.prep4czxi(type,zmin,zmax,nran=10,indir=lssdir,ver=version,outdir=os.environ['CSCRATCH']+'/cz/',ranwt1=False)
subprocess.run(['chmod','+x','czpc.sh'])
subprocess.run('./czpc.sh')
xt.calcxi_dataCZ(type,zmin,zmax,ver=version)
