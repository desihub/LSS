import os
import sys
import argparse
import subprocess
import numpy as np

from   astropy.io    import ascii
from   astropy.table import Table

home  = os.environ['HOME']
sys.path.append('{}/DESI'.format(home))

from   ddp           import tmr_DDP1, tmr_DDP2, tmr_DDP3
from   delta8_limits import delta8_tier, d8_limits
from   findfile      import findfile

parser = argparse.ArgumentParser(description='Generate Summary Stats')
parser.add_argument('-s', '--survey', help='Select survey', default='gama')
parser.add_argument('-v', '--version', help='Select version', default='GAMA4')

args    = parser.parse_args()
survey  = args.survey
version = args.version

fpath            = findfile(ftype='ddp', survey=survey, version=version)
dat              = Table.read(fpath)
names            = ['ZMIN', 'ZMAX', 'VZ', 'DENS']

tmr_DDPs         = np.array([tmr_DDP1, tmr_DDP2, tmr_DDP3])
tmr_ngold_ratios = np.array([1.002, 0.591, 0.097])
tmr_ddp_ratios   = np.array([1., 0.589, 0.097])
tmr_Nd8          = np.array([2.18, 2.31, 2.72, 9.48, 16.1, 17.3, 16.2, 7.21, 7.57])

result   = Table()
rows     = []

print('\n\n')

for ddp, tmr_DDP, tmr_ngold_ratio, tmr_ddp_ratio in zip(np.arange(1, 4, 1), tmr_DDPs, tmr_ngold_ratios, tmr_ddp_ratios):
    row = [ddp, tmr_DDP[0], tmr_DDP[1]]
    
    for col in names:
        row += [dat.meta['DDP{}_{}'.format(ddp, col)]]
        
    row += [dat.meta['DDP{}{}'.format(ddp, 'ZLIMS_NGAL')]]
    row += [1. * dat.meta['DDP{}_NGAL'.format(ddp)] / dat.meta['GOLD_NGAL']]
    row += [tmr_ngold_ratio]
        
    row = tuple(row)         
    rows.append(row)

names  = ['DDP', 'MIN_M', 'MAX_M'] + names + ['ZLIMS_NGAL', 'N_NGOLD', 'N_NREF_TMR']
result = Table(rows=rows, names=names)

result['ZLIMS_NGAL']      = result['ZLIMS_NGAL'] / 10**3
result['VZ']              = result['VZ'] / 10**6
result['DENS']            = result['DENS'] / 10**-3
result['N_GAL/N_GAL_MAX'] = 1. * result['ZLIMS_NGAL'] / max(result['ZLIMS_NGAL'])

for name in result.dtype.names:
    result[name]          = np.round(result[name], 3)

opath = 'tables/Tab2.fits'
result.write(opath, format='fits', overwrite=True)
    
result.rename_column('MIN_M', '$M_{Min.}$')
result.rename_column('MAX_M', '$M_{Max.}$')
result.rename_column('ZLIMS_NGAL', '$N_{GAL} / 10^3$')
result.rename_column('VZ', '$V_{DDP}$ / 10^6$')
result.rename_column('DENS', '$\\rho_{DDP} / 10^{-3}$')
result.rename_column('N_NGOLD', '$N/N_{GOLD}$')
result.rename_column('N_NREF_TMR', '$N/N_{REF}_{TMR}$')
result.rename_column('N_GAL/N_GAL_MAX', '$N_{GAL}/N_{GALMAX}$')

result.pprint()

# https://arxiv.org/pdf/1409.4681.pdf
ascii.write(result, 'tables/Tab2.tex', Writer=ascii.Latex, latexdict=ascii.latex.latexdicts['AA'], overwrite=True)

rows = []

rpath = findfile(ftype='randoms_bd_ddp_n8', survey=survey, version=version, field='G9')
rand  = Table.read(rpath)

print('\n\n')

for idx in np.arange(9):    
    nd8 = 0 
    
    for field in ['G9', 'G12', 'G15']:
        
        fpath = findfile(ftype='ddp_n8_d0', survey=survey, version=version, field=field, utier=idx)
        dat   = Table.read(fpath)
                
        nd8  += len(dat) / 1.e3
        
    rows.append(('d{}'.format(idx), d8_limits[idx][0], d8_limits[idx][1], nd8, rand.meta['DDP1_d{}_VOLFRAC'.format(idx)]))
    
print('\n\n')

# Generate Table 3 of McNaught-Roberts (2014).
result = Table(rows=rows, names=['Label', 'Min_{d8}', 'Max_{d8}', 'N_{d8} [1e3]', 'fd8']) #, 'N_{d8}/N_{max}'])

# TODO: CHECK THE MATHS BY HAND
result['N_{d8} / N_{max}'] = result['N_{d8} [1e3]'] / max(result['N_{d8} [1e3]'])

# BUG: result['TMR N/N_{max}'] = tmr_Nd8 / tmr_Nd8[5]
result['TMR N/N_{max}'] = tmr_Nd8 / tmr_Nd8[4]

for col in result.itercols():
    if col.info.dtype.kind == 'f':        
        np.around(col, decimals=3, out=col)

result.pprint()

opath = 'tables/Tab3.fits'
result.write(opath, format='fits', overwrite=True)

# https://arxiv.org/pdf/1409.4681.pdf
ascii.write(result, 'tables/Tab3.tex', Writer=ascii.Latex, latexdict=ascii.latex.latexdicts['AA'], overwrite=True)

# subprocess.run('cd {}/tables; pdflatex compile_tables.tex'.format(os.environ['CODE_ROOT']), shell=True, check=True)

print('\n\nDone.\n\n')
