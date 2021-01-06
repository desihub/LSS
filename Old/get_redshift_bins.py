from astropy.table import Table, join
from astropy.io import fits
import numpy as np

z_data = fits.open('/project/projectdirs/desi/datachallenge/quicksurvey2016/output/dark/0/zcat.fits')[1].data
mask = z_data['SPECTYPE']=='GALAXY'
z_gal = z_data[mask]

mtl_data = fits.open('/project/projectdirs/desi/datachallenge/quicksurvey2016/output/dark/0/mtl.fits')[1].data
mask = mtl_data['DESI_TARGET']!=4
gal_data = mtl_data[mask]

joined_table = join(z_gal, gal_data, keys='TARGETID')

mask_join = joined_table['DESI_TARGET']==2

masked_table = joined_table[mask_join]

bins = np.linspace(0, 3, 500)
bin_size = bins[1]-bins[0]
digitized = np.digitize(masked_table['Z'], bins)
max_bin = np.zeros(len(digitized))
min_bin = np.zeros(len(digitized))

#masked_table.groups.indices
masked_table['red bin no.'] = digitized
for i in range(1,len(digitized)):
    max_bin[i-1] = masked_table['red bin no.'][i-1]*bin_size
    min_bin[i-1] = (masked_table['red bin no.'][i-1]-1)*bin_size
masked_table['red bin max'] =max_bin
masked_table['red bin min'] =min_bin
grouped_table = masked_table.group_by('red bin no.')
index_of_group = grouped_table.groups.indices
bin_count = [(index_of_group[i]-index_of_group[i-1]) for i in range(1, len(index_of_group))]
norm_bin_count = np.divide(bin_count, sum(bin_count))

red_min = [(grouped_table['red bin min'][index_of_group[i-1]]) for i in range(1, len(index_of_group))]
red_max = [(grouped_table['red bin max'][index_of_group[i-1]]) for i in range(1, len(index_of_group))]

t = Table([norm_bin_count, red_min, red_max], names=('number', 'red_min', 'red_max'))
t.write('ELG_nz_zcat0.csv', format='csv')
