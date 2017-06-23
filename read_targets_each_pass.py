import os
from astropy.io import fits
from astropy.table import Table, vstack, unique
from astropy.table import Column, join
import numpy as np

tile_info_file = '/project/projectdirs/desi/software/edison/desimodel/master/data/footprint/desi-tiles.fits'
pass_data = fits.open(tile_info_file)[1].data
mask_desi = pass_data['IN_DESI']==1
pass_desi =pass_data[mask_desi]
mask_dark = pass_desi['PROGRAM']=='DARK'
pass_desi_dark = pass_desi[mask_dark]
len_pass_file = pass_desi_dark.shape[0]
t = Table( names=('TILEID', 'ELG', 'LRG', 'QSO'))

i=0
for filename in os.listdir(os.getcwd()):
    start = filename.find('tile_') + 5
    end = filename.find('.fits', start)
    tile_id =filename[start:end].lstrip("0")
    tile_data = fits.open(filename)[1].data
    mask_ELG = tile_data['DESI_TARGET']==2
    elg_t = tile_data[mask_ELG]
    elg_id = np.array(elg_t['TARGETID'])
    if i==0:
        new_elg_array = elg_id
        no_ELG = len(new_elg_array)
        old_ELG = no_ELG
    else:
        new_elg_array= np.concatenate((new_elg_array, elg_id), axis=0)
        elg_unique = np.unique(new_elg_array)
        no_ELG = len(elg_unique)-old_ELG
        old_ELG = len(elg_unique)
    mask_LRG = tile_data['DESI_TARGET']==1
    lrg_t = tile_data[mask_LRG]
    lrg_id = np.array(lrg_t['TARGETID'])
    if i==0:
        new_lrg_array =lrg_id
        no_LRG = len(new_lrg_array)
        old_LRG= no_LRG
    else:
        new_lrg_array= np.concatenate((new_lrg_array, lrg_id), axis=0)
        lrg_unique = np.unique(new_lrg_array)
        no_LRG = len(lrg_unique)-old_LRG
        old_LRG= len(lrg_unique)
    mask_QSO = tile_data['DESI_TARGET']==4
    qso_t = tile_data[mask_QSO]
    qso_id = np.array(qso_t['TARGETID'])
    if i==0:
        new_qso_array =qso_id
        no_QSO = len(new_qso_array)
        old_QSO= no_QSO
    else:
        new_qso_array= np.concatenate((new_qso_array, qso_id), axis=0)
        qso_unique = np.unique(new_qso_array)
        no_QSO = len(qso_unique)-old_QSO
        old_QSO= len(qso_unique)
    t.add_row([tile_id ,no_ELG,no_LRG, no_QSO])
    i+=1
nt = unique(t, keys=['TILEID'], keep='last')
joined_table = join(nt, pass_desi, keys='TILEID')
new_table = joined_table.group_by('PASS')
ink =new_table.groups.indices

new_tab_survey_dark=new_table[ink[0]:ink[1]].group_by('PROGRAM')
sumvals = new_tab_survey_dark.groups.aggregate(np.sum)
avvals = new_tab_survey_dark.groups.aggregate(np.mean)
output_tab_D = Table([sumvals['ELG'], sumvals['LRG'], sumvals['QSO'], sumvals['PASS'],sumvals['PROGRAM']], names=('ELG', 'LRG', 'QSO', 'PASS', 'PROGRAM'))
for i in range(1,len(ink)-1):
    new_tab_survey_dark=new_table[ink[i]:ink[i+1]].group_by('PROGRAM')
    sumvals = new_tab_survey_dark.groups.aggregate(np.sum)
    avvals = new_tab_survey_dark.groups.aggregate(np.mean)
    output_tab_D.add_row([sumvals['ELG'], sumvals['LRG'], sumvals['QSO'], avvals['PASS'], sumvals['PROGRAM']])
print(output_tab_D)

         

mtl_data = fits.open('/global/project/projectdirs/desi/datachallenge/quicksurvey2017/output/dark/0/mtl.fits')[1].data
