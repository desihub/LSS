from desitarget import mtl
import sys
import glob
import os
import errno
from LSS.SV3 import altmtltools as amtl
from astropy.table import Table, vstack

def create_dirs(value):
    if not os.path.exists(value):
        try:
            os.makedirs(value, 0o755)
            print('made ' + value)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise



#python /pscratch/sd/a/acarnero/codes/LSS/scripts/mock_tools/run_mtl_ledger.py $DESI_ROOT/survey/catalogs/Y1/mocks/SecondGenMocks/$SeconGenVer/forFA$j.fits $DESI_ROOT/survey/catalogs/Y1/mocks/SecondGenMocks/$SeconGenVer/altmtl$j DARK

par=True
addextra=False

arg1 = sys.argv[1].replace('global','dvs_ro') #Input mock
arg2 = sys.argv[2] #Output path
obscon = sys.argv[3] #DARK or BRIGHT


initledger_path = os.path.join(arg2, 'initled')

altmtl_path = os.path.join(arg2, 'Univ000')
print('Running initial ledgers')
if par:
    mtl.make_ledger(arg1, initledger_path, obscon=obscon.upper(), numproc=32)
else:
    mtl.make_ledger(arg1, initledger_path, obscon=obscon.upper())

print('Creating list of tiles to be processed by AltMTL mock production')

path = os.path.join(initledger_path, 'main', obscon.lower())

ff = glob.glob(os.path.join(path, 'mtl-{obscon}-hp-*.ecsv'.format(obscon=obscon.lower())))


dd=[]

for f in ff:
    dd.append(int(f.split('hp-')[-1].split('.ecsv')[0]))
tosave = ','.join(map(str, sorted(dd)))

savepath = os.path.join(initledger_path, 'hpxlist_{obscon}.txt'.format(obscon = obscon.lower()))

ff = open(savepath, 'w')
ff.write(tosave)
ff.close()

print('saving list of HP ledgers in '+savepath)
path_to_altmtl = os.path.join(altmtl_path, 'main', obscon.lower())

print('Copying initial ledgers to altmtl directory ', path_to_altmtl)

create_dirs(path_to_altmtl)

os.system('cp %s %s' % (os.path.join(path,'*'), path_to_altmtl))

print('Creating tileTracker file and tilestatus file')


startDateShort = 19990101
endDate='20240418' #2024-04-18T00:00:00+00:00'
#20240418

os.system('cp /global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummit_v100/aux_data/mainsurvey-DARKobscon-TileTracker.ecsv %s' % os.path.join(altmtl_path, 'mainsurvey-DARKobscon-TileTracker.ecsv'))

'''
if ('T' in endDate) & ('-' in endDate):
    endDateShort = int(endDate.split('T')[0].replace('-', ''))
else:
    endDateShort = int(endDate)

amtl.makeTileTracker(altmtl_path, survey = 'main', obscon = 'DARK', startDate = startDateShort, endDate = endDateShort, overwrite = True)
'''

ztilefile = '/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv'

ztilefn = ztilefile.split('/')[-1]

if not os.path.isfile(os.path.join(altmtl_path, ztilefn)):
    amtl.processTileFile(ztilefile, os.path.join(altmtl_path, ztilefn), None, endDate)

if addextra:
    extratiles = Table.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4/aux_data/extra_{obscon}.ecsv'.format(obscon = obscon.lower()), format='ascii.ecsv')

    tileref = extratiles['TILEID'][-1]
#print(tileref)
    input_track = os.path.join(altmtl_path, 'mainsurvey-{OBSCON}obscon-TileTracker.ecsv'.format(OBSCON = obscon.upper()))
    tiles = Table.read(input_track, format='ascii.ecsv')
    tiles.meta['amtldir'] = altmtl_path

    if tiles['TILEID'][-1] != tileref:
        print('Adding extra tiles')
        newtable = vstack([tiles, extratiles])
        newtable.meta = tiles.meta
        newtable.write(input_track, overwrite=True)
    else:
        print('It has already been merged')
