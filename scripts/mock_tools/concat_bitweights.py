from astropy.table import Table,vstack
import numpy as np

def concatenateBWFiles(BWDir, hpList, survey = 'main', obscon = 'dark', OFName = '{0}bw-{1}-allTiles.fits', skipFailures = False, overwrite = False):
    BWBaseFN = BWDir + '{0}/{1}/{0}bw-{1}-hp-{2:d}.fits'
    assert(len(hpList) > 1)
    AllBWFiles = Table.read(BWBaseFN.format(survey, obscon, hpList[0]), hdu = 1)
    notCompleted = []    
    for hp in hpList[1:]:
        print(BWBaseFN.format(survey, obscon, hp))
        try:
            BWTemp = Table.read(BWBaseFN.format(survey, obscon, hp), hdu = 1)
        except Exception as e:
            if skipFailures:
                notCompleted.append(hp)
                continue
            else:
                raise(e)
                
        AllBWFiles = vstack([AllBWFiles, BWTemp])
#    print(', '.join(map(str, notCompleted)))
    AllBWFiles.write(BWDir + '{0}/{1}/'.format(survey, obscon) + OFName.format(survey, obscon), format = 'fits', overwrite = overwrite)




hpL_file = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1/altmtl1/initled/hpxlist_dark.txt'
HPList = np.array(open(hpL_file, 'r').readlines()[0].split(',')).astype(int)

BW_dir = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1/altmtl1_R64/BitweightFiles/'


concatenateBWFiles(BW_dir, HPList, skipFailures=False, overwrite=True)
