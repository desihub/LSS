import os

path = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS'
for i in range(6,25):
#    os.system('cp -R %s/altmtl%d/initled/main %s/altmtl%d/Univ000/.' %(path, i, path, i))
    os.system('cp %s/altmtl5/Univ000/*.ecsv %s/altmtl%d/Univ000/.' %(path, path, i))
