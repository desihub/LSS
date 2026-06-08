import os
for i in range(7,25):
    print(i)
    os.system('cp -R /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS/altmtl%d/initled/main /global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummitBGS/altmtl%d/Univ000/.' %(i,i))

