#!/bin/bash
#

python /pscratch/sd/a/acarnero/codes/LSS/bin/runAltMTLParallel.py --altMTLBaseDir=/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1//altmtl1_R128_ELGs/ --obscon=DARK --survey=main --verbose --mock --targfile=/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v3_1/forFA1_toBit.fits --multiDate --universe=$1
