OBSCON=DARK  #BRIGHT

SeconGenVer=AbacusSummit_v4_1 #AbacusSummitBGS_v2
for j in {1..3}
do
echo $j
python /pscratch/sd/a/acarnero/codes/LSS/scripts/mock_tools/run_mtl_ledger.py $DESI_ROOT/survey/catalogs/Y1/mocks/SecondGenMocks/EZmock/FFA/FFA_temp/SecondGenMocks/EZmock/forFA/forFA$j.fits $DESI_ROOT/survey/catalogs/Y1/mocks/SecondGenMocks/ALTMTL_EZmock/altmtl$j $OBSCON
done
