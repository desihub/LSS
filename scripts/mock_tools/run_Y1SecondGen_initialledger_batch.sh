SeconGenVer=AbacusSummit_v3_1 #AbacusSummit
for j in {0..24}
do
#j=0
echo $j
#echo $DESI_ROOT/survey/catalogs/Y1/mocks/SecondGenMocks/$SeconGenVer/forFA$j.fits $DESI_ROOT/survey/catalogs/Y1/mocks/SecondGenMocks/$SeconGenVer/altmtl$j/initled
python /pscratch/sd/a/acarnero/codes/LSS/scripts/mock_tools/run_mtl_ledger.py $DESI_ROOT/survey/catalogs/Y1/mocks/SecondGenMocks/$SeconGenVer/forFA$j.fits $DESI_ROOT/survey/catalogs/Y1/mocks/SecondGenMocks/$SeconGenVer/altmtl$j/initled DARK
#python /pscratch/sd/a/acarnero/codes/LSS/scripts/mock_tools/run_mtl_ledger.py $DESI_ROOT/survey/catalogs/Y1/mocks/SecondGenMocks/$SeconGenVer/forFA$j.fits $DESI_ROOT/survey/catalogs/Y1/mocks/SecondGenMocks/$SeconGenVer/altmtl$j/initled DARK
done
