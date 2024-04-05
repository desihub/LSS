OBSCON=DARK  #BRIGHT

SeconGenVer=AbacusSummit_v4_1 #AbacusSummitBGS_v2
for j in {0..24}
do
echo $j
python /pscratch/sd/a/acarnero/codes/LSS/scripts/mock_tools/run_mtl_ledger.py $DESI_ROOT/survey/catalogs/Y1/mocks/SecondGenMocks/$SeconGenVer/forFA$j.fits $DESI_ROOT/survey/catalogs/Y1/mocks/SecondGenMocks/$SeconGenVer/altmtl$j $OBSCON
done
