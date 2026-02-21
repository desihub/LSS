
source /global/common/software/desi/desi_environment.sh main

export LSSDIR=$HOME
export SYSNETDIR=$HOME/desicode
export LSSBASE=/global/cfs/cdirs/desi/survey/catalogs/
PYTHONPATH=$PYTHONPATH:$LSSDIR/LSS


srun python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type QSO --prepsysnet y --imsys_zbin y --fulld n --survey Y1 --verspec iron --version $1 &

srun python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type BGS_BRIGHT-21.5 --prepsysnet y --imsys_zbin y --fulld n --survey Y1 --verspec iron --version $1 &

wait


srun $LSSDIR/LSS/scripts/run_sysnet.sh N QSO0.8_1.3 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &
srun $LSSDIR/LSS/scripts/run_sysnet.sh S QSO0.8_1.3 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &

srun $LSSDIR/LSS/scripts/run_sysnet.sh N QSO1.3_2.1 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &
srun $LSSDIR/LSS/scripts/run_sysnet.sh S QSO1.3_2.1 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &

srun $LSSDIR/LSS/scripts/run_sysnet.sh N QSO2.1_3.5 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &
srun $LSSDIR/LSS/scripts/run_sysnet.sh S QSO2.1_3.5 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &

srun $LSSDIR/LSS/scripts/run_sysnet.sh N BGS_BRIGHT-21.50.1_0.4 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &
srun $LSSDIR/LSS/scripts/run_sysnet.sh S BGS_BRIGHT-21.50.1_0.4 true false 1024 0.003 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &

wait


srun $LSSDIR/LSS/scripts/run_sysnet.sh N QSO0.8_1.3 false true 1024 0.01 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &
srun $LSSDIR/LSS/scripts/run_sysnet.sh S QSO0.8_1.3 false true 1024 0.01 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &

srun $LSSDIR/LSS/scripts/run_sysnet.sh N QSO1.3_2.1 false true 1024 0.01 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &
srun $LSSDIR/LSS/scripts/run_sysnet.sh S QSO1.3_2.1 false true 1024 0.01 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &

srun $LSSDIR/LSS/scripts/run_sysnet.sh N QSO2.1_3.5 false true 1024 0.01 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &
srun $LSSDIR/LSS/scripts/run_sysnet.sh S QSO2.1_3.5 false true 1024 0.01 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &

srun $LSSDIR/LSS/scripts/run_sysnet.sh N BGS_BRIGHT-21.50.1_0.4 false true 1024 0.01 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &
srun $LSSDIR/LSS/scripts/run_sysnet.sh S BGS_BRIGHT-21.50.1_0.4 false true 1024 0.01 dnnp pnll $1 $LSSBASE/Y1/LSS/iron/LSScats/ &

wait


srun python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type QSO --add_sysnet y --imsys_zbin y --fulld n --survey Y1 --verspec iron --version $1 &

srun python scripts/main/mkCat_main.py --basedir /global/cfs/cdirs/desi/survey/catalogs/ --type BGS_BRIGHT-21.5 --add_sysnet y --imsys_zbin y --fulld n --survey Y1 --verspec iron --version $1 &

wait


srun python scripts/validation/validation_improp_full.py --tracer BGS_BRIGHT-21.5 --version $1 --weight_col WEIGHT_SN &

srun python scripts/validation/validation_improp_full.py --tracer QSO --version $1 --weight_col WEIGHT_SN &

wait