#!/bin/bash

srun -n 1 python xirunpc.py --tracer ELG --survey SV3 --nran 1 --njack 10
srun -n 1 python xirunpc.py --tracer ELG LRG --survey SV3 --nran 1 --njack 10
srun -n 1 python xirunpc.py --tracer QSO --survey SV3 --nran 1 --njack 10 --bin_type log
srun -n 1 python xirunpc.py --tracer QSO --survey DA02 --verspec guadalupe --bin_type log --corr_type rppi --nran 1 --njack 10
srun -n 1 python xirunpc.py --tracer ELG --survey DA02 --verspec guadalupe --bin_type lin --corr_type rppi --nran 1 --njack 10 --region DN --zlim highz --weight_type angular --version 2
srun -n 1 python xirunpc.py --tracer ELG --survey SV3 --verspec everest --bin_type lin --corr_type rppi --nran 1 --njack 10 --region N --zlim highz --weight_type bitwise_angular --version 2.1 --split_ran_above 30


srun -n 1 python xirunpc.py --tracer ELG --survey SV3 --verspec everest --bin_type lin --corr_type smu --nran 1 --njack 10 --region N --weight_type bitwise_angular --version 2.1 --split_ran_above 30

srun -n 2 python xirunpc.py --tracer ELG --survey SV3 --verspec everest --bin_type lin --corr_type smu --nran 1 --njack 10 --region N --weight_type bitwise_angular --version 2.1 --split_ran_above 30

srun -n 64 python pkrun.py --tracer ELG QSO --survey SV3 --verspec everest --nran 1 --region N --weight_type bitwise_angular --version 2.1
srun -n 64 python pkrun.py --tracer ELG QSO --survey SV3 --verspec everest --nran 1 --region N --version 2.1 --nmesh 512

python recon.py --tracer ELG --survey DA02 --verspec guadalupe --version 2 --nran 4 --algorithm MG --convention reciso


srun -n 1 python xirunpc.py --tracer ELG_LOPnotqso --survey DA02 --verspec guadalupe --bin_type lin --corr_type smu --nran 4 --njack 0 --region N --weight_type RF_FKP --version 2 --split_ran_above 30

srun -n 2 python xirunpc.py --tracer ELG_LOPnotqso --survey DA02 --verspec guadalupe --bin_type lin --corr_type smu --nran 4 --njack 0 --region N --weight_type RF_FKP --version 2 --split_ran_above 30