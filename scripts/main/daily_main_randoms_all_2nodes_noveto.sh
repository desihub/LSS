#!/bin/bash

source /project/projectdirs/desi/software/desi_environment.sh master
PYTHONPATH=$PYTHONPATH:$HOME/LSS/py

./daily_main_randoms_type_2nodes_noveto.sh ELG n
./daily_main_randoms_type_2nodes_noveto.sh ELG_LOP n
./daily_main_randoms_type_2nodes_noveto.sh ELG_LOP y
./daily_main_randoms_type_2nodes_noveto.sh LRG n
./daily_main_randoms_type_2nodes_noveto.sh QSO n
./daily_main_randoms_type_2nodes_noveto.sh BGS_ANY n
./daily_main_randoms_type_2nodes_noveto.sh BGS_BRIGHT n
