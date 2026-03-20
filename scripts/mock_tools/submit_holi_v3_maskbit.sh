#!/bin/bash

#submit batch scripts for maskbits in 25 mock increments

for i in {450..975..25}; do
    let maxi=$i+25
    echo qsolrgelg$maxi.sbatch
    #python $HOME/LSScode/LSS/scripts/mock_tools/prep_bitmask.py --qsomin $i --qsomax $maxi --lrgmin $i --lrgmax $maxi --elgmin $i --elgmax $maxi --subname qsolrgelg$maxi
    #sbatch $HOME/BRICKMASKcode/brickmask/qsolrgelg$maxi.sbatch
done

