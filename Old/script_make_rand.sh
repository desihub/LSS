#!/bin/bash
cd /scratch1/scratchdirs/angela/rand5b/
tail -n +2  ELG*.csv >> allELG_0.csv
awk -F "\"*,\"*" '{print $9, $10}' allELG_0.csv > radec_ELG.csv
tail -n +2 radec_ELG.csv  > radec.csv
