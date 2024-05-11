#!/bin/bash
#python scripts/readwrite_pixel_bitmask.py --tracer lrg --input $1 --cat_type Y1EZmock
python scripts/mock_tools/ffa2clus_fast.py --mockver AbacusSummit_v4_2 --realization $1 --outloc prod --tracer QSO --overwrite y
python scripts/mock_tools/ffa2clus_fast.py --mockver AbacusSummit_v4_2 --realization $1 --outloc prod --tracer ELG_LOP --overwrite y
