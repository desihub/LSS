#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import logging

import numpy as np
from astropy.table import Table, vstack
from matplotlib import pyplot as plt

from pycorr import TwoPointCorrelationFunction, TwoPointEstimator, KMeansSubsampler, utils, setup_logging

njack = '60'
trs = ['ELG_LOPnotqso','QSO','LRG','BGS_BRIGHT']
bsl = [1,2,4,5,10]
dirxi = '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/LSScats/test/xi/smu/'
xit = 'poles'
for tr in trs:
    if tr == 'ELG_LOPnotqso':
        zws = ['0.8_1.6','0.8_1.1','1.1_1.6']
    if tr == 'QSO':
        zws = ['0.8_1.1','0.8_2.1lowz','1.1_1.6','1.6_2.1','2.1_3.5','0.8_3.5']
    if tr == 'LRG':
        zws = ['0.4_0.6','0.6_0.8','0.8_1.1','0.4_1.1']
    if tr == 'BGS_BRIGHT':
        zws = ['0.1_0.3','0.3_0.5','0.1_0.5']
    for zw in zws:
        result_N = pycorr.TwoPointCorrelationFunction.load(dirxi+'allcounts_'+tr+'_N_'+zw+'_default_FKP_lin_njack'+njack+'.npy')
        result_S = pycorr.TwoPointCorrelationFunction.load(dirxi+'allcounts_'+tr+'_S_'+zw+'_default_FKP_lin_njack'+njack+'.npy')
        result_NS = result_N.normalize() + result_S.normalize()
        fn = dirxi+'allcounts_'+tr+'_NScomb_'+zw+'_default_FKP_lin_njack'+njack+'.npy'
        result_NS.save(fn)
        for bs = bsl:
            rebinned = result_NS[:(result_NS.shape[0]//bs)*bs:bs]
            fn_txt = dirxi+'xi'+xit+'_'+tr+'_NScomb_'+zw+'_default_FKP_lin_njack'+njack+'.txt'
            rebinned.save_txt(fn_txt, ells=(0, 2, 4))


