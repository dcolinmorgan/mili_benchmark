#!/usr/bin/python
# python mili_benchmark/src/python/test.py
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import os,gc,glob
# import jax.numpy as np
# import s3fs
import seaborn as sns
from psutil import *
import scipy.io
from scipy.stats import zscore
import datetime
from sklearn import metrics

import netZooPy
from netZooPy.panda import Panda

# unzip bench.zip -d bench
# unzip slop100_bench.zip -d bench_100

def test_panda(motif):#,size):
    panda_objC=netZooPy.panda.Panda(expression_file='data/Hugo_exp1_lcl.txt',
                          motif_file='data/bench_100/'+motif,
                          ppi_file='data/ppi2015_freezeCellLine.txt',
                          computing='cpu',modeProcess='intersection',save_memory=False,save_tmp=False,
                          precision='single',keep_expression_matrix=False)

    if not os.path.exists('data/bench_100/output'):
        os.makedirs('data/bench_100/output/')
    panda_objC.export_panda_results.to_csv('data/bench_100/output/'+motif,sep='\t',header=True,index=False)
    return panda_objC.export_panda_results


test_panda('slop100_max_pwm_uniqMotif_BS_Chmax_refseq.txt')
test_panda('slop100_max_wgbs_uniqMotif_BS_Chmax_refseq.txt')

test_panda('slop100_min_pwm_uniqMotif_BS_Chmax_refseq.txt')
test_panda('slop100_min_wgbs_uniqMotif_BS_Chmax_refseq.txt')

test_panda('slop100_mean_pwm_uniqMotif_BS_Chmax_refseq.txt')
test_panda('slop100_mean_wgbs_uniqMotif_BS_Chmax_refseq.txt')



