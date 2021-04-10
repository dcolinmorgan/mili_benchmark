#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import os, glob
import pandas as pd
import sys
from sklearn import metrics
import getopt
import tests
import multiprocessing as mp# to run FIMO in parallel
from functools import partial

"""
Example:
  source /proj/relibs/relib00/conda/bin/activate
  source activate mypy3
  python mili_benchmark/src/python/maxHits.py
"""

def maxHit(i,chunk_size):
	print(i)
	data=pd.read_csv('../../d/tmp/redmo/bench/slop100_Motif_BS_Ch.txt',sep='\t',header=None,nrows=chunk_size,skiprows=chunk_size*i)
	data1=data.groupby([0,1,2,10]).max()
	data1[3]=data1[3].round(4)
	data1[[3,4,8,9]].to_csv('../../d/tmp/redmo/bench/slop100_Motif_BS_Chmax2.txt',sep='\t',index=True,header=False,mode='a')



chunk_size=1000000
# for i in range(np.round(3505594261/chunk_size).astype('int')): ## will NOT be the same length, despite the want
xx=range(np.round(3505594261/chunk_size).astype('int'))

pool = mp.Pool(20)
res  = pool.map(partial(maxHit,chunk_size=chunk_size),xx)


