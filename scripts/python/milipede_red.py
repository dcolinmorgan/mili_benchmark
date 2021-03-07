import matplotlib.pyplot as plt
import numpy as np
import os, glob
import pandas as pd
# import multiprocessing as mp# to run FIMO in parallel
# from functools import partial
# import mili_benchmark.scripts.python.milipede_red_help as parallelMOTIF

# pool = mp.Pool(mp.cpu_count())

def parallelMOTIF(index,traces):
    # print(i*chunk_size)
    tf=traces[index]
    TF=os.path.basename(tf)
    data=pd.read_csv(tf,sep='\t',header=None)
    data1=data.groupby([0,1,2,10]).max()
    data1[3]=data1[3].round(4)
    data1[[3,4,8,9]].to_csv('../../d/tmp/redmo/bench/alt/'+TF+'_red.txt',sep='\t',index=True,header=False,mode='a')


traces=glob.glob('../../d/tmp/redmo/bench/alt/*')
indices = [i for i, s in enumerate(traces) if 'motif' not in s]

# results = [pool.apply(parallelMOTIF, args=[index,traces]) for index in indices]
# pool.close() 

for index in indices:
    parallelMOTIF(index,traces)