#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import os, glob
import pandas as pd
import sys

# indir = sys.argv[1] 
# outdir = sys.argv[2] 
# import multiprocessing as mp# to run FIMO in parallel
# from functools import partial
# import mili_benchmark.scripts.python.milipede_red_help as parallelMOTIF

"""
Example:
  source /proj/relibs/relib00/conda/bin/activate
  source activate mypy3
  python mili_benchmark/scripts/python/milipede_red.py -i ../../d/tmp/redmo/bench/250 -o test
"""

def main(argv):
    #Create variables
    indir = None
    try:
        opts, args = getopt.getopt(argv, 'hi:o:', ['indir=','outdir='])
    except getopt.GetoptError as err:
        print(str(err))  # will print something like "option -a not recognized"
        print(__doc__)
        return 2

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(__doc__)
            return 0
        elif opt in ('-i', '--indir'):
            indir = arg
        elif opt in ('-o', '--outdir'):
            outdir = arg
        else:
            print('Unknown option', opt)
            return 1


# pool = mp.Pool(mp.cpu_count())

def milipede_red(index,indir,outdir):
    # print(i*chunk_size)

    tf=traces[index]
    TF=os.path.basename(tf)
    TF=TF.split('.')[0]
    data=pd.read_csv(tf,sep='\t',usecols=[0,1,2,3,4,8,12,13,14],
                     names=["chr", "start", "end",'weight','TF','gene','depth','W1','ChIPTF'])
    data=data[data['W1']!='.']
    data['W1']=pd.to_numeric(data['W1'])
    data['W1']=1-(data.W1/100)

    # data1=data.groupby(['chr','start','end','gene']).agg({'weight':'max',"W1":'max',"ChIPTF":'max'})

    data.ChIPTF[data.ChIPTF!='0']=1
    data['ChIPTF']=pd.to_numeric(data['ChIPTF'])
    if len(data)>0:
        data1=data.groupby(['chr','start','end','TF','gene']).agg({'weight':'mean',"W1":'mean',"ChIPTF":'max'})
    
        fpr, tpr, thresholds = metrics.roc_curve(data1.ChIPTF, data1.weight)
        roc_auc=metrics.auc(fpr, tpr)
        fpr, tpr, thresholds = metrics.roc_curve(data1.ChIPTF,data1.W1)
        roc_auc2=metrics.auc(fpr, tpr)
        data1.to_csv(indir+'/'+outdir+'/alt_motif.txt',header=True,sep='\t',mode='a')
        column = TF, roc_auc, roc_auc2,len(data),len(data1)

        np.transpose(pd.DataFrame((column))).to_csv(indir+'/'+outdir+'/alt_results.txt',mode='a',header=['tf','pwm_auroc','W1_auroc','full_len','agg_len'])

def milipede_red_plot(inputA,outdir):
    meltboxA=pd.read_csv(inputA+'/'+outdir+'/alt_results.txt',sep=',')
    meltboxA=meltboxA.dropna(how='any')

    total_reads=pd.read_csv('data/MotifPipeline/compare/total_reads.txt',sep=' ',names=['count','motif'])
    TFs=pd.read_csv('mili_benchmark/data/Homo_sapiens_motifinfo.txt',sep='\t',header=None,names=['motif','tf'])
    total2=total_reads.merge(TFs)

    meltboxA=meltboxA.merge(total2)
    meltboxA['perc_full']=np.divide(meltboxA['agg_len'].values.astype('float'),
              meltboxA['count'].values.astype('float'))
    meltboxA['bin']=pd.cut(meltboxA['perc_full'], 6)

    goldboxA=pd.melt(meltboxA,id_vars=['tf','bin'],value_vars=['pwm_auroc','W1_auroc'])

    g=plt.figure(figsize=(8, 5))
    # g=sns.boxplot(x='variable',y=goldboxA.value.astype('float'),hue='bin', data=goldboxA)
    g=sns.violin(x='variable',y=goldboxA.value.astype('float'),hue='bin',
             split=False,scale="count",inner="quartile",data=goldboxA)
    g.set(ylim=(0, 1))
    g.set(ylim=(0, 1))
    g.set_xticklabels(g.get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.savefig(inputA+'/'+outdir+'/boxplot.png')

traces=glob.glob(indir)
indices = [i for i, s in enumerate(traces) if '_' not in s]
for index in indices:
	milipede_red(index,indir,outdir)

milipede_red_plot(indir,outdir)


# results = [pool.apply(parallelMOTIF, args=[index,traces]) for index in indices]
# pool.close() 





