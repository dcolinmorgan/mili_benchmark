import matplotlib.pyplot as plt
import numpy as np
import os, glob
import pandas as pd
import multiprocessing as mp# to run FIMO in parallel
from functools import partial
import mili_benchmark.scripts.python.milipede_red_help as parallelMOTIF


def main(argv):
    #Create variables
    indir = None
    try:
        opts, args = getopt.getopt(argv, 'hi:', ['help', 'indir='])
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
        else:
            print('Unknown option', opt)
            return 1


pool = mp.Pool(mp.cpu_count())

def milipede_red(index,indir):
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
        data1.to_csv(indir+'/test/alt_motif.txt',header=True,sep='\t',mode='a')
        column = TF, roc_auc, roc_auc2,len(data),len(data1)

        np.transpose(pd.DataFrame((column))).to_csv(indir+'/test/alt_results.txt',mode='a',header=['tf','pwm_auroc','W1_auroc','full_len','agg_len'])



traces=glob.glob(indir)
indices = [i for i, s in enumerate(traces) if '_red' not in s]
for index in indices:
	milipede_red(index,indir)


# results = [pool.apply(parallelMOTIF, args=[index,traces]) for index in indices]
# pool.close() 





