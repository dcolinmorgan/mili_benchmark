
import glob, os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path
import warnings
from sklearn import metrics
import scipy as sp

warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')
# import matplotlib.backends.backend_pdf

# agg_meth='median'
# buffer='500' -- inherited from outdir

def buffer_distr_comp(agg_meth='mean',outdir='data/MotifPipeline/sthlm_motif_50_QCbeta/red/'):
    i=0
    traces= glob.glob('data/MotifPipeline/sthlm_motif_0_QCbeta/*')
    traces = list(filter(lambda file: os.stat(file).st_size > 0, traces))
    buffer=outdir.split('/')[2]
    buffer=buffer.split('_')[2]
    pwm_sp=[]
    wg_sp=[]
    me_sp=[]
    Path(outdir+'/'+agg_meth).mkdir(parents=True, exist_ok=True)
    # indices = [i for i, s in enumerate(traces) if 'sthlm_auroc_meltbox' in s]
    for trace in (traces):

        try:
            data=pd.read_csv(trace,sep='\t',usecols=[0,1,2,3,4,8,9,10,11,15,16],
                             names=["chr", "start", "end",'weight','pval',"hits1",'W1','array','region','gene','ChIPTF'],engine='python')
            data.weight=(data.weight-data.weight.min())/(data.weight.max()-data.weight.min())
            data['wgbs']=data.W1/100
            data.wgbs=1-data.wgbs
            data.array=1-data.array
            # print(data)
        
        
        
            cell=os.path.basename(trace).split('_')[0]
            TF=os.path.basename(trace).split('_')[1]
            data = data.groupby([data.chr,data.start,data.end]).agg({'weight':agg_meth,"wgbs":agg_meth,"ChIPTF":'max',"array":agg_meth,'hits1':'mean','W1':'mean'})
            
            trace2=os.path.basename(trace)#.split('/')[4]
            com_data=pd.read_csv(outdir+'/'+trace2,sep='\t',usecols=[0,1,2,3,4,8,9,10,11,15,16],
                                 names=["chr", "start", "end",'weight','pval',"hits1",'W1','array','region','gene','ChIPTF'],engine='python')
            com_data.weight=(data.weight-data.weight.min())/(data.weight.max()-data.weight.min())
            com_data['wgbs']=com_data.W1/100
            com_data.wgbs=1-com_data.wgbs
            com_data.array=1-com_data.array
#           print(com_data)

            com_data = com_data.groupby([com_data.chr,com_data.start,com_data.end]).agg({'weight':agg_meth,"wgbs":agg_meth,"ChIPTF":'max',"array":agg_meth,'hits1':'mean','W1':'mean'})

            data2=data.merge(com_data,on=['chr','start','end'])
            pwm=(sp.stats.spearmanr(data2.weight_x,data2.weight_y)[0])
            wg=(sp.stats.spearmanr(data2.wgbs_x,data2.wgbs_y)[0])
            me=(sp.stats.spearmanr(data2.array_x,data2.array_y)[0])
            i=i+1
            column= pwm,wg,me
#             print(pwm_sp,wg_sp,me_sp)
            np.transpose(pd.DataFrame((column))).to_csv(outdir+'/'+agg_meth+'/'+buffer+'_'+agg_meth+'_sthlm_distr_comp.txt',mode='a',index=False,header=False)
#             print('running distribution comparison between 0 and '+buffer+' for '+TF+' in '+cell)
        #         return pwm_sp
            pwm_sp.append(pwm)
            wg_sp.append(wg)
            me_sp.append(me)
            column= pwm_sp,wg_sp,me_sp
        except:
            pass
            print('missing data for taking '+agg_meth+' of '+TF+' in '+cell)
    return column

