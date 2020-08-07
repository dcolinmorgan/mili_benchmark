import glob, os
import numpy as np
from matplotlib import colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path
import warnings
from sklearn import metrics
from collections import Counter

warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')
# import matplotlib.backends.backend_pdf
table=[]
tab=[]

def predRegion(indir='data/MotifPipeline/sthlm_motif1010',outdir='data/MotifPipeline/test/sthlm_motif1010/',cell=None,TF=N
one):
    traces= glob.glob(indir+'/*')
    traces = list(filter(lambda file: os.stat(file).st_size > 0, traces))
    WWW=3
    if traces is None:
        print('incompatible intersection directory')
        return
    
    Path(outdir).mkdir(parents=True, exist_ok=True)
    # pdf = matplotlib.backends.backend_pdf.PdfPages("output.pdf")
    # X=np.random.randint(0,high=len(traces))

    for jac,trace in enumerate(traces):
        if (os.path.isfile(trace) & os.stat(trace).st_size != 0):

        #     try:
        # trace='data/MotifPipeline/sthlm_motif_5_QCbeta/red/A549_ATF3'
            data=pd.read_csv(trace,sep='\t',usecols=[0,1,2,3,4,8,9,10,11,15,16],names=["chr", "start", "end",'weight','pval',"hits1",'W1','array','region','gene','ChIPTF'])
            Col1=os.path.basename(trace).split('_')[0] #cell
            Col2=os.path.basename(trace).split('_')[1] #TF
            #         except:
            #             pass
            data.ChIPTF=data.ChIPTF.replace('.',0)
            data.ChIPTF[(data.gene==Col2)]=1
            data=data[(data.ChIPTF==0)|(data.ChIPTF==1)]
            data.ChIPTF=pd.to_numeric(data.ChIPTF)
            if np.sum(data.ChIPTF)>5:
                try:
            #         data=pd.read_csv(trace,sep='\t',usecols=[0,1,2,3,4,8,9,10,11,15,16],names=["chr", "start", "end",'weight','pval',"hits1",'W1','array','region','gene','ChIPTF'])
            #         Col1=os.path.basename(trace).split('_')[0] #cell
            #         Col2=os.path.basename(trace).split('_')[1] #TF


                    data.weight=(data.weight-data.weight.min())/(data.weight.max()-data.weight.min())
                    data['wgbs']=data.W1/100
                    data.wgbs=1-data.wgbs
                    data.array=1-data.array
                    data['cell']=Col1#os.path.basename(trace).split('_')[0]
                    data['TF']=Col2
                    aaa=(data[(data.ChIPTF-data.wgbs)>=-.5]['region'])#FP
                    bbb=(data[(data.ChIPTF-data.wgbs)<-.5]['region'])
                    ccc=(data[(data.ChIPTF-data.wgbs)>=.5]['region'])#FN
                    ddd=(data[(data.ChIPTF-data.wgbs)<.5]['region'])

                    # from collections import Counter
                    aa = dict(Counter(aaa))
                    bb = dict(Counter(bbb))
                    cc = dict(Counter(ccc))
                    dd = dict(Counter(ddd))

                #         column = Col1, Col2, Col3, Col4, Col10, Col5, Col6, Col11,Col7,Col8,Col12,Col9,Col13,Col14,Col15,Col16
                    column = data.cell[1],data.TF[1],aa['N_Shore'],aa['S_Shore'],aa['OpenSea'],aa['N_Shelf'],aa['Island'],bb['N_Shore'],bb['S_Shore'],bb['OpenSea'],bb['N_Shelf'],bb['Island'],
                    cc['N_Shore'],cc['S_Shore'],cc['OpenSea'],cc['N_Shelf'],cc['Island'],dd['N_Shore'],dd['S_Shore'],dd['OpenSea'],dd['N_Shelf'],dd['Island']
                #             table.append(column)
                    print("calculating chromosomal region stats for "+Col2+' in '+Col1)

                    np.transpose(pd.DataFrame((column))).to_csv(outdir+'/sthlm_regionPRE_indv.txt',mode='a',header=None)#,header=['cell','tf','mauroc','wauroc','meauroc','maupr','waupr','meaupr','size','max_hit','mo_length'])
                except:
                    pass

