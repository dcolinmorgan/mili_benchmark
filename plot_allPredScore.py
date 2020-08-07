import re, glob, os, collections, scipy
import numpy as np
# import dask
from matplotlib import colors as mcolors
from IPython.display import Image
# import networkx as nx
import matplotlib.pyplot as plt
from datetime import datetime, date
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from sklearn import metrics
from pathlib import Path
# import sklearn.svm as svm
from scipy import stats

def plot_allPredScore(directory,agg_meth='mean',method='auroc',region=None,depth=None):
	allbox = pd.DataFrame()
	traces= glob.glob('data/MotifPipeline/**/red/test/'+region+'/'+depth+'/*.txt',recursive = True)
	indices = [i for i, s in enumerate(traces) if 'sthlm_'+method+'_meltbox' in s]


	for jac in (indices):
	    trace=traces[jac]
	    buffer=(trace).split('/')[2]
	    buffer=(buffer).split('_')[2]
	    meltbox=pd.read_csv(trace,sep='\t')
	    meltbox['buffer']=buffer
	    allbox=pd.concat([allbox,meltbox],axis=0)
	allbox['cell_buff']=allbox['data - cell line']+'_'+allbox['buffer'].astype(str)
	from pathlib import Path
	outdir='data/MotifPipeline/compare'
	Path(outdir).mkdir(parents=True, exist_ok=True)
	# plt.figure(figsize=(12, 5))
	plt.xticks(rotation=30)

	cells=['A549','GM12878', 'HeLa', 'HepG2', 'K562','SKNSH']
	tests=['pwm','me','wg']
	for cell in cells:
	    for test in tests:
	        allbox2=allbox[allbox['cell_buff'].str.contains(pat=cell)]
	        allbox2=allbox2[allbox2['cell_buff'].str.contains(pat=test)]
	        allbox2.buffer=(allbox2.buffer).astype(int)
	        allbox2=allbox2.sort_values(by='buffer')
	        g=sns.boxplot(x='cell_buff',y=method.upper(), data=allbox2)
	        g=sns.swarmplot(x='cell_buff',y=method.upper(), data=allbox2,size=2, color=".3", linewidth=0)
	        g.set(ylim=(0, 1))
	        g.set_xticklabels(g.get_xticklabels(), rotation=45, horizontalalignment='right')


	        plt.savefig(outdir+"/sthlm_"+cell+test+"_allbox_buff.png",dpi=300,bbox_inches = "tight")
	        plt.close()
	return