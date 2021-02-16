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
# import sklearn.svm as svm
from scipy import stats

def plot_predScore(indir,outdir):

	# table2=table
	# table.to_csv('data/MotifPipeline/PRE_overall.txt')
	distances=pd.read_table(indir+agg_meth+'/sthlm_PRE_overall.txt',sep=',',usecols=[1,2,9,10,11,12,13,14,15,16,17],names=['cell','TF','size','counts','length',
		'mean_dist_array','mean_dist_wg','mean_dist_pwm','hits_array','hits_wg','hits_pwm'])
	distances=distances[distances.cell!='0']

	# if measure=='auroc':
	# 	aurocs.drop(['mo_aupr','wgbs_aupr','me_aupr','rand_mo_aupr','rand_wgbs_aupr','rand_me_aupr','tr_wg_aupr','tr_me_aupr'],axis=1,inplace=True)
	# elif measure=='aupr':
	# 	aurocs.drop(['mo_auroc','wgbs_auroc','me_auroc','rand_mo_auroc','rand_wgbs_auroc','rand_me_auroc','tr_wg_auroc','tr_me_auroc'],axis=1,inplace=True)


	###PLOT
	# measure="wg_auroc"
	# measure='auroc'

	##sum across cells and sort
	heat=aurocs.pivot_table(index=['TF'], columns='cell')
	# heat=heat[[measure]]
	heat['mean']=np.nanmean(heat,axis=1)
	heat['count']=heat.isnull().sum(axis=1)
	heat['count']=(21-heat['count'])
	# heat=heat.sort_values(by=['count'],ascending=True)
	# heat=heat.sort_values(by=['mean'],ascending=False)
	heat['weight']=heat['mean']*(heat['count'])
	heat=heat.sort_values(by=['weight'],ascending=False)
	# heat=heat[heat.index.str.contains(pat='CEB')]

	# heat=heat[heat['count']==5]
	# heat=heat[heat['mean']<0.1]
	# del heat.iloc[:,0]
	heat=heat.dropna(axis=1, how='all')
	del heat['count']
	del heat['mean']
	del heat['weight']
	# heat['sum']=heat.isnull().sum(axis=0)
	# heat=heat.sort_values(by=['sum'],ascending=True)
	# del heat['sum']
	heat=pd.DataFrame(heat.to_records())

	heat.columns=['TF','me-A549','me-GM12878','me-HeLa','me-HepG2','me-K562','me-SKNSH',
	             'pwm-A549','pwm-GM12878','pwm-HeLa','pwm-HepG2','pwm-K562','pwm-SKNSH',
	             'wg-A549','wg-GM12878','wg-HeLa','wg-HepG2','wg-K562','wg-SKNSH']

	# heat=heat[['TF','pwm-A549','pwm-GM12878','pwm-HeLa','pwm-HepG2','pwm-K562','pwm-SKNSH',
	#           'me-A549','me-GM12878','me-HeLa','me-HepG2','me-K562','me-SKNSH',
	#              'wg-A549','wg-GM12878','wg-HeLa','wg-HepG2','wg-K562','wg-SKNSH']]


	heat=heat[['TF','pwm-A549','me-A549','wg-A549',
	           'pwm-GM12878','me-GM12878','wg-GM12878',
	           'pwm-HeLa', 'me-HeLa','wg-HeLa',
	           'pwm-HepG2','me-HepG2','wg-HepG2',
	           'pwm-K562','me-K562','wg-K562',
	           'pwm-SKNSH','me-SKNSH','wg-SKNSH']]
	heat55=heat
	heat=heat.set_index('TF')

	plt.figure(figsize=(12, 30))
	# grid_kws = {"height_ratios": (60,12), "hspace": .3}

	# f, (ax, cbar_ax) = plt.subplots(2)#, gridspec_kw=grid_kws)
	# rdgn = sns.diverging_palette(h_neg=200, h_pos=20, s=99, l=55, sep=3, as_cmap=True)
	# rdgn=sns.diverging_palette(220, 20, sep=20, as_cmap=True)
	rdgn=sns.color_palette("RdBu", 100)
	ax = sns.heatmap(heat, center=.5,linewidth=.1,cmap=rdgn, cbar=False)#cbar_ax=cbar_ax,cbar_kws={"orientation": "horizontal"})
	mappable = ax.get_children()[0]
	plt.colorbar(mappable, ax = [ax],orientation = 'horizontal',pad=.04) #.02 with 12x60, .03 for 12x40
	plt.xticks(rotation=30)
	# plt.imshow(heat, cmap='hot', interpolation='nearest')
	# ax= sns.heatmap(heat, center=.5,cmap=rdgn,linewidth=.1)
	# ax= sns.heatmap(heat, linewidth=.1,cmap="YlGnBu")

	plt.savefig(outdir+"sthlm_"+measure+"_heat.png",dpi=300,bbox_inches = "tight")
	plt.show

	print([(heat55['pwm-A549'].dropna().shape),(heat55['pwm-GM12878'].dropna().shape),(heat55['pwm-HeLa'].dropna().shape),(heat55['pwm-HepG2'].dropna().shape),(heat55['pwm-K562'].dropna().shape),(heat55['pwm-SKNSH'].dropna().shape)])

	box=heat55
	# box.columns=['TF','A549','GM12878','H1','HeLa','HepG2','K562','SKNSH']
	meltbox=pd.melt(box,id_vars=['TF'])
	# del box.TF
	if measure=='auroc':
		meltbox.columns=['TF','data - cell line','AUROC']
		plt.figure(figsize=(12, 5))
		plt.xticks(rotation=30)
		g=sns.boxplot(x='data - cell line',y='AUROC', data=meltbox)
		g=sns.swarmplot(x='data - cell line',y='AUROC', data=meltbox,
		              size=2, color=".3", linewidth=0)
	elif measure=='aupr':
		meltbox.columns=['TF','data - cell line','AUPR']
		plt.figure(figsize=(12, 5))
		plt.xticks(rotation=30)
		g=sns.boxplot(x='data - cell line',y='AUPR', data=meltbox)
		g=sns.swarmplot(x='data - cell line',y='AUPR', data=meltbox,
		              size=2, color=".3", linewidth=0)
	g.set(ylim=(0, 1))
	plt.savefig(outdir+"sthlm_"+measure+"_box.png",dpi=300,bbox_inches = "tight")
	# meltbox.to_csv(outdir+measure+"meltbox.txt")

	plt.show
	# Tweak the visual presentation
	# ax.xaxis.grid(True)
	# ax.set(ylabel="")
	sns.despine(trim=True, left=True)
	# plt.show()


	t0a, p0a = stats.ttest_ind(box['pwm-A549'].dropna(),box['wg-A549'].dropna())
	t0b, p0b = stats.ttest_ind(box['pwm-A549'].dropna(),box['me-A549'].dropna())
	t1a, p1a = stats.ttest_ind(box['pwm-GM12878'].dropna(),box['wg-GM12878'].dropna())
	t1b, p1b = stats.ttest_ind(box['pwm-GM12878'].dropna(),box['me-GM12878'].dropna())
	t2a, p2a = stats.ttest_ind(box['pwm-HeLa'].dropna(),box['wg-HeLa'].dropna())
	t2b, p2b = stats.ttest_ind(box['pwm-HeLa'].dropna(),box['me-HeLa'].dropna())
	t3a, p3a = stats.ttest_ind(box['pwm-HepG2'].dropna(),box['wg-HepG2'].dropna())
	t3b, p3b = stats.ttest_ind(box['pwm-HepG2'].dropna(),box['me-HepG2'].dropna())
	t4a, p4a = stats.ttest_ind(box['pwm-K562'].dropna(),box['wg-K562'].dropna())
	t4b, p4b = stats.ttest_ind(box['pwm-K562'].dropna(),box['me-K562'].dropna())
	t5a, p5a = stats.ttest_ind(box['pwm-SKNSH'].dropna(),box['wg-SKNSH'].dropna())
	t5b, p5b = stats.ttest_ind(box['pwm-SKNSH'].dropna(),box['me-SKNSH'].dropna())


	# initialise data of lists. 
	ttest = {'WG ttest':[t0a,t1a,t2a,t3a,t4a,t5a], 'WG pvalue':[p0a,p1a,p2a,p3a,p4a,p5a], 'Me ttest':[t0b,t1b,t2b,t3b,t4b,t5b],'Me pvalue':[p0b,p1b,p2b,p3b,p4b,p5b]} 
	  
	# Creates pandas DataFrame. 
	df_ttest = pd.DataFrame(ttest, index =['A549','GM12878', 'HeLa', 'HepG2', 'K562','SKNSH']) 

	df_ttest.to_csv(outdir+agg_meth+"/sthlm_"+measure+"_ttest.txt")
	meltbox.to_csv(outdir+agg_meth+"/sthlm_"+measure+"_meltbox.txt",sep='\t')

	return




