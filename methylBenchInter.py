import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
from matplotlib_venn import venn3
from matplotlib import pyplot as plt

def methylBenchInter(cell):

	wgbs=pd.read_csv('data/MotifPipeline/ENCODE/wgbsin/'+cell+'both.txt',usecols=[0,1,2],sep='\t',names=['chr','start','end'])
	if cell=='A549':
		cellline='A-549'
	elif cell=='HeLa':
		cellline='HeLa-S3'
	elif cell=='HepG2':
		cellline='Hep-G2'
	elif cell=='K562':
		cellline='K-562'
	elif cell=='SKNSH':
		cellline='SK-N-SH'
	else:
		cellline=cell
	chip=pd.read_csv('data/MotifPipeline/remap/'+cellline+'.bed',usecols=[0,1,2],sep='\t',names=['chr','start','end'])
	me=pd.read_csv('data/MotifPipeline/ENCODE/methyl_array/crossQC/'+cellline+'_shufMeArrayHG38b.txt',usecols=[0,1,2],sep='\t',names=['chr','start','end'])

	mechip = pd.merge(me, chip, left_on=['chr','start','end'], right_on = ['chr','start','end'])
	# del mechip['type_x'], mechip['type_y']
	mewgbs = pd.merge(me, wgbs, left_on=['chr','start','end'], right_on = ['chr','start','end'])
	# del mewgbs['type_x'], mewgbs['type_y']
	chipwgbs = pd.merge(chip, wgbs, left_on=['chr','start','end'], right_on = ['chr','start','end'])
	# del chipwgbs['type_x'], chipwgbs['type_y']
	all = pd.merge(mewgbs, chip, left_on=['chr','start','end'], right_on = ['chr','start','end'])
	# del all['type_x'], all['type_y']

	ual=len(all.drop_duplicates())
	umc=len(mechip.drop_duplicates())-ual
	umw=len(mewgbs.drop_duplicates())-ual
	ucw=len(chipwgbs.drop_duplicates())-ual

	meme=len(me.drop_duplicates())-umc-umw-ual
	chch=len(chip.drop_duplicates())-umc-ucw-ual
	wgwg=len(wgbs.drop_duplicates())-ucw-umw-ual



	plt.figure(figsize=(5,5))
	v3 = venn3(subsets = {'100':np.log(meme),
	                      '010':np.log(chch),
	                      '110':np.log(umc),
	                      '001':np.log(wgwg),
	                      '101':np.log(umw),
	                      '011':np.log(ucw),
	                      '111':np.log(ual)},
	                    set_labels = ('', '', ''))

	# plt.figure( figsize =(7,7) )
	v3.get_patch_by_id('100').set_color('tab:orange')
	v3.get_patch_by_id('010').set_color('deeppink')
	v3.get_patch_by_id('001').set_color('tab:blue')
	# try:
	# v3.get_patch_by_id('110').set_color('orange')
	# # except AttributeError:
	# #     pass
	# v3.get_patch_by_id('101').set_color('purple')
	# v3.get_patch_by_id('011').set_color('green')
	# v3.get_patch_by_id('111').set_color('grey')

	v3.get_label_by_id('100').set_text('B. meArray\n'+str(meme))
	v3.get_label_by_id('010').set_text('D. ChIP\n'+str(chch))
	v3.get_label_by_id('001').set_text('A. WGBS\n'+str(wgwg))
	v3.get_label_by_id('110').set_text(str(umc))
	v3.get_label_by_id('101').set_text(str(umw))
	v3.get_label_by_id('011').set_text(str(ucw))
	v3.get_label_by_id('111').set_text(str(ual))
	# plt.title("Overlap Venn diagram (areas log scaled)")
	# plt.tight_layout()
	# plt.figure(figsize=(18,18))
	# plt.title(cell)
	for text in v3.subset_labels:
	    text.set_fontsize(15)

	plt.savefig(cell+'_venn.png')

	return v3# plt.show()

