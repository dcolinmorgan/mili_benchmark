from __future__ import print_function

import sys, os, glob,re 
# sys.path.insert(1,'../panda')
from netZooPy.panda.panda import Panda
import numpy as np
import collections
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import netZooPy
from datetime import datetime, date
timestr = time.strftime("%Y%m%d-%H%M")

class null_panda(Panda):
    '''Null model to compare meaningful panda models to.'''
  def __init__(self,data,design,inits,computation='gpu'):
    
  FC=data.iloc[:,data.columns.str.contains('0_1')]
  FNC=data.iloc[:,data.columns.str.contains('0_0')]
  MC=data.iloc[:,data.columns.str.contains('1_1')]
  MNC=data.iloc[:,data.columns.str.contains('1_0')]

  for i in range(1,inits):
    ## four random group exact same size as initial cases and controls
    j=np.random.choice(data.shape[1], size=(MC.shape[1]+MNC.shape[1], ),replace=False) ## total males
    nullA=data.iloc[:,j.tolist()] ##reduce to random set
    null2=data.drop(nullA.columns,axis=1) ## create second group without these

    k=np.random.choice(j,MNC.shape[1],replace=False) ## male controls
    nullB=data.iloc[:,k.tolist()] ##remove controls from male pop
    nullA=nullA.drop(nullB.columns,axis=1) ## remove these from male cases

    l=np.random.choice(null2.shape[1], size=(FC.shape[1]+FNC.shape[1], ),replace=False) ##84 total females
    nullC=null2.iloc[:,l.tolist()] ##subset leftover from first null to 85
    m=np.random.choice(l[l<FC.shape[1]+1+FNC.shape[1]],FNC.shape[1],replace=False) ##25 controls
    nullD=nullC.iloc[:,m.tolist()] ##take 25 for female random control
    nullC=nullC.drop(nullD.columns,axis=1)

    negCont(nullA,nullB,i)
    negCont(nullC,nullD,i+1)

    print('saved iteration='+str(i))

  def negCont(nullA,nullB,i):
    panda_nullA = Panda(nullA, motif_data, ppi_data,computing='gpu',precision='single',save_tmp=False,save_memory = True, remove_missing=False, keep_expression_matrix = False,modeProcess = 'intersection')
    # np.save('drive/My Drive/Panda_nullA.npy',panda_nullA.export_panda_results)
    panda_nullB = Panda(nullB, motif_data, ppi_data,computing='gpu',precision='single',save_tmp=False,save_memory = True, remove_missing=False, keep_expression_matrix = False,modeProcess = 'intersection')
    # np.save('drive/My Drive/Panda_nullB.npy',panda_nullB.export_panda_results)

    subset_panda_resultsA = panda_nullA.export_panda_results.sort_values(by=['force'], ascending=False)
    subset_panda_resultsB = panda_nullB.export_panda_results.sort_values(by=['force'], ascending=False)
    del subset_panda_resultsA['motif'],subset_panda_resultsB['motif']
    del panda_nullA, panda_nullB

    ##aggregate TF link sum and zscore
    TFsum_A=subset_panda_resultsA.groupby(['tf']).agg({'force':sum})
    TFsum_B=subset_panda_resultsB.groupby(['tf']).agg({'force':sum})
    TFsum_B['force']=zscore(TFsum_B['force'])
    TFsum_A['force']=zscore(TFsum_A['force'])

    ##aggregate gene link sum and zscore
    geneSum_A=subset_panda_resultsA.groupby(['gene']).agg({'force':sum})
    geneSum_B=subset_panda_resultsB.groupby(['gene']).agg({'force':sum})
    geneSum_A['force']=zscore(geneSum_A['force'])
    geneSum_B['force']=zscore(geneSum_B['force'])

    ##merge
    sumTFdiff=TFsum_A.merge(TFsum_B, left_on=['tf'], right_on=['tf'])
    sumGenediff=geneSum_A.merge(geneSum_B, left_on=['gene'], right_on=['gene'])

    del [[subset_panda_resultsA, subset_panda_resultsB]]
    gc.collect()
    subset_panda_resultsA=pd.DataFrame()
    subset_panda_resultsB=pd.DataFrame()

    ##diff
    sumTFdiff['diff']=sumTFdiff.force_x-sumTFdiff.force_y
    # sumTFdiff = sumTFdiff.sort_values(by ='diff' )
    del sumTFdiff['force_x'], sumTFdiff['force_y']
    sumTFdiff=sumTFdiff.sort_index(axis=0)
    sumTFdiff=sumTFdiff.T
    sumGenediff['diff']=sumGenediff.force_x-sumGenediff.force_y
    # sumGenediff = sumGenediff.sort_values(by ='diff' )
    del sumGenediff['force_x'], sumGenediff['force_y']
    sumGenediff=sumGenediff.sort_index(axis=0)
    sumGenediff=sumGenediff.T
    if i==1:
      sumTFdiff.to_csv('dTF_zscore_rand_'+timestr+'.txt',sep='\t',header=True,index=None)
      sumGenediff.to_csv('gene_zscore_rand_'+timestr+'.txt',sep='\t',header=True,index=None)
    elif i!=1:
      sumTFdiff.to_csv('TF_zscore_rand_'+timestr+'.txt',sep='\t',header=None,index=None,mode='a')
      sumGenediff.to_csv('gene_zscore_rand_'+timestr+'.txt',sep='\t',header=None,index=None,mode='a')
    
    del sumTFdiff, sumGenediff
    gc.collect()

