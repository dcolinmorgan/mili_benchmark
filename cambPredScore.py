
import glob, os
import numpy as np
from matplotlib import colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path
import warnings
from sklearn import metrics



warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')
# import matplotlib.backends.backend_pdf

    
def cambPredScore(indir='../../pc/redmo/data/MotifPipeline/camb_motif_pipeline_gamma5'): #,outdir='../../pc/redmo/data/MotifPipeline/camb_motif_pipeline_gamma5'):

    traces = [log for log in glob.glob(indir+'/*') if not os.path.isdir(log)]
    aurocs = pd.DataFrame()
    os.makedirs(indir+'/test')

    for i,trace in enumerate(traces):
        # print(trace)
    #     plt.subplot(2, 2, 3)
        combo=trace.split('/')[7]
        TF=combo.split('_')[1]
        cell=combo.split('_')[0]
    #     data=pd.read_csv(trace,sep='\t',names=['chr','start','end','pwm','CG','ChIPTF','array','location'])
        data=pd.read_csv(trace,sep='\t',usecols=[0,1,2,3,9,10,11],names=["chr", "start", "end",'pwm',"W1",'array','ChIPTF'])

        Col1=os.path.basename(trace).split('_')[0] #cell
        Col2=os.path.basename(trace).split('_')[1] #TF
    # except:
    #     pass

        table2=[]
        tbl=[]
        tmpTBL2=[]
        tmpTBL=[]
        table3=[]
        tblC=[]
        tmpTBLC=[]
        data.pwm=(data.pwm-data.pwm.min())/(data.pwm.max()-data.pwm.min())
        data=data.fillna(0)
        data.ChIPTF=data.ChIPTF.replace('.',0)
        # print(data.shape)
        data.ChIPTF[(data.ChIPTF==Col2)]=1
        data=data[(data.ChIPTF==0)|(data.ChIPTF==1)]
        data.ChIPTF=pd.to_numeric(data.ChIPTF)

        if np.sum(data.ChIPTF)>5: # and np.sum(data_nonCG.ChIPTF)>5:

            plt.plot([0, 1], [0, 1], 'k--')
            fpr, tpr, thresholds = metrics.roc_curve(data.ChIPTF, data.pwm)
            roc_auc=metrics.auc(fpr, tpr)
            plt.plot(fpr, tpr,
                     label='pwm (area = {0:0.2f})'
                         ''.format(roc_auc),
                     color='b', linestyle=':', linewidth=4)

            fpr2, tpr2, thresholds = metrics.roc_curve(data.ChIPTF, 1-(data.W1/100))
            roc_auc2=metrics.auc(fpr2, tpr2)
            plt.plot(fpr2, tpr2,
                     label='wgbs (area = {0:0.2f})'
                           ''.format(roc_auc2),
                     color='r', linestyle=':', linewidth=4)

            fpr, tpr, thresholds = metrics.roc_curve(data.ChIPTF, 1-(data.array/100))
            roc_auc3=metrics.auc(fpr, tpr)
            plt.plot(fpr, tpr,
                     label='array (area = {0:0.2f})'
                         ''.format(roc_auc),
                     color='g', linestyle=':', linewidth=4)



            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title('AUROC for '+TF+' in '+cell)
            plt.legend(loc="best")

            plt.show()

            Col1=cell
            Col2=TF
            Col3=roc_auc #motif auroc
            Col4=roc_auc2
            Col5=roc_auc3
            column = [Col1, Col2, Col3, Col4,Col5]
            column=np.transpose(pd.DataFrame((column)))
            column.to_csv(indir+'/test/camb_PRE_window_methyl.txt',mode='a',header=False,index=False)
            aurocs=pd.concat([aurocs,column],axis=0)
            print("AUROC calculated for "+Col1+" in "+Col2)

    aurocs.columns=['cell','TF','pwm','wgbs','array']
    ##sum across cells and sort
    heat=aurocs.pivot_table(index=['TF'], columns='cell',aggfunc='first')# heat=heat[[measure]]
    heat['mean']=np.nanmean(heat,axis=1)
    heat['count']=heat.isnull().sum(axis=1)

    heat['weight']=heat['mean']*(heat['count'])
    heat=heat.sort_values(by=['weight'],ascending=False)

    heat=heat.dropna(axis=1, how='all')
    del heat['count']
    del heat['mean']
    del heat['weight']

    heat=pd.DataFrame(heat.to_records())

    heat.columns=['TF','me-A549','me-GM12878','me-HeLa','me-HepG2','me-K562','me-SKNSH',
                  'pwm-A549','pwm-GM12878','pwm-HeLa','pwm-HepG2','pwm-K562','pwm-SKNSH',
                'wg-A549','wg-GM12878','wg-HeLa','wg-HepG2','wg-K562','wg-SKNSH']

    heat=heat[['TF','pwm-A549','me-A549','wg-A549',
               'pwm-GM12878','me-GM12878','wg-GM12878',
               'pwm-HeLa', 'me-HeLa','wg-HeLa',
               'pwm-HepG2','me-HepG2','wg-HepG2',
               'pwm-K562','me-K562','wg-K562',
               'pwm-SKNSH','me-SKNSH','wg-SKNSH']]
    heat55=heat
    heat=heat.set_index('TF')

    plt.figure(figsize=(12, 30))
    rdgn=sns.color_palette("RdBu", 100)
    ax = sns.heatmap(heat, center=.5,linewidth=.1,cmap=rdgn, cbar=False)#cbar_ax=cbar_ax,cbar_kws={"orientation": "horizontal"})
    mappable = ax.get_children()[0]
    plt.colorbar(mappable, ax = [ax],orientation = 'horizontal',pad=.04) #.02 with 12x60, .03 for 12x40
    plt.xticks(rotation=30)

    plt.savefig(indir+"/test/camb_auroc_heat.png",dpi=300,bbox_inches = "tight")
    plt.show

    print([(heat55['pwm-A549'].dropna().shape),(heat55['pwm-GM12878'].dropna().shape),(heat55['pwm-HeLa'].dropna().shape),(heat55['pwm-HepG2'].dropna().shape),(heat55['pwm-K562'].dropna().shape),(heat55['pwm-SKNSH'].dropna().shape)])

    box=heat55
    # box.columns=['TF','A549','GM12878','H1','HeLa','HepG2','K562','SKNSH']
    meltbox=pd.melt(box,id_vars=['TF'])
    # del box.TF
    # if measure=='auroc':
    meltbox.columns=['TF','data - cell line','AUROC']
    plt.figure(figsize=(12, 5))
    plt.xticks(rotation=30)
    g=sns.boxplot(x='data - cell line',y='AUROC', data=meltbox)
    g=sns.swarmplot(x='data - cell line',y='AUROC', data=meltbox,
                  size=2, color=".3", linewidth=0)

    g.set(ylim=(0, 1))
    plt.savefig(indir+"/test/camb_auroc_box.png",dpi=300,bbox_inches = "tight")
    # meltbox.to_csv(outdir+measure+"meltbox.txt")

    plt.show

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

    df_ttest.to_csv(indir+"/test/camb_auroc_ttest.txt")
    meltbox.to_csv(indir+"/test/camb_auroc_meltbox.txt",sep='\t')



# def camb_allPredScore(,outdir='../../pc/redmo/data/MotifPipeline/motif_window_analysis/',method='auroc'):
    indir='../../pc/redmo/data/MotifPipeline/**'
    outdir='../../pc/redmo/data/MotifPipeline'
    allbox = pd.DataFrame()
    traces= glob.glob(indir+'/tmp/*.txt',recursive = True)
    indices = [i for i, s in enumerate(traces) if 'camb_auroc_meltbox' in s]


    for jac in (indices):
        trace=traces[jac]
        buffer=(trace).split('/')[2]
        buffer=(buffer).split('_')[2]
        meltbox=pd.read_csv(trace,sep='\t')
        meltbox['buffer']=buffer
        allbox=pd.concat([allbox,meltbox],axis=0)
    allbox['cell_buff']=allbox['data - cell line']+'_'+allbox['buffer'].astype(str)
    from pathlib import Path
    # outdir='data/MotifPipeline/compare'
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


            plt.savefig(outdir+"camb_allbox_buff.png",dpi=300,bbox_inches = "tight")
            plt.close()
    return
