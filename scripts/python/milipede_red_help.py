def parallelMOTIF(index,traces):
    # print(i*chunk_size)
    tf=traces[index]
    TF=os.path.basename(tf)
    data=pd.read_csv(tf,sep='\t',header=None)
    data1=data.groupby([0,1,2,10]).max()
    data1[3]=data1[3].round(4)
    data1[[3,4,8,9]].to_csv('../../d/tmp/redmo/bench/alt/'+TF+'_red.txt',sep='\t',index=True,header=False,mode='a')

