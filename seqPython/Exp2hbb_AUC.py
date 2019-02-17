# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 20:47:08 2019

@author: 51164
"""

## 多个k值的 以全部数据集求权重和特征
import ReadData
import Sequence
import Similarity
import time
import markov
import numpy as np

from sklearn.metrics import roc_auc_score

## 返回列表的第2个 元素，即相似度
def getSim(Lis):
    return Lis[1]

if __name__=="__main__":
    rd=ReadData.ReadData()
#    datasets,pos,neg=rd.getData2("fly_blastoderm")
#    name="human_muscle"
    name="HBB"

#    name="pns"
    datasets,pos,neg=rd.getData2(name)
    print(name)

    Sim=Similarity.Similarity()
    kstart=2
    kend=6
    sizePo=len(pos)
    top=int(sizePo*(sizePo-1))
    flag=True
#    top=3321
    sq=Sequence.Sequence()
    
    ## 计算权重
    start = time.process_time()
    posweight = sq.getMulWeight_mid(pos,kstart,kend)
    negweight = sq.getMulWeight_mid(neg,kstart,kend)
    end = time.process_time()
    print("正权重个数",len(posweight))
    print("负权重个数",len(negweight))
    print("权重计算时间：",(end-start))

    print("--------------------------MulD2-----------------------")
    print("kstart=",kstart," kend=",kend)
    start = time.process_time()
    ## 计算MulD2特征
    posD2DicLis=sq.getMulCount_suf(pos,kstart,kend,pos) 
    negD2DicLis=sq.getMulCount_suf(neg,kstart,kend,neg)
    
    
    posweightLis=[]
    negweightLis=[]
    posD2Lis=[]
    negD2Lis=[]
    #####将特征与 权重理顺： 正权重
    for key in sorted(posweight):
        posweightLis.append(posweight[key])
        
     #####将特征与 权重理顺： 负权重
    for key in sorted(negweight):
        negweightLis.append(negweight[key])
      
    ## 正特征
    for i in range(len(posD2DicLis)):
        tmp=[]
        for key in sorted(posweight):
            tmp.append(posD2DicLis[i][key])
        posD2Lis.append(tmp)
        
    ## 负特征
    for i in range(len(negD2DicLis)):
        tmp=[]
        for key in sorted(negweight):
            tmp.append(negD2DicLis[i][key])
        negD2Lis.append(tmp)
    
    ## 构成基因对,并设置标识
    posPair=[]
    negPair=[]
    for i in range(len(posD2Lis)):
        for j in range(i+1,len(posD2Lis)):
            posSim=Sim.getMulD2WeightSim3Lis(posD2Lis[i],posD2Lis[j],posweightLis)
            negSim=Sim.getMulD2WeightSim3Lis(negD2Lis[i],negD2Lis[j],negweightLis)
            posPair.append([1,posSim])
            negPair.append([0,negSim])
    Pair=posPair+negPair
    pSim=sorted(Pair,key=getSim,reverse=True)
    print("数据集大小：",len(pSim))
    PP=np.array(pSim)
    label=PP[:,0]
    pred=PP[:,1]    
    auc=roc_auc_score(label, pred)
    print(auc)
    end = time.process_time()
    print("程序运行时间：",(end-start))
    
    
    
    


    for r in range(0,2):
        start=time.process_time()
        ma=markov.Markov()
        kmer_propos=ma.get_Mulk_Mul_kmer_Pro(pos,kstart,kend,r)
        kmer_proneg=ma.get_Mulk_Mul_kmer_Pro(neg,kstart,kend,r)
        end=time.process_time()
        print("计算马尔可夫概率时间：",(end-start))
        print("--------------------------MulD2s-----------------------")
        print("kstart=",kstart," kend=",kend," r=",r)
        start = time.process_time()
        ## 计算MulD2特征
        start1=time.process_time()
        posD2Lis=sq.getD2SMulCount_nonorm(pos,kstart,kend,r,flag,kmer_propos) 
        negD2Lis=sq.getD2SMulCount_nonorm(neg,kstart,kend,r,flag,kmer_proneg)
        end1=time.process_time()
        print("特征计算时间",(end1-start1))
        ## 构成基因对,并设置标识
        posPair=[]
        negPair=[]
        for i in range(len(posD2Lis)):
            for j in range(i+1,len(posD2Lis)):
                posSim=Sim.getMulD2sWeightSim2(posD2Lis[i],posD2Lis[j],posweight)
                negSim=Sim.getMulD2sWeightSim2(negD2Lis[i],negD2Lis[j],negweight)
                posPair.append([1,posSim])
                negPair.append([0,negSim])
                
        Pair=posPair+negPair
        pSim=sorted(Pair,key=getSim,reverse=True)
        print("数据集大小：",len(pSim))
        PP=np.array(pSim)
        label=PP[:,0]
        pred=PP[:,1]    
        auc=roc_auc_score(label, pred)
        print(auc)
        end = time.process_time()
        print("程序运行时间：",(end-start)) 
        
        
        print("--------------------------MulD2star-----------------------")
        print("kstart=",kstart," kend=",kend," r=",r)
        start = time.process_time()
        ## 计算MulD2特征
        posD2Lis=sq.getD2StarMulCount_mid(pos,kstart,kend,r,flag,kmer_propos) 
        negD2Lis=sq.getD2StarMulCount_mid(neg,kstart,kend,r,flag,kmer_proneg)
        ## 构成基因对,并设置标识
        posPair=[]
        negPair=[]
        for i in range(len(posD2Lis)):
            for j in range(i+1,len(posD2Lis)):
                posSim=Sim.getMulD2starWeightSim2(posD2Lis[i],posD2Lis[j],posweight)
                negSim=Sim.getMulD2starWeightSim2(negD2Lis[i],negD2Lis[j],negweight)
                posPair.append([1,posSim])
                negPair.append([0,negSim])
        Pair=posPair+negPair
        pSim=sorted(Pair,key=getSim,reverse=True)
        print("数据集大小：",len(pSim))
        PP=np.array(pSim)
        label=PP[:,0]
        pred=PP[:,1]    
        auc=roc_auc_score(label, pred)
        print(auc)
        end = time.process_time()
        print("程序运行时间：",(end-start))