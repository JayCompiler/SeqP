# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 09:59:01 2018

@author: Yzi

"""

import ReadData
import Sequence
import Similarity
from sklearn import metrics
import time
import markov
if __name__=="__main__":
    rd=ReadData.ReadData()
    query,dataset,label=rd.getData("dataset1")
    sim=Similarity.Similarity()
    sq=Sequence.Sequence()
    datasets=list.copy(dataset)
    datasets.append(query[0])
    kstart=2
    kend=2
    flag=False
    start=time.process_time()
    weight = sq.getMulWeight(datasets,kstart,kend)
    end=time.process_time()
    print("计算权重时间：",(end-start))
    print(len(weight))
    print("--------------MulD2Weight相似度方法------------")
    
   
    
    start = time.process_time()
    pred=[]
    for data in dataset:
        si=sim.getMulD2WeightSim(query[0],data,kstart,kend,datasets,weight)
        pred.append(si)
#    print(pred)
    fpr, tpr, thresholds = metrics.roc_curve(label, pred,pos_label=1)
    auc=metrics.auc(fpr, tpr)  
    
    end = time.process_time()
    print("程序运行时间：",(end-start))
    print("kstart=",kstart,",kend=",kend," eur.auc=:")
    print(auc)
    
    
    for r in range(0,3):
        ## 计算kmer概率
        start=time.process_time()
        ma=markov.Markov()
        kmer_pro=ma.get_Mulk_Mul_kmer_Pro(datasets,kstart,kend,r)
        end=time.process_time()
        print("计算马尔可夫概率时间：",(end-start))
        
        
        print("--------------MulD2S相似度方法------------")
        start = time.process_time()
        pred=[]
        for data in dataset:
            slis=[]
            slis.append(query[0])
            slis.append(data)
            si=sim.getMulD2sWeightSim(query[0],data,kstart,kend,r,flag,datasets,weight,kmer_pro)
            pred.append(si)
        fpr, tpr, thresholds = metrics.roc_curve(label, pred,pos_label=1)
        auc=metrics.auc(fpr, tpr)    
        end = time.process_time()
        print("程序运行时间：",(end-start))
        print("kstart=",kstart,",kend=",kend,"r=",r," d2sMul.auc=:")
        print(auc)
        print("--------------MulD2star相似度方法------------")
        start = time.process_time()
        pred=[]
        for data in dataset:
            slis=[]
            slis.append(query[0])
            slis.append(data)
            si=sim.getMulD2starWeightSim(query[0],data,kstart,kend,r,flag,datasets,weight,kmer_pro)
            pred.append(si)
        fpr, tpr, thresholds = metrics.roc_curve(label, pred,pos_label=1)
        auc=metrics.auc(fpr, tpr)  
        end = time.process_time()
        print("程序运行时间：",(end-start))
        print("kstart=",kstart,",kend=",kend,"r=",r," d2starMul.auc=:")
        print(auc)
    
    