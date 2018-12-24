# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 10:08:16 2018

@author: Yzi
"""
###引入训练集求概率

import ReadData
import Sequence
import Similarity
from sklearn import metrics
import time
import markov
import random
if __name__=="__main__":
    rd=ReadData.ReadData()
    query,dataset,label=rd.getData("dataset1")
    sim=Similarity.Similarity()
    ## 随机打乱
    randnum = random.randint(0,100)
    random.seed(randnum)
    random.shuffle(dataset)
    random.seed(randnum)
    random.shuffle(label)
    ### 划分 训练集和测试集
    trainset=[]
    trainlabel=[]
    testset=[]
    testlabel=[]
    poscount=0
    negcount=0
#    counttrain=0
#    counttest=0
    
    for i in range(len(dataset)):
        if (label[i]==1 and poscount<15) or (label[i]==0 and negcount<14):
            trainset.append(dataset[i])
            trainlabel.append(label[i])
#            counttrain=counttrain+1
            if label[i]==1:
                poscount=poscount+1
            else:
                negcount=negcount+1
        else:
            testset.append(dataset[i])
            testlabel.append(label[i])
    trainSets=list.copy(trainset)
    trainSets.append(query[0])
    testSets=list.copy(testset)
    testSets.append(query[0])

    sq=Sequence.Sequence()
    datasets=list.copy(dataset)
    datasets.append(query[0])
#    print(label)
    kstart=2
    kend=5
    flag=False
    start=time.process_time()
    ## 改成 训练集得权重
    weight = sq.getMulWeight_mid(trainSets,kstart,kend)
    testkmerweight=sq.getMulSeqKerSet(testSets,kstart,kend)
    for key in testkmerweight:
        if key in weight.keys():
            testkmerweight[key]=weight[key]
    weight=testkmerweight
    end=time.process_time()
    print("计算权重时间：",(end-start))
    print(len(weight))
    print("--------------MulD2Weight相似度方法------------")
    
   
    
    start = time.process_time()
    pred=[]
    ## 变成测试集
    for data in testset:
        si=sim.getMulD2WeightSim(query[0],data,kstart,kend,testSets,weight)
        pred.append(si)
    
#    修改成testeLabel
    fpr, tpr, thresholds = metrics.roc_curve(testlabel,pred,pos_label=1)
    auc=metrics.auc(fpr, tpr)  
    
    end = time.process_time()
    print("程序运行时间：",(end-start))
    print("kstart=",kstart,",kend=",kend," eur.auc=:")
    print(auc)
    
    
    for r in range(0,3):
        ## 计算kmer概率
        start=time.process_time()
        ma=markov.Markov()
        ## 修改 成训练集得到的概率
        kmer_pro=ma.get_Mulk_Mul_kmer_Pro(trainSets,kstart,kend,r)
        kmer_test=sq.getMulSeqKerSet(testSets,kstart,kend)
        for key in kmer_test:
            if key in kmer_pro.keys():
                kmer_test[key]=kmer_pro[key]
        kmer_pro=kmer_test
        end=time.process_time()
        print("计算马尔可夫概率时间：",(end-start))
        
        
        print("--------------MulD2S相似度方法------------")
        start = time.process_time()
        pred=[]
        for data in testset:
            slis=[]
            slis.append(query[0])
            slis.append(data)
            si=sim.getMulD2sWeightSim(query[0],data,kstart,kend,r,flag,testSets,weight,kmer_pro)
            pred.append(si)
        fpr, tpr, thresholds = metrics.roc_curve(testlabel, pred,pos_label=1)
        auc=metrics.auc(fpr, tpr)    
        end = time.process_time()
        print("程序运行时间：",(end-start))
        print("kstart=",kstart,",kend=",kend,"r=",r," d2sMul.auc=:")
        print(auc)
        
        
        
        print("--------------MulD2star相似度方法------------")
        start = time.process_time()
        pred=[]
        for data in testset:
            slis=[]
            slis.append(query[0])
            slis.append(data)
            si=sim.getMulD2starWeightSim(query[0],data,kstart,kend,r,flag,testSets,weight,kmer_pro)
            pred.append(si)
        fpr, tpr, thresholds = metrics.roc_curve(testlabel, pred,pos_label=1)
        auc=metrics.auc(fpr, tpr)  
        end = time.process_time()
        print("程序运行时间：",(end-start))
        print("kstart=",kstart,",kend=",kend,"r=",r," d2starMul.auc=:")
        print(auc)
    
    