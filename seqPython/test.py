# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 14:42:05 2018

@author: Yzi
"""
import ReadData
import Sequence
import Similarity
from sklearn import metrics
import numpy as np
if __name__=="__main__":

     a=np.zeros((3,4))
     print(a)
     b=np.ones((3,3))
     print(b)
     for i in range(3):
         pass
         
#    rd=ReadData.ReadData()
#    query,dataset,label=rd.getData("dataset1")
#    sim=Similarity.Similarity()
#    sq=Sequence.Sequence()
#    datasets=list.copy(dataset)
#    datasets.append(query[0])
# 
#    for k in range(2,11):
#        ## 欧式距离
#        kemrset,kmersetdic=sq.getSeqKerSet(datasets,k)    
#        freqLis,freq=sq.getSeqfreq(datasets,k,kmersetdic)
#        weight = sq.getWeight(freqLis)
#         ## D2距离
#        print("--------------D2Weight相似度方法------------")
#        pred=[]
#        for data in dataset:
#            si=sim.getD2WeightSim(query[0],data,k,datasets,weight)
#            pred.append(si)
#        fpr, tpr, thresholds = metrics.roc_curve(label, pred,pos_label=1)
#        auc=metrics.auc(fpr, tpr)    
#        print("k=",k," eur.auc=:")
#        print(auc)