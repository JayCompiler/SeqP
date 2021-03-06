# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 10:56:32 2018

@author: Yzi
"""
import ReadData
import Sequence
import Similarity
import markov
from sklearn import metrics
if __name__=="__main__":
    rd=ReadData.ReadData()
    query,dataset,label=rd.getData("dataset1")
    sim=Similarity.Similarity()
    sq=Sequence.Sequence()
    datasets=list.copy(dataset)
#    print("---------------------------------")
#    print(datasets)
    datasets.append(query[0])
    

    for k in range(2,11):
        ## 欧式距离
        kemrset,kmersetdic=sq.getSeqKerSet(datasets,k)
        print("--------------欧式相似度方法------------")
        pred=[]
        for data in dataset:
            si=sim.getEuSim(query[0],data,k)
            pred.append(si)
        fpr, tpr, thresholds = metrics.roc_curve(label, pred,pos_label=1)
        auc=metrics.auc(fpr, tpr)    
        print("k=",k," eur.auc=:")
        print(auc)
        
        ## manhattan距离
        print("--------------manhattan相似度方法------------")
        pred=[]
        for data in dataset:
            si=sim.getmanhattanSim(query[0],data,k)
            pred.append(si)
        fpr, tpr, thresholds = metrics.roc_curve(label, pred,pos_label=1)
        auc=metrics.auc(fpr, tpr)    
        print("k=",k," eur.auc=:")
        print(auc)
        
         ## KLD距离
        print("--------------KLD相似度方法------------")
        pred=[]
        for data in dataset:
            si=sim.getKLDSim(query[0],data,k)
            pred.append(si)
        fpr, tpr, thresholds = metrics.roc_curve(label, pred,pos_label=1)
        auc=metrics.auc(fpr, tpr)    
        print("k=",k," eur.auc=:")
        print(auc)
        
        ## pcc距离
        print("--------------pcc相似度方法------------")
        pred=[]
        for data in dataset:
            si=sim.getpccSim(query[0],data,k)
            pred.append(si)
        fpr, tpr, thresholds = metrics.roc_curve(label, pred,pos_label=1)
        auc=metrics.auc(fpr, tpr)    
        print("k=",k," eur.auc=:")
        print(auc)
         
        ## 余弦距离
        print("--------------余弦相似度方法------------")
        pred=[]
        for data in dataset:
            si=sim.getcosineSim(query[0],data,k)
            pred.append(si)
        fpr, tpr, thresholds = metrics.roc_curve(label, pred,pos_label=1)
        auc=metrics.auc(fpr, tpr)    
        print("k=",k," eur.auc=:")
        print(auc)
        
         ## chebyshev距离
        print("--------------chebyshev相似度方法------------")
        pred=[]
        for data in dataset:
            si=sim.getchebyshevSim(query[0],data,k)
            pred.append(si)
        fpr, tpr, thresholds = metrics.roc_curve(label, pred,pos_label=1)
        auc=metrics.auc(fpr, tpr)    
        print("k=",k," eur.auc=:")
        print(auc)
        
         ## D2距离
        print("--------------D2相似度方法------------")
        pred=[]
        for data in dataset:
            si=sim.getD2Sim(query[0],data,k)
            pred.append(si)
        fpr, tpr, thresholds = metrics.roc_curve(label, pred,pos_label=1)
        auc=metrics.auc(fpr, tpr)    
        print("k=",k," eur.auc=:")
        print(auc)
        
        
       
        for r in range(0,3):
             ## 计算马尔可夫概率
            ma=markov.Markov()
            kmer_pro=ma.get_Mul_kmer_Pro(datasets,k,r)
            
            print("--------------D2S相似度方法------------")
            pred=[]
            for data in dataset:
                slis=[]
                slis.append(query[0])
                slis.append(data)
#                kemrset,kmersetdic=sq.getSeqKerSet(slis,k)
#                kemrset,kmersetdic=sq.getSeqKerSet(datasets,k)
#                si=sim.getD2sSim(query[0],data,k,r,False,kmersetdic)
                si=sim.getD2sSim(query[0],data,k,r,True,kmersetdic,kmer_pro)

                pred.append(si)
            fpr, tpr, thresholds = metrics.roc_curve(label, pred,pos_label=1)
            auc=metrics.auc(fpr, tpr)    
            print("k=",k,"r=",r," eur.auc=:")
            print(auc)
            print("--------------D2star相似度方法------------")
            pred=[]
            for data in dataset:
                slis=[]
                slis.append(query[0])
                slis.append(data)
#                kemrset,kmersetdic=sq.getSeqKerSet(slis,k)
#                kemrset,kmersetdic=sq.getSeqKerSet(datasets,k)
#                si=sim.getD2starSim(query[0],data,k,r,False,datasets,kmersetdic)
                si=sim.getD2starSim(query[0],data,k,r,True,datasets,kmersetdic,kmer_pro)

                pred.append(si)
            fpr, tpr, thresholds = metrics.roc_curve(label, pred,pos_label=1)
            auc=metrics.auc(fpr, tpr)    
            print("k=",k,"r=",r," auc=:")
            print(auc)
        
        