# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 15:19:58 2019
第二个实验的经典算法的auc
@author:Yzi
"""

#import random

from deap import base
from deap import creator
from deap import tools
#import numpy as np
import ReadData
import Sequence
from sklearn.metrics import roc_auc_score
import Distance
import markov


name="human_muscle"
#name="fly_blastoderm"
#name="human_HBB"
#k=6
## 获取数据集 整个数据集，正数据集，负数据集 都是序列，没有标签
#datasets,pos,neg=rd.getData2("fly_blastoderm")
  ## 获取数据 kmer 从2-->6
#print("---------",name,"------------","k=",k)
rd=ReadData.ReadData()

datasets,pos,neg=rd.getData2(name)
## 正例数据集条数
possize=len(pos)
## 合并两个list
datasets=pos+neg
sizePo=len(pos)

top=int(sizePo*(sizePo-1))
sq=Sequence.Sequence()

flag=True
    ## 获取整个数据集的kmer集合,dict
    
    ## 获取字典集合
d2set,d2dic=sq.getSeqKerSet(datasets,2)
d3set,d3dic=sq.getSeqKerSet(datasets,3)
d4set,d4dic=sq.getSeqKerSet(datasets,4)
d5set,d5dic=sq.getSeqKerSet(datasets,5)
d6set,d6dic=sq.getSeqKerSet(datasets,6)


## d2 以整个数据集计算 计算频率 不需要标准化，用于计算单个k的其他距离定义方法
d2freLis,d2arf=sq.getSeqfreq(datasets,2,d2dic)
d3freLis,d3arf=sq.getSeqfreq(datasets,3,d3dic)
d4freLis,d4arf=sq.getSeqfreq(datasets,4,d4dic)
d5freLis,d5arf=sq.getSeqfreq(datasets,5,d5dic)
d6freLis,d6arf=sq.getSeqfreq(datasets,6,d6dic)




### 将2-6 合成一个字典,并计算权重

size=len(d2freLis)
## 标准化后合并
freqLis=[None]*size
for i in range(len(d2freLis)):
    freqLis[i]=dict(sq.normdata_max_min(d2freLis[i]),**(sq.normdata_max_min(d3freLis[i])))
    freqLis[i]=dict(freqLis[i],**(sq.normdata_max_min(d4freLis[i])))
    freqLis[i]=dict(freqLis[i],**(sq.normdata_max_min(d5freLis[i])))
    freqLis[i]=dict(freqLis[i],**(sq.normdata_max_min(d6freLis[i])))


    ## d2 特征 count 所有数据集--------------------------------
d2countLis,d2arc=sq.getSeqCount(datasets,2,d2dic)
d3countLis,d3arc=sq.getSeqCount(datasets,3,d3dic)
d4countLis,d4arc=sq.getSeqCount(datasets,4,d4dic)
d5countLis,d5arc=sq.getSeqCount(datasets,5,d5dic)
d6countLis,d6arc=sq.getSeqCount(datasets,6,d6dic)






pairPoslist=[]
poslabel=[]
pairNeglist=[]
neglabel=[]
for i in range(possize):
    # 装填正例
    for j in range(i+1,possize):
        tpos=[i,j,1]
        pairPoslist.append(tpos)
        poslabel.append(1)
# 装填负例
for i in range(possize):
    # 装填正例
    for j in range(i+1,possize):
        tneg=[i+possize,j+possize,0]
        pairNeglist.append(tneg)
        neglabel.append(0)

testLis=pairPoslist+pairNeglist        
testlabel=poslabel+neglabel

if __name__ == "__main__":
    
    dis =Distance.Distance()
    print("数据集大小：",top)
    for k in range(2,7):
        ## eu:
       
        feature="d"+str(k)+"freLis"
        print("----------------EU-----------------")
        print("k=",k)
        sim=[]
        for i in range(len(testLis)):
            tmp=1/(dis.EuD_feature(eval(feature)[testLis[i][0]],eval(feature)[testLis[i][1]]))
            sim.append(tmp)
        auc=roc_auc_score(testlabel, sim)  
        print(auc)
        
        print("----------------Ma-----------------")
        print("k=",k)
        sim=[]
        for i in range(len(testLis)):
            tmp=1/(dis.manhattan_feature(eval(feature)[testLis[i][0]],eval(feature)[testLis[i][1]]))
            sim.append(tmp)
        auc=roc_auc_score(testlabel, sim)  
        print(auc)
        
        
        
        print("----------------kld-----------------")
        print("k=",k)
        sim=[]
        for i in range(len(testLis)):
            tmp=1/(dis.KLD_feature(eval(feature)[testLis[i][0]],eval(feature)[testLis[i][1]]))
            sim.append(tmp)
        auc=roc_auc_score(testlabel, sim)  
        print(auc)
        
        
        print("----------------pcc-----------------")
        print("k=",k)
        sim=[]
        for i in range(len(testLis)):
            tmp=1/(dis.pcc_feature(eval(feature)[testLis[i][0]],eval(feature)[testLis[i][1]]))
            sim.append(tmp)
        auc=roc_auc_score(testlabel, sim)  
        print(auc)
        
    
        print("----------------cosine-----------------")
        print("k=",k)
        sim=[]
        for i in range(len(testLis)):
            tmp=1/(dis.cosine_feature(eval(feature)[testLis[i][0]],eval(feature)[testLis[i][1]]))
            sim.append(tmp)
        auc=roc_auc_score(testlabel, sim)  
        print(auc)
        
        
        print("----------------D2-----------------")
        print("k=",k)
        feature="d"+str(k)+"countLis"
        
        sim=[]
        for i in range(len(testLis)):
            tmp=1/(dis.getD2_feature(eval(feature)[testLis[i][0]],eval(feature)[testLis[i][1]]))
            sim.append(tmp)
        auc=roc_auc_score(testlabel, sim)  
        print(auc)
        
        
        for r in range(0,2):     
              ## 获取kmer 概率集合
            mar=markov.Markov()
            ## 获得概率
            kmer_pro=mar.get_Mulk_Mul_kmer_Pro(datasets,k,k,r)
            dicp="d"+str(k)+"dic"
            #######d2s d2star特征--------------------
            d2scountLis=sq.getD2SMul_SeqCount(datasets,k,r,True,eval(dicp),kmer_pro)
            d2starcountLis=sq.getD2Star_Mul_seq_Count(datasets,k,r,True,eval(dicp),kmer_pro)
            
            print("----------------D2s-----------------")
            print("k=",k)
            feature="d"+str(2)+"scountLis"
            
            
            sim=[]
            for i in range(len(testLis)):
                tmp=1/(dis.getD2S_feature(eval(feature)[testLis[i][0]],eval(feature)[testLis[i][1]]))
                sim.append(tmp)
            auc=roc_auc_score(testlabel, sim)  
            print(auc)
            
            
            print("----------------D2star-----------------")
            print("k=",k)
            feature="d"+str(2)+"starcountLis"
            
            sim=[]
            for i in range(len(testLis)):
                tmp=1/(dis.getD2Star_feature(eval(feature)[testLis[i][0]],eval(feature)[testLis[i][1]]))
                sim.append(tmp)
            auc=roc_auc_score(testlabel, sim)  
            print(auc)
        
        