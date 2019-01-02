# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 09:13:50 2018

@author: Yzi
"""

import random

from deap import base
from deap import creator
from deap import tools
import numpy as np
import ReadData
import Sequence
from sklearn.metrics import roc_auc_score
import Distance


#name="human_muscle"
#name="fly_blastoderm"
name="human_HBB"

## 获取数据集 整个数据集，正数据集，负数据集 都是序列，没有标签
#datasets,pos,neg=rd.getData2("fly_blastoderm")
  ## 获取数据 kmer 从2-->6
print("---------",name,"------------")
rd=ReadData.ReadData()

datasets,pos,neg=rd.getData2(name)
## 正例数据集条数
possize=len(pos)
## 合并两个list
datasets=pos+neg


sq=Sequence.Sequence()

flag=True
    ## 获取整个数据集的kmer集合,dict
    
    ## 获取字典集合
d2set,d2dic=sq.getSeqKerSet(datasets,2)
d3set,d3dic=sq.getSeqKerSet(datasets,3)
d4set,d4dic=sq.getSeqKerSet(datasets,4)
d5set,d5dic=sq.getSeqKerSet(datasets,5)
d6set,d6dic=sq.getSeqKerSet(datasets,6)


## d2 以整个数据集计算 计算频率
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
    freqLis[i]=dict(sq.normdata(d2freLis[i]),**(sq.normdata(d3freLis[i])))
    freqLis[i]=dict(freqLis[i],**(sq.normdata(d4freLis[i])))
    freqLis[i]=dict(freqLis[i],**(sq.normdata(d5freLis[i])))
    freqLis[i]=dict(freqLis[i],**(sq.normdata(d6freLis[i])))

# 计算权重
weight=sq.getWeight(freqLis)


    ## d2 特征 count 所有数据集
d2countLis,d2arc=sq.getSeqCount(datasets,2,d2dic)
d3countLis,d3arc=sq.getSeqCount(datasets,3,d3dic)
d4countLis,d4arc=sq.getSeqCount(datasets,4,d4dic)
d5countLis,d5arc=sq.getSeqCount(datasets,5,d5dic)
d6countLis,d6arc=sq.getSeqCount(datasets,6,d6dic)



## 标准化过程： 标准化 count
for i in range(len(datasets)):
    d2countLis[i]=sq.normdata(d2countLis[i])
    d3countLis[i]=sq.normdata(d3countLis[i])
    d4countLis[i]=sq.normdata(d4countLis[i])
    d5countLis[i]=sq.normdata(d5countLis[i])
    d6countLis[i]=sq.normdata(d6countLis[i])
    
    
## 标准化后合并
countLis=[None]*size
for i in range(len(d2freLis)):
    countLis[i]=dict(d2freLis[i],**(d3freLis[i]))
    countLis[i]=dict(freqLis[i],**(sq.normdata(d4freLis[i])))
    countLis[i]=dict(sq.normdata(freqLis[i]),**(sq.normdata(d5freLis[i])))
    countLis[i]=dict(sq.normdata(freqLis[i]),**(sq.normdata(d6freLis[i])))

    

## 构造训练集合和测试集合
# 构造对集合
pairPoslist=[]
pairNeglist=[]
for i in range(possize):
    # 装填正例
    for j in range(i+1,possize):
        tpos=[i,j,1]
        pairPoslist.append(tpos)
# 装填负例
for i in range(possize):
    # 装填正例
    for j in range(i+1,possize):
        tneg=[i+possize,j+possize,0]
        pairNeglist.append(tneg)

print(pairPoslist)
print(len(pairPoslist))
print(pairNeglist)
print(len(pairNeglist))

## 产生随机种子，并将数据集打乱
#randnum = random.randint(0,100)
#random.seed(randnum)
## 打乱正负序列集 
random.shuffle(pos)
random.shuffle(neg)


## 划分训练集合和测试集合  并设置标签
trainLis=[]
trainlabel=[]
testLis=[]
testlabel=[]
trainPosSize=int(0.8*len(pairPoslist))
for i in range(len(pairPoslist)):
    if i<=trainPosSize:
        trainLis.append(pairPoslist[i])
        trainlabel.append(1)
        trainLis.append(pairNeglist[i])
        trainlabel.append(0)
    else:
        testLis.append(pairPoslist[i])
        testlabel.append(1)
        testLis.append(pairNeglist[i])
        testlabel.append(0)



    


    









