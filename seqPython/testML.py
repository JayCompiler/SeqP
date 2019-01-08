# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 14:38:22 2019

@author: Yzi
"""

import random

from deap import base
from deap import creator
from deap import tools
#import numpy as np
import ReadData
import Sequence
from sklearn.metrics import roc_auc_score
import Distance


#name="human_muscle"
#name="fly_blastoderm"
name="human_HBB"
k=2
## 获取数据集 整个数据集，正数据集，负数据集 都是序列，没有标签
#datasets,pos,neg=rd.getData2("fly_blastoderm")
  ## 获取数据 kmer 从2-->6
print("---------",name,"------------","k=",k)
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


## d2 以整个数据集计算 计算频率,用于计算单个k的其他距离定义方法
d2freLis,d2arf=sq.getSeqfreq(datasets,2,d2dic)
d3freLis,d3arf=sq.getSeqfreq(datasets,3,d3dic)
d4freLis,d4arf=sq.getSeqfreq(datasets,4,d4dic)
d5freLis,d5arf=sq.getSeqfreq(datasets,5,d5dic)
d6freLis,d6arf=sq.getSeqfreq(datasets,6,d6dic)


### 将2-6 合成一个字典,并计算权重----------------------------------------------------------------

size=len(d2freLis)
## 标准化后合并
freqLis=[None]*size
for i in range(len(d2freLis)):
    freqLis[i]=dict(sq.normdata_max_min(d2freLis[i]),**(sq.normdata_max_min(d3freLis[i])))
    freqLis[i]=dict(freqLis[i],**(sq.normdata_max_min(d4freLis[i])))
    freqLis[i]=dict(freqLis[i],**(sq.normdata_max_min(d5freLis[i])))
    freqLis[i]=dict(freqLis[i],**(sq.normdata_max_min(d6freLis[i])))

# 计算权重
weightMax=sq.getWeight(freqLis)


    ## d2 特征 count 所有数据集
d2countLis,d2arc=sq.getSeqCount(datasets,2,d2dic)
d3countLis,d3arc=sq.getSeqCount(datasets,3,d3dic)
d4countLis,d4arc=sq.getSeqCount(datasets,4,d4dic)
d5countLis,d5arc=sq.getSeqCount(datasets,5,d5dic)
d6countLis,d6arc=sq.getSeqCount(datasets,6,d6dic)

## 将 字典 离散 二值化
for i in range(len(d2countLis)):
    for key in dict.keys(d2countLis[0]):
        if d2countLis[i][key]>0:
            d2countLis[i][key]=1
#        if d2countLis[i][key]!=0 and d2countLis[i][key]!=1:
#            print(d2countLis[i][key])
for i in range(len(d3countLis)):
    for key in dict.keys(d3countLis[0]):
        if d3countLis[i][key]>0:
            d3countLis[i][key]=1
for i in range(len(d4countLis)):
    for key in dict.keys(d4countLis[0]):
        if d4countLis[i][key]>0:
            d4countLis[i][key]=1
for i in range(len(d5countLis)):
    for key in dict.keys(d5countLis[0]):
        if d5countLis[i][key]>0:
            d5countLis[i][key]=1
for i in range(len(d6countLis)):
    for key in dict.keys(d6countLis[0]):
        if d6countLis[i][key]>0:
            d6countLis[i][key]=1


## 使用同或运算构造特征数据集合
    ## 标准化后合并
countLis=[None]*size
for i in range(len(d2freLis)):
    countLis[i]=dict(d2countLis[i],**(d3countLis[i]))
    countLis[i]=dict(countLis[i],**(d4countLis[i]))
    countLis[i]=dict(countLis[i],**(d5countLis[i]))
    countLis[i]=dict(countLis[i],**(d6countLis[i]))


pairPoslist=[]
pairNeglist=[]
for i in range(possize):
    # 装填正例 同或操作后 保存特征
    for j in range(i+1,possize):
        tpos=dict.copy(countLis[0])
        for key in dict.keys(countLis[0]):
            tpos[key]=(countLis[i][key]^countLis[j][key])^1
        pairPoslist.append(tpos)
        
# 装填负例
for i in range(possize):
    # 装填负例
    for j in range(i+1,possize):       
        tneg=dict.copy(countLis[0])
        for key in dict.keys(countLis[0]):
            tpos[key]=(countLis[i+possize][key]^countLis[j+possize][key])^1
        pairNeglist.append(tneg)
        
## 将lis里面的dict转换成list
pairPoslistTo=[]
pairNeglistTo=[]
for i in range(len(pairPoslist)):
    tmplis=[]
    ## 正例子转换
    for key in sorted(countLis[0]):
        tmplis.append(pairPoslist[i][key])
    pairPoslistTo.append(tmplis)
    ## 负例子转换
    tmplis=[]
    for key in sorted(countLis[0]):
        tmplis.append(pairNeglist[i][key])
    pairNeglistTo.append(tmplis)

print(len(pairNeglistTo))
print(len(pairPoslistTo))      
  


### 产生随机种子，并将数据集打乱
##randnum = random.randint(0,100)
##random.seed(randnum)
### 打乱正负序列集 
#random.shuffle(pairPoslist)
#random.shuffle(pairNeglist)
#
### 划分训练集合和测试集合  并设置标签
#trainLis=[]
#trainlabel=[]
#testLis=[]
#testlabel=[]
#trainPosSize=int(0.8*len(pairPoslist))
#for i in range(len(pairPoslist)):
#    if i<=trainPosSize:
#        trainLis.append(pairPoslist[i])
#        trainlabel.append(1)
#        trainLis.append(pairNeglist[i])
#        trainlabel.append(0)
#    else:
#        testLis.append(pairPoslist[i])
#        testlabel.append(1)
#        testLis.append(pairNeglist[i])
#        testlabel.append(0)
#
