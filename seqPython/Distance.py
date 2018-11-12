# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 16:04:51 2018

@author: Yzi
"""

import Sequence
import math
import numpy as np
import sys
class Distance:
    ## 欧式距离
    def EuD(self,seqA,seqB,k):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        kmerSet,dic =Sq.getSeqKerSet(seqLis,k)
        lisFea,freq=Sq.getSeqfreq(seqLis,k,dic)
        su=0.0
        for key in kmerSet:
            su=(lisFea[0][key]-lisFea[1][key])**2+su
        return math.sqrt(su)
    ## 曼哈顿距离
    def manhattan(self,seqA,seqB,k):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        kmerSet,dic =Sq.getSeqKerSet(seqLis,k)
        lisFea,freq=Sq.getSeqfreq(seqLis,k,dic)
        su=0.0
        for key in kmerSet:
            su=abs(lisFea[0][key]-lisFea[1][key])+su
        return su
    
     ## KLD 相对熵距离
    def KLD(self,seqA,seqB,k):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        kmerSet,dic =Sq.getSeqKerSet(seqLis,k)
        lisFea,freq=Sq.getSeqfreq(seqLis,k,dic)
        k1=0.0
        k2=0.0
        for j in range(freq.shape[1]):
            k1=k1+freq[0,j]*math.log2(freq[0,j]/freq[1,j])
            k2=k2+freq[1,j]*math.log2(freq[1,j]/freq[0,j])
        return (k1+k2)/2
        
        ## KLD 相对熵距离
    def pcc(self,seqA,seqB,k):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        kmerSet,dic =Sq.getSeqKerSet(seqLis,k)
        lisFea,freq=Sq.getSeqfreq(seqLis,k,dic)
        #计算平均值
        meanA=0.0
        meanB=0.0
        for i in range(freq.shape[1]):
            meanA=freq[0,i]+meanA
            meanB=freq[1,i]+meanB
        meanA=meanA/freq.shape[1]
        meanB=meanB/freq.shape[1]
        #计算协方差
        cov=0.0
        for i in range(freq.shape[1]):
            cov=cov+(freq[0,i]-meanA)*(freq[1,i]-meanB)
        cov=cov/(freq.shape[1]-1)
        # 计算方差
        stA=np.std(freq[0,:],ddof=1)
        stB=np.std(freq[1,:],ddof=1)
        #计算pcc
        pcc=cov/(stA*stB)
        return abs(1/pcc)
    
    def cosine(self,seqA,seqB,k):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        kmerSet,dic =Sq.getSeqKerSet(seqLis,k)
        lisFea,freq=Sq.getSeqfreq(seqLis,k,dic)
        #计算余弦距离
        LA=0.0
        LB=0.0
        su=0.0
        for i in range(freq.shape[1]):
            su=su+freq[0,i]*freq[1,i]
            LA=freq[0,i]**2+LA
            LB=freq[1,i]**2+LB
        cosine=su/(math.sqrt(LA)*math.sqrt(LB))
        return abs(1/cosine)
    
    def chebyshev(self,seqA,seqB,k):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        kmerSet,dic =Sq.getSeqKerSet(seqLis,k)
        lisFea,freq=Sq.getSeqfreq(seqLis,k,dic)
        #计算切比雪夫距离
        ma=-sys.maxsize-1
        for i in range(freq.shape[1]):
            ma=max(ma,abs(freq[0,i]-freq[1,i]))
        return ma
    
    def D2(self,seqA,seqB,k):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        kmerSet,dic =Sq.getSeqKerSet(seqLis,k)
        lisFea,count=Sq.getSeqCount(seqLis,k,dic)
        #计算D2
        su=0.0
        for i in range(count.shape[1]):
            su=su+count[0,i]*count[1,i]
        return su
    
    
    
    
if __name__=="__main__":
    seqA="ATCG"
    seqB="TCAG"
    distance =Distance()
    su=distance.EuD(seqA,seqB,2)
    print(su)
    su=distance.manhattan(seqA,seqB,2)
    print(su)
    print(math.log2(11))
    su=distance.KLD(seqA,seqB,2)
    print(su)
    su=distance.pcc(seqA,seqB,2)
    print(su)
    su=distance.cosine(seqA,seqB,2)
    print(su)
    su=distance.chebyshev(seqA,seqB,2)
    print(su)
    su=distance.D2(seqA,seqB,2)
    print(su)