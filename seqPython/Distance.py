# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 16:04:51 2018

@author: Yzi
"""

import Sequence
import math
import numpy as np
import sys
import markov 
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
    
    def getD2(self,seqA,seqB,k):
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
        return 1/(su+np.spacing(1))
    
    # 输入参数为 序列A，B，kmer长度k，马尔可夫阶数r,是否用全局马尔可夫flag 标志，True表示全局，False表示单条
    # kmerset 表示多少条链参与比较
    def getD2S(self,seqA,seqB,k,r,flag,kmersetdic,kmer_pro):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        kmerSet,dic =Sq.getSeqKerSet(seqLis,k)
        if flag==False:
            lisFeaA=Sq.getD2SCount(seqA,seqLis,k,r,flag,dic)
            lisFeaB=Sq.getD2SCount(seqB,seqLis,k,r,flag,dic)
        else:
            lisFeaA=Sq.getD2SCount(seqA,seqLis,k,r,flag,kmersetdic,kmer_pro)
            lisFeaB=Sq.getD2SCount(seqB,seqLis,k,r,flag,kmersetdic,kmer_pro)
        #计算D2S
        su=0.0
        for key in dict.keys(lisFeaA):
            su=su+(lisFeaA[key]*lisFeaB[key])/math.sqrt(lisFeaA[key]**2+lisFeaB[key]**2)
        return 1/(su+np.spacing(1))
    
    # 输入参数为 序列A，B，kmer长度k，马尔可夫阶数r,是否用全局马尔可夫flag 标志，True表示全局，False表示单条
    # kmerset 表示多少条链参与比较
    def getD2Star(self,seqA,seqB,k,r,flag,sequences,kmersetdic,kmer_pro):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        kmerSet,dic =Sq.getSeqKerSet(seqLis,k)
        # 获取kmer概率
        Ma=markov.Markov()
#        kmerPA={}
#        kmerPB={}
        if flag==False:
            lisFeaA=Sq.getD2SCount(seqA,seqLis,k,r,flag,dic)
            lisFeaB=Sq.getD2SCount(seqB,seqLis,k,r,flag,dic)
#            kmerPA=Ma.get_Single_kmer_Pro(seqA,seqLis,k,r)
#            kmerPB=Ma.get_Single_kmer_Pro(seqB,seqLis,k,r)
        else:
            lisFeaA=Sq.getD2SCount(seqA,seqLis,k,r,flag,kmersetdic,kmer_pro)
            lisFeaB=Sq.getD2SCount(seqB,seqLis,k,r,flag,kmersetdic,kmer_pro)
#            kmerPA=Ma.get_Mul_kmer_Pro(seqA,sequences,k,r)
#            kmerPB=Ma.get_Mul_kmer_Pro(seqB,sequences,k,r)
        #计算D2Star
        su=0.0
        lenA=len(seqA)
        lenB=len(seqB)
        for key in dict.keys(lisFeaA):
            su=su+(lisFeaA[key]*lisFeaB[key])/math.sqrt(lenA*kmer_pro[key]*lenB*kmer_pro[key])
        return 1/(su+np.spacing(1))
    
    ## 带有权重的D2
    def getD2Weight(self,seqA,seqB,k,sequences,weight):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        kmerSet,dic =Sq.getSeqKerSet(sequences,k)
        lisFea,count=Sq.getSeqCount(seqLis,k,dic)
        #计算D2
        su=0.0
        for key in lisFea[0]:
            su=su+lisFea[0][key]*lisFea[1][key]*weight[key]
        return 1/(su+np.spacing(1))
    
    ## 带权重的D2s
    def getD2SWeight(self,seqA,seqB,k,r,flag,kmersetdic,weight,kmer_pro):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        kmerSet,dic =Sq.getSeqKerSet(seqLis,k)
        if flag==False:
            lisFeaA=Sq.getD2SCount(seqA,seqLis,k,r,flag,dic)
            lisFeaB=Sq.getD2SCount(seqB,seqLis,k,r,flag,dic)
        else:
            lisFeaA=Sq.getD2SCount(seqA,seqLis,k,r,flag,kmersetdic,kmer_pro)
            lisFeaB=Sq.getD2SCount(seqB,seqLis,k,r,flag,kmersetdic,kmer_pro)
        #计算D2S
        su=0.0
        for key in dict.keys(lisFeaA):
            su=su+(lisFeaA[key]*lisFeaB[key])/math.sqrt(lisFeaA[key]**2+lisFeaB[key]**2)*weight[key]
        return 1/(su+np.spacing(1))
    
    def getD2StarWeight(self,seqA,seqB,k,r,flag,sequences,kmersetdic,weight,kmer_pro):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        kmerSet,dic =Sq.getSeqKerSet(seqLis,k)
        # 获取kmer概率
        Ma=markov.Markov()
        kmerPA={}
        kmerPB={}
        if flag==False:
            lisFeaA=Sq.getD2SCount(seqA,seqLis,k,r,flag,dic)
            lisFeaB=Sq.getD2SCount(seqB,seqLis,k,r,flag,dic)
            kmerPA=Ma.get_Single_kmer_Pro(seqA,seqLis,k,r)
            kmerPB=Ma.get_Single_kmer_Pro(seqB,seqLis,k,r)
        else:
            lisFeaA=Sq.getD2SCount(seqA,seqLis,k,r,flag,kmersetdic,kmer_pro)
            lisFeaB=Sq.getD2SCount(seqB,seqLis,k,r,flag,kmersetdic,kmer_pro)
#            kmerPA=Ma.get_Mul_kmer_Pro(seqA,sequences,k,r)
#            kmerPB=Ma.get_Mul_kmer_Pro(seqB,sequences,k,r)
        #计算D2Star
        su=0.0
        lenA=len(seqA)
        lenB=len(seqB)
        for key in dict.keys(lisFeaA):
            su=su+(lisFeaA[key]*lisFeaB[key])/math.sqrt(lenA*kmer_pro[key]*lenB*kmer_pro[key])*weight[key]
        return 1/(su+np.spacing(1))
    
    
     ## 多个k值权重的D2
    def getMulD2Weight(self,seqA,seqB,kstart,kend,sequences,weight):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        lisFea=Sq.getMulCount(seqLis,kstart,kend,sequences)
#        print(lisFea)
        #计算D2
        su=0.0
        for key in lisFea[0]:
            su=su+lisFea[0][key]*lisFea[1][key]*weight[key]
        return 1/(su+np.spacing(1))
    
     ## 多个k值权重的D2
    def getMulD2Weight2(self,feaA,feaB,weight):
        su=0.0
        for key in dict.keys(feaA):
            su=su+feaA[key]*feaB[key]*weight[key]
        return 1/(su+np.spacing(1))
    
    
    ## 多个k值的带权重的D2s
    def getMulD2SWeight(self,seqA,seqB,kstart,kend,r,flag,sequences,weight,kmer_pro):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        lisFeaA=Sq.getD2SMulCount(seqA,sequences,kstart,kend,r,flag,kmer_pro)
        lisFeaB=Sq.getD2SMulCount(seqB,sequences,kstart,kend,r,flag,kmer_pro)
        
        #计算D2S
        su=0.0
        for key in dict.keys(lisFeaA):
#            su=su+(lisFeaA[key]*lisFeaB[key])/math.sqrt(lisFeaA[key]**2+lisFeaB[key]**2)*weight[key]
            su=su+(lisFeaA[key]*lisFeaB[key])*weight[key]

        return 1/(su+np.spacing(1))
    
     ## 多个k值的带权重的D2s  
    def getMulD2SWeight2(self,feaA,feaB,weight):
        #计算D2S
        su=0.0
        for key in dict.keys(feaA):
            if feaA[key]**2+feaB[key]**2!=0:
#                su=su+(feaA[key]*feaB[key])/math.sqrt(feaA[key]**2+feaB[key]**2)*weight[key]
                su=su+(feaA[key]*feaB[key])*weight[key]

        return 1/(su+np.spacing(1))
    
    
    def getMulD2StarWeight2(self,feaA,feaB,weight):
        su=0.0
        for key in dict.keys(feaA):
            su=su+(feaA[key]*feaB[key])*weight[key]
            
        return 1/(su+np.spacing(1))
    
    
    def getMulD2StarWeight(self,seqA,seqB,kstart,kend,r,flag,sequences,weight,kmer_pro):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)        
        # 获取 关键字集合 字典dic
        Sq=Sequence.Sequence()
        lisFeaA=Sq.getD2StarMulCount(seqA,sequences,kstart,kend,r,flag,kmer_pro)
        lisFeaB=Sq.getD2StarMulCount(seqB,sequences,kstart,kend,r,flag,kmer_pro)
        su=0.0
        for key in dict.keys(lisFeaA):
            su=su+(lisFeaA[key]*lisFeaB[key])*weight[key]
        return 1/(su+np.spacing(1))
    

    
    
if __name__=="__main__":
    seqA="ATCCATCGCAATATCTCTATGGGAACACTCATCATCTACTATCTAGAGGAG"
#    seqB="ATCCATCGCAATATCTCTATGGGAACACTCATCATCTACTATCTAGAGGAG"
    seqB="ATCCATCGATCTGCCTACATCAGGGCTACAGCTCTGCATTGCAGCAGTCATCAGCAT"
    slis=[]
    slis.append(seqA)
    slis.append(seqB)
    k=5
    r=1
    Seq=Sequence.Sequence()
    kmerset,kmersetdic=Seq.getSeqKerSet(slis,k)
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
    su=distance.getD2(seqA,seqB,2)
    print(su)
    print("-----------")
    su=distance.getD2S(seqA,seqB,k,r,False,kmersetdic)
    print(su)
    su=distance.getD2Star(seqA,seqB,k,r,False,slis,kmersetdic)
    print(su)