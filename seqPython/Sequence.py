# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 13:29:52 2018
@author: Yzi

"""
import numpy as np
import markov as mk
class Sequence:
    chDic=["A","T","C","G"]
    ## 求 序列kmer集合：返回集合和字典  序列长度n tc=O(n)
    def getSingleSeqKerSet(self,sequence,k):
        if k>len(sequence):
            print("k的长度大于序列长度，不合法，请重新设置")
        if k<=0:
            print("k的设定小于等于0，不合法")
            return 
        kmerlist=[];
        for i in range(len(sequence)-k+1):
            kmerlist.append(sequence[i:i+k])
        kmerSet=set(kmerlist)
        dic=dict.fromkeys(kmerSet,0)
        return kmerSet,dic
    
    ## 求数据集合序列kmer集合：集合转换成字典返回  序列长度n，序列个数m tc=O(mn)
    def getSeqKerSet(self,sequences,k):
        kmerSet=set()
        for i in sequences:
            tmpset,no=self.getSingleSeqKerSet(i,k)
            #集合的并集
            kmerSet=kmerSet|tmpset
        #初始化字典
        dic=dict.fromkeys(kmerSet,0)
        return kmerSet,dic
    
    ## 求数据集合序列kmer集合：集合转换成字典返回  序列长度n，序列个数m tc=O((kend-ksatrt)*mn)
    def getMulSeqKerSet(self,sequences,kstart,kend):
        kmerset,dic=self.getSeqKerSet(sequences,kstart)
        for k in range(kstart+1,kend+1):
            tmpset,tmpdic=self.getSeqKerSet(sequences,k)
            dic=dict(dic,**tmpdic)
        return dic
    
    ## 统计词计数   输入序列list，k和字典  序列长度n，序列个数m  tc=O(mn)
    def getSeqCount(self,sequences,k,dic):
        lis=[]
        for sequence in sequences:
            # 字典的浅复制
            tmpdic=dict.copy(dic)
            for i in range(len(sequence)-k+1): 
                tmpdic[sequence[i:i+k]]= tmpdic[sequence[i:i+k]]+1
            lis.append(tmpdic)
        # 初始化数组
        count=np.zeros(shape=(len(sequences),len(dic)))
        for i in range(len(sequences)):
            co=0
            for dc in dic.keys():
                count[i][co]=lis[i][dc]
                co=co+1
        return lis,count+np.spacing(1)
    
    
     # 获取多个k 的频数 时间复杂度 为O((kend-kstart)*mn)
    def getMulCount(self,sequences,kstart,kend):
         kmerset,dic=self.getSeqKerSet(sequences,kstart)
         countLis,count=self.getSeqCount(sequences,kstart,dic)
         for k in range(kstart+1,kend+1):
            kmerset,dic=self.getSeqKerSet(sequences,k)
            tmpcountLis,count=self.getSeqCount(sequences,k,dic)
            for i in range(len(countLis)):
                countLis[i]=dict(countLis[i],**(tmpcountLis[i]))
         return countLis
    
    ## 统计kmer 在所有数据出现的次数,返回字典 tc=O(mn)
    def getkmerCount(self,sequences,dic,k):
        lis,count=self.getSeqCount(sequences,k,dic)
        resultdic =dict.copy(lis[0])
        for i in range(1,len(lis)):
            for key in lis[i].keys():
                resultdic[key]=resultdic[key]+lis[i][key]
        return resultdic
                
    ## 为序列去掉末尾元素 tc=O(n)
    def trimSequence(self,sequences):
        result=list.copy(sequences)
        for i in range(len(sequences)):
            result[i]=sequences[i][0:-1]
        return result
    
    ## 统计词频 输入序列集合，k和字典 序列长度n，序列个数m  tc=O(mn)
    def getSeqfreq(self,sequences,k,dic):
        lis=[]
        for sequence in sequences:
            # 字典的浅复制
            tmpdic=dict.copy(dic)
            #词计数
            for i in range(len(sequence)-k+1): 
                tmpdic[sequence[i:i+k]]= tmpdic[sequence[i:i+k]]+1
            # 词频
            for i in tmpdic.keys(): 
                tmpdic[i]=tmpdic[i]/(len(sequence)-k+1)
            lis.append(tmpdic)
          # 初始化数组
        freq=np.zeros(shape=(len(sequences),len(dic)))
        for i in range(len(sequences)):
            co=0
            for dc in dic.keys():
                freq[i][co]=lis[i][dc]
                co=co+1
        return lis,freq+np.spacing(1)
    
    ## 统计D2S kmersetdic 表示所有序列所有的kmer集合 O(m*n*k)
    def getD2SCount(self,sequence,sequences,k,r,flag,kmersetdic):
        ses=[]
        ses.append(sequence)
        # 获得词统计
        lis,count=self.getSeqCount(ses,k,kmersetdic)
        ma=mk.Markov()
        prodic={}
        if flag==False:
            prodic=ma.get_Single_kmer_Pro(sequence,sequences,k,r)
        else:
            prodic=ma.get_Mul_kmer_Pro(sequence,sequences,k,r)
        n=len(sequence)
        for key in prodic:
            if lis[0][key]==0:
                prodic[key]=0
            else:
                prodic[key]=lis[0][key]-n*prodic[key]
        return self.addfloat(prodic)
    
    
    ## 统计多个k值 D2S kmersetdic 表示所有序列所有的kmer集合 O((kend-kstart)*m*n*k)
    def getD2SMulCount(self,sequence,sequences,kstart,kend,r,flag):
        kmerset,kmersetdic=self.getSeqKerSet(sequences,kstart)
        prodic =self.getD2SCount(sequence,sequences,kstart,r,flag,kmersetdic)
        for k in range(kstart+1,kend+1):
            kmerset,kmersetdic=self.getSeqKerSet(sequences,k)
            tmpprodic =self.getD2SCount(sequence,sequences,k,r,flag,kmersetdic)
            prodic=dict(prodic,**tmpprodic)
        return prodic
    
    # 获取多个k 的频率
    def getMulFreq(self,sequences,kstart,kend):
         kmerset,dic=self.getSeqKerSet(sequences,kstart)
         freqLis,freq=self.getSeqfreq(sequences,kstart,dic)
         for k in range(kstart+1,kend+1):
            kmerset,dic=self.getSeqKerSet(sequences,k)
            tmpfreqLis,freq=self.getSeqfreq(sequences,k,dic)
            for i in range(len(freqLis)):
                freqLis[i]=dict(freqLis[i],**(tmpfreqLis[i]))
         return freqLis
    
    ## 获取权重  freqLis表示频率矩阵 时间复杂度O(n*m^2)
    def getWeight(self,freqLis):
        # 获得字典
        resultdic = dict.copy(freqLis[0]) 
        for key in resultdic:
            tmp=0.0
            for j in range(len(freqLis)):
                for k in range(len(freqLis)):
                    tmp=abs(freqLis[j][key]-freqLis[k][key])+tmp
            resultdic[key]=tmp
        # 计算权重    
        suf=0.0
        for key in resultdic:
            suf=resultdic[key]+suf
        for key in resultdic:
            resultdic[key]=resultdic[key]/suf
        return resultdic
    
    
    
    ## 获取多个k的权重  freqLis表示频率矩阵 时间复杂度O((kend-kstart)*n*m^2)
    def getMulWeight(self,sequences,kstart,kend):
        #获取第一个频率合集合字典 
        kmerset,dic=self.getSeqKerSet(sequences,kstart)
        freqLis,freq=self.getSeqfreq(sequences,kstart,dic)
        for k in range(kstart+1,kend+1):
            kmerset,dic=self.getSeqKerSet(sequences,k)
            tmpfreqLis,freq=self.getSeqfreq(sequences,k,dic)
            for i in range(len(freqLis)):
                freqLis[i]=dict(freqLis[i],**(tmpfreqLis[i]))
        resultdic=self.getWeight(freqLis)
        return resultdic
        
    # 平滑数据
    def addfloat(self,dic):
        tmp=dict.copy(dic)
        for key in dict.keys(tmp):
            tmp[key]=tmp[key]+np.spacing(1)
        return tmp
         

if __name__ =="__main__":
    sequece="ATCTCTAAAGGGA"
    sequences=["ATCCATA","CGACCCC","ACTAA"]
    k=4
    r=1
    ks=2
    ke=4
    SequenceTest=Sequence()
    wei=SequenceTest.getMulWeight(sequences,ks,ke)
    print(wei)
#    print(SequenceTest.trimSequence(sequences))
    # 字典集合
#    kmerset,dic=SequenceTest.getSeqKerSet(sequences,k)
#    lis,freq=SequenceTest.getSeqfreq(sequences,k,dic)
#    print(lis)
#    print("--------------")
#    prodic=SequenceTest.getD2SCount(sequences[0],sequences,k,r,False,dic)
#    print(prodic)