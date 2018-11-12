# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 13:29:52 2018
@author: Yzi

"""
import numpy as np
class Sequence:
    chDic=["A","T","C","G"]
    ## 求 序列kmer集合：返回集合和字典  序列长度n tc=O(n)
    def getSingleSeqKerSet(self,sequence,k):
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
    
    ## 统计kmer 在所有数据出现的次数,返回字典 tc=O(mn)
    def getkmerCount(self,sequences,dic,k):
        lis,count=self.getSeqCount(sequences,k,dic)
        resultdic =dict.copy(lis[0])
        for i in range(1,len(lis)):
            for key in lis[i].keys():
                resultdic[key]=resultdic[key]+lis[i][key]
        return resultdic
                
    ## 为序列去掉末尾元素
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

            
if __name__ =="__main__":
    sequece="ATCTCTAAAGGGA"
    sequences=["ATCG","CGATC","ACTGTT"]
   
    SequenceTest=Sequence()
    print(SequenceTest.trimSequence(sequences))
    # 字典集合
    kmerset,dic=SequenceTest.getSeqKerSet(sequences,3)
    lis,freq=SequenceTest.getSeqfreq(sequences,3,dic)
    print(lis)
    #初始化
