# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 13:29:52 2018
@author: Yzi

"""
import numpy as np
import markov as mk
import Statistic
import math
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
    
    
     ## 统计词计数   输入序列list，k和字典  序列长度n，序列个数m  tc=O(mn)，带标准化的count
    def getSeqCount_pre(self,sequences,k,dic):
        lis=[]
        for sequence in sequences:
            # 字典的浅复制
            tmpdic=dict.copy(dic)
            for i in range(len(sequence)-k+1): 
                tmpdic[sequence[i:i+k]]= tmpdic[sequence[i:i+k]]+1
            lis.append(self.normdata(tmpdic))   ### 多加了一句标准化
        # 初始化数组
        count=np.zeros(shape=(len(sequences),len(dic)))
        for i in range(len(sequences)):
            co=0
            for dc in dic.keys():
                count[i][co]=lis[i][dc]
                co=co+1
        return lis,count+np.spacing(1)
    
    
  
    
    
   ## max_min 标准化数据
    def normdata_max_min(self,dic):
        ma=-10000000
        mi=10000000
        for key in dict.keys(dic):
            if dic[key]>ma:
                ma=dic[key]
            if dic[key]<mi:
                mi=dic[key]
        for key in dict.keys(dic):
            dic[key]=(dic[key]-mi)/(ma-mi)
        return dic  
    
    
    
    
    
    
    
    
    ## z-score 标准化数据
    def normdata(self,dic):
        mea=Statistic.mean(dic)
        sd1=Statistic.sd(dic)
        for key in dict.keys(dic):
            dic[key]=(dic[key]-mea)/sd1
        return dic
    
     # 获取多个k 的频数 时间复杂度 为O((kend-kstart)*mn)  中间count
    def getMulCount(self,seqLis,kstart,kend,sequences): 
         kmerset,dic=self.getSeqKerSet(sequences,kstart)
         countLis,count=self.getSeqCount(seqLis,kstart,dic)
         ## 时间复杂度为O(mn)
         for i in range(len(countLis)):
             countLis[i]=self.normdata(countLis[i])
             
         for k in range(kstart+1,kend+1):
            kmerset,dic=self.getSeqKerSet(sequences,k)
            tmpcountLis,count=self.getSeqCount(seqLis,k,dic)
            for i in range(len(tmpcountLis)):
                tmpcountLis[i]=self.normdata(tmpcountLis[i])
            for i in range(len(countLis)):
                countLis[i]=dict(countLis[i],**(tmpcountLis[i]))
         return countLis
     
         # 获取多个k 的频数 时间复杂度 为O((kend-kstart)*mn) 不标准化, 最后norm
    def getMulCount_nonorm(self,seqLis,kstart,kend,sequences):
        kmerset,dic=self.getSeqKerSet(sequences,kstart)
        countLis,count=self.getSeqCount(seqLis,kstart,dic)
         ## 时间复杂度为O(mn)
#         for i in range(len(countLis)):
#             countLis[i]=self.normdata(countLis[i])
             
        for k in range(kstart+1,kend+1):
            kmerset,dic=self.getSeqKerSet(sequences,k)
            tmpcountLis,count=self.getSeqCount(seqLis,k,dic)
#            for i in range(len(tmpcountLis)):
#                tmpcountLis[i]=self.normdata(tmpcountLis[i])
        for i in range(len(countLis)):
            countLis[i]=dict(countLis[i],**(tmpcountLis[i]))
        return countLis
     
         # 获取多个k 的频数 时间复杂度 为O((kend-kstart)*mn) 不标准化, 最后norm
    def getMulCount_suf(self,seqLis,kstart,kend,sequences):
         kmerset,dic=self.getSeqKerSet(sequences,kstart)
         countLis,count=self.getSeqCount(seqLis,kstart,dic)
         ## 时间复杂度为O(mn)
#         for i in range(len(countLis)):
#             countLis[i]=self.normdata(countLis[i])
             
         for k in range(kstart+1,kend+1):
            kmerset,dic=self.getSeqKerSet(sequences,k)
            tmpcountLis,count=self.getSeqCount(seqLis,k,dic)
#            for i in range(len(tmpcountLis)):
#                tmpcountLis[i]=self.normdata(tmpcountLis[i])
            for i in range(len(countLis)):
                countLis[i]=dict(countLis[i],**(tmpcountLis[i]))
            countLis[i]=self.normdata(countLis[i])
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
    def getD2SCount(self,sequence,sequences,k,r,flag,kmersetdic,kmer_pro):
        ses=[]
        ses.append(sequence)
        # 获得词统计
        lis,count=self.getSeqCount(ses,k,kmersetdic)
        
        ma=mk.Markov()
        prodic={}
        if flag==False:
            prodic=ma.get_Single_kmer_Pro(sequence,sequences,k,r)
#        else:
##            prodic=ma.get_Mul_kmer_Pro(sequence,sequences,k,r)
#            prodic=kmer_pro
        n=len(sequence)
        for key in lis[0].keys():
            if lis[0][key]==0:
                prodic[key]=0
            else:
                prodic[key]=lis[0][key]-n*kmer_pro[key]
        return self.addfloat(prodic)
    
    
    
    ## 统计所有序列的D2scount 
    def getD2SMul_SeqCount(self,sequences,k,r,flag,kmersetdic,kmer_pro):
        Lisdic=[]
        for sequence in sequences:
            tmpdic=self.getD2SCount(sequence,sequences,k,r,flag,kmersetdic,kmer_pro)
            Lisdic.append(tmpdic)
        return Lisdic
    
    
    
    
      ## 统计D2S kmersetdic 表示所有序列所有的kmer集合 O(m*n*k) 先标准化count
    def getD2SCount_pre(self,sequence,sequences,k,r,flag,kmersetdic,kmer_pro):
        ses=[]
        ses.append(sequence)
        # 获得词统计
        lis,count=self.getSeqCount_pre(ses,k,kmersetdic)
 
        ma=mk.Markov()
        prodic={}
        if flag==False:
            prodic=ma.get_Single_kmer_Pro(sequence,sequences,k,r)
#    
        n=len(sequence)
        for key in lis[0].keys():
            if lis[0][key]==0:
                prodic[key]=0
            else:
                prodic[key]=lis[0][key]-n*kmer_pro[key]
        return self.addfloat(prodic)
    
    
    
    
    

                    

    ## 统计D2Star kmersetdic 表示所有序列所有的kmer集合 O(m*n*k)
    def getD2StarCount(self,sequence,sequences,k,r,flag,kmersetdic,kmer_pro):
        n=len(sequence)
        ses=[]
        ses.append(sequence)
        # 获得词统计
        lis,count=self.getSeqCount(ses,k,kmersetdic)
        ma=mk.Markov()
        prodic={}
        if flag==False:
            prodic=ma.get_Single_kmer_Pro(sequence,sequences,k,r)
#        else:
##            prodic=ma.get_Mul_kmer_Pro(sequence,sequences,k,r)
#            prodic=kmer_pro
        n=len(sequence)
        for key in lis[0].keys():
            if lis[0][key]==0 :
                prodic[key]=0
            else:
                prodic[key]=(lis[0][key]-n*kmer_pro[key])/math.sqrt(n*(kmer_pro[key]+np.spacing(1)))
        return self.addfloat(prodic)
    
    ## 统计所有序列的 d2star 
    def getD2Star_Mul_seq_Count(self,sequences,k,r,flag,kmersetdic,kmer_pro):
        dicLis=[]
        for sequence in sequences:
            tmpdic=self.getD2StarCount(sequence,sequences,k,r,flag,kmersetdic,kmer_pro)
            dicLis.append(tmpdic)
        return dicLis
            
    
    
     ## 统计D2Star kmersetdic 表示所有序列所有的kmer集合 O(m*n*k)
    def getD2StarCount_pre(self,sequence,sequences,k,r,flag,kmersetdic,kmer_pro):
        n=len(sequence)
        ses=[]
        ses.append(sequence)
        # 获得词统计
        lis,count=self.getSeqCount_pre(ses,k,kmersetdic)
        ma=mk.Markov()
        prodic={}
        if flag==False:
            prodic=ma.get_Single_kmer_Pro(sequence,sequences,k,r)
#        else:
##            prodic=ma.get_Mul_kmer_Pro(sequence,sequences,k,r)
#            prodic=kmer_pro
        n=len(sequence)
        for key in lis[0].keys():
            if lis[0][key]==0 :
                prodic[key]=0
            else:
                prodic[key]=(lis[0][key]-n*kmer_pro[key])/math.sqrt(n*(kmer_pro[key]+np.spacing(1)))
        return self.addfloat(prodic)
    
    
    ## 统计D2Star kmersetdic 表示所有序列所有的kmer集合 O(m*n*k)
    def getD2StarCount_suf(self,sequence,sequences,k,r,flag,kmersetdic,kmer_pro):
        n=len(sequence)
        ses=[]
        ses.append(sequence)
        # 获得词统计
        lis,count=self.getSeqCount(ses,k,kmersetdic)
        ma=mk.Markov()
        prodic={}
        if flag==False:
            prodic=ma.get_Single_kmer_Pro(sequence,sequences,k,r)
#        else:
##            prodic=ma.get_Mul_kmer_Pro(sequence,sequences,k,r)
#            prodic=kmer_pro
        n=len(sequence)
        for key in lis[0].keys():
            if lis[0][key]==0 :
                prodic[key]=0
            else:
                prodic[key]=(lis[0][key]-n*kmer_pro[key])/math.sqrt(n*(kmer_pro[key]+np.spacing(1)))
        return self.addfloat(prodic)
    
    
    
    
    ## 统计多个k值 D2S kmersetdic 表示所有序列所有的kmer集合 O((kend-kstart)*m*n*k) 
    # 单条标准化后面的特征
    def getD2SMulCount(self,sequence,sequences,kstart,kend,r,flag,kmer_pro):
        kmerset,kmersetdic=self.getSeqKerSet(sequences,kstart)
        prodic =self.getD2SCount(sequence,sequences,kstart,r,flag,kmersetdic,kmer_pro)
        prodic =self.normdata(prodic)
        for k in range(kstart+1,kend+1):
            kmerset,kmersetdic=self.getSeqKerSet(sequences,k)
            tmpprodic =self.getD2SCount(sequence,sequences,k,r,flag,kmersetdic,kmer_pro)
            tmpprodic =self.normdata(tmpprodic)
            prodic=dict(prodic,**tmpprodic)
        return prodic
    
    
    ## 中标准化count
    def getD2SMulCount_mid(self,sequences,kstart,kend,r,flag,kmer_pro):
        res=[]
        ## 每一条序列的多个k值
        for sequence in sequences:
            kmerset,kmersetdic=self.getSeqKerSet(sequences,kstart)
            prodic =self.getD2SCount(sequence,sequences,kstart,r,flag,kmersetdic,kmer_pro)
            prodic =self.normdata(prodic)
            for k in range(kstart+1,kend+1):
                kmerset,kmersetdic=self.getSeqKerSet(sequences,k)
                tmpprodic =self.getD2SCount(sequence,sequences,k,r,flag,kmersetdic,kmer_pro)
                tmpprodic =self.normdata(tmpprodic)
                prodic=dict(prodic,**tmpprodic)
            res.append(prodic)
        return res
    
    ## 不norm
    def getD2SMulCount_nonorm(self,sequences,kstart,kend,r,flag,kmer_pro):
        res=[]
        ## 每一条序列的多个k值
        for sequence in sequences:
            kmerset,kmersetdic=self.getSeqKerSet(sequences,kstart)
            prodic =self.getD2SCount(sequence,sequences,kstart,r,flag,kmersetdic,kmer_pro)
#            prodic =self.normdata(prodic)
            for k in range(kstart+1,kend+1):
                kmerset,kmersetdic=self.getSeqKerSet(sequences,k)
                tmpprodic =self.getD2SCount(sequence,sequences,k,r,flag,kmersetdic,kmer_pro)
#                tmpprodic =self.normdata(tmpprodic)
                prodic=dict(prodic,**tmpprodic)
#            prodic=self.normdata(prodic)
            res.append(prodic)
        return res
    
      ## 后norm
    def getD2SMulCount_suf(self,sequences,kstart,kend,r,flag,kmer_pro):
        res=[]
        ## 每一条序列的多个k值
        for sequence in sequences:
            kmerset,kmersetdic=self.getSeqKerSet(sequences,kstart)
            prodic =self.getD2SCount(sequence,sequences,kstart,r,flag,kmersetdic,kmer_pro)
#            prodic =self.normdata(prodic)
            for k in range(kstart+1,kend+1):
                kmerset,kmersetdic=self.getSeqKerSet(sequences,k)
                tmpprodic =self.getD2SCount(sequence,sequences,k,r,flag,kmersetdic,kmer_pro)
#                tmpprodic =self.normdata(tmpprodic)
                prodic=dict(prodic,**tmpprodic)
            prodic=self.normdata(prodic)
            res.append(prodic)
        return res
    
    
    ## 先标准化count
    def getD2SMulCount_pre(self,sequences,kstart,kend,r,flag,kmer_pro):
        res=[]
        for sequence in sequences:
            kmerset,kmersetdic=self.getSeqKerSet(sequences,kstart)
            prodic =self.getD2SCount_pre(sequence,sequences,kstart,r,flag,kmersetdic,kmer_pro)
#            prodic =self.normdata(prodic)    去掉一句 标准化
            for k in range(kstart+1,kend+1):
                kmerset,kmersetdic=self.getSeqKerSet(sequences,k)
                tmpprodic =self.getD2SCount_pre(sequence,sequences,k,r,flag,kmersetdic,kmer_pro)
#                tmpprodic =self.normdata(tmpprodic)    去掉一句 标准化
                prodic=dict(prodic,**tmpprodic)
            res.append(prodic)
        return res
                   


     ## 统计多个k值 D2Star kmersetdic 表示所有序列所有的kmer集合 O((kend-kstart)*m*n*k)
     # 单条
    def getD2StarMulCount(self,sequence,sequences,kstart,kend,r,flag,kmer_pro):
        kmerset,kmersetdic=self.getSeqKerSet(sequences,kstart)
        prodic =self.getD2StarCount(sequence,sequences,kstart,r,flag,kmersetdic,kmer_pro)
        prodic =self.normdata(prodic)
        for k in range(kstart+1,kend+1):
            kmerset,kmersetdic=self.getSeqKerSet(sequences,k)
            tmpprodic =self.getD2StarCount(sequence,sequences,k,r,flag,kmersetdic,kmer_pro)
            tmpprodic =self.normdata(tmpprodic)
            prodic=dict(prodic,**tmpprodic)
        return prodic
    
    
    
     ## 统计多个k值 D2Star kmersetdic 表示所有序列所有的kmer集合 O((kend-kstart)*m*n*k) 先标准化count
     # 单条序列
    def getD2StarMulCount_preSingle(self,sequence,sequences,kstart,kend,r,flag,kmer_pro):
        kmerset,kmersetdic=self.getSeqKerSet(sequences,kstart)
        prodic =self.getD2StarCount_pre(sequence,sequences,kstart,r,flag,kmersetdic,kmer_pro)
#        prodic =self.normdata(prodic)
        for k in range(kstart+1,kend+1):
            kmerset,kmersetdic=self.getSeqKerSet(sequences,k)
            tmpprodic =self.getD2StarCount_pre(sequence,sequences,k,r,flag,kmersetdic,kmer_pro)
#            tmpprodic =self.normdata(tmpprodic)
            prodic=dict(prodic,**tmpprodic)
        return prodic
    
  
     # 多条序列 中mid
    def getD2StarMulCount_mid(self,sequences,kstart,kend,r,flag,kmer_pro):
        res=[]
        for sequence in sequences:
            kmerset,kmersetdic=self.getSeqKerSet(sequences,kstart)
            prodic =self.getD2StarCount(sequence,sequences,kstart,r,flag,kmersetdic,kmer_pro)
            prodic =self.normdata(prodic)
            for k in range(kstart+1,kend+1):
                kmerset,kmersetdic=self.getSeqKerSet(sequences,k)
                tmpprodic =self.getD2StarCount(sequence,sequences,k,r,flag,kmersetdic,kmer_pro)
                tmpprodic =self.normdata(tmpprodic)
                prodic=dict(prodic,**tmpprodic)
            res.append(prodic)
        return res
    
    ## 单个markov模型，无norm
    def getD2StarMulCount_nonorm(self,sequences,kstart,kend,r,flag,kmer_pro):
        res=[]
        for sequence in sequences:
            kmerset,kmersetdic=self.getSeqKerSet(sequences,kstart)
            ## 可以更改为后置getD2StarCount_suf，前置getD2StarCount_pre，中间getD2StarCount
            prodic =self.getD2StarCount_suf(sequence,sequences,kstart,r,flag,kmersetdic,kmer_pro)
#            prodic =self.normdata(prodic)
            for k in range(kstart+1,kend+1):
                kmerset,kmersetdic=self.getSeqKerSet(sequences,k)
                ## 可以更改为后置getD2StarCount_suf，前置getD2StarCount_pre，中间getD2StarCount
                tmpprodic =self.getD2StarCount(sequence,sequences,k,r,flag,kmersetdic,kmer_pro)
#                tmpprodic =self.normdata(tmpprodic)
                prodic=dict(prodic,**tmpprodic)
#            prodic=self.normdata(prodic)
            res.append(prodic)
        return res
    
     ## 单个markov模型，无norm
    def getD2StarMulCount_pre(self,sequences,kstart,kend,r,flag,kmer_pro):
        res=[]
        for sequence in sequences:
            kmerset,kmersetdic=self.getSeqKerSet(sequences,kstart)
            ## 可以更改为后置getD2StarCount_suf，前置getD2StarCount_pre，中间getD2StarCount
            prodic =self.getD2StarCount_pre(sequence,sequences,kstart,r,flag,kmersetdic,kmer_pro)
#            prodic =self.normdata(prodic)
            for k in range(kstart+1,kend+1):
                kmerset,kmersetdic=self.getSeqKerSet(sequences,k)
                ## 可以更改为后置getD2StarCount_suf，前置getD2StarCount_pre，中间getD2StarCount
                tmpprodic =self.getD2StarCount_pre(sequence,sequences,k,r,flag,kmersetdic,kmer_pro)
#                tmpprodic =self.normdata(tmpprodic)
                prodic=dict(prodic,**tmpprodic)
#            prodic=self.normdata(prodic)
            res.append(prodic)
        return res
    
        ## 单个markov模型，后norm
    def getD2StarMulCount_suf(self,sequences,kstart,kend,r,flag,kmer_pro):
        res=[]
        for sequence in sequences:
            kmerset,kmersetdic=self.getSeqKerSet(sequences,kstart)
            ## 可以更改为后置getD2StarCount_suf，前置getD2StarCount_pre，中间getD2StarCount
            prodic =self.getD2StarCount_suf(sequence,sequences,kstart,r,flag,kmersetdic,kmer_pro)
#            prodic =self.normdata(prodic)
            for k in range(kstart+1,kend+1):
                kmerset,kmersetdic=self.getSeqKerSet(sequences,k)
                ## 可以更改为后置getD2StarCount_suf，前置getD2StarCount_pre，中间getD2StarCount
                tmpprodic =self.getD2StarCount_suf(sequence,sequences,k,r,flag,kmersetdic,kmer_pro)
#                tmpprodic =self.normdata(tmpprodic)
                prodic=dict(prodic,**tmpprodic)
            prodic=self.normdata(prodic)
            res.append(prodic)
        return res
    
    
    # 获取多个k 的频率
    def getMulFreq_mid(self,sequences,kstart,kend):
         kmerset,dic=self.getSeqKerSet(sequences,kstart)
         freqLis,freq=self.getSeqfreq(sequences,kstart,dic)
         for i in range(len(freqLis)):
             freqLis[i]=self.normdata(freqLis[i])
         for k in range(kstart+1,kend+1):
            kmerset,dic=self.getSeqKerSet(sequences,k)
            tmpfreqLis,freq=self.getSeqfreq(sequences,k,dic)
            for i in range(len(tmpfreqLis)):
                tmpfreqLis[i]=self.normdata(tmpfreqLis[i])
            for i in range(len(freqLis)):
                freqLis[i]=dict(freqLis[i],**(tmpfreqLis[i]))
         return freqLis
     
        # 获取多个k 的频率
    def getMulFreq_nonorm(self,sequences,kstart,kend):
         kmerset,dic=self.getSeqKerSet(sequences,kstart)
         freqLis,freq=self.getSeqfreq(sequences,kstart,dic)
         for k in range(kstart+1,kend+1):
            kmerset,dic=self.getSeqKerSet(sequences,k)
            tmpfreqLis,freq=self.getSeqfreq(sequences,k,dic)
            for i in range(len(freqLis)):
                freqLis[i]=dict(freqLis[i],**(tmpfreqLis[i]))
         return freqLis
     
        
     #后置
    def getMulFreq_suf(self,sequences,kstart,kend):
         kmerset,dic=self.getSeqKerSet(sequences,kstart)
         freqLis,freq=self.getSeqfreq(sequences,kstart,dic)
         for k in range(kstart+1,kend+1):
            kmerset,dic=self.getSeqKerSet(sequences,k)
            tmpfreqLis,freq=self.getSeqfreq(sequences,k,dic)
            for i in range(len(freqLis)):
                freqLis[i]=dict(freqLis[i],**(tmpfreqLis[i]))
        ## 后置
         for i in range(len(freqLis)):
            freqLis[i]=self.normdata(freqLis[i])
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
        suf=np.spacing(1)
        for key in resultdic:
            suf=resultdic[key]+suf
        for key in resultdic:
            resultdic[key]=resultdic[key]/suf
        return resultdic
#    def getWeightFromSeq(self,sequences):
        
    
    ## O(M*n)  方差
    def getWeight2(self,freqLis):
        # 获得字典
        resultdic = dict.copy(freqLis[0]) 
        for key in resultdic:
            mea=0.0
            va=0.0
            ##求平均值
            for j in range(len(freqLis)):
                mea=mea+freqLis[j][key]
            mea=mea/len(freqLis)
            
            for j in range(len(freqLis)):
                va=(freqLis[j][key]-mea)**2+va
            va=va/len(freqLis)
            sd=math.sqrt(va)
            resultdic[key]=sd
        # 计算权重    
        suf=0.0
        for key in resultdic:
            suf=resultdic[key]+suf
        for key in resultdic:
            resultdic[key]=resultdic[key]/suf
        return resultdic
    
    
    ## 获取多个k的权重  freqLis表示频率矩阵 时间复杂度O((kend-kstart)*n*m^2)
    ## 中norm
    def getMulWeight_mid(self,sequences,kstart,kend):
        #获取第一个频率合集合字典 
        freqLis=self.getMulFreq_mid(sequences,kstart,kend)
        resultdic=self.getWeight(freqLis)
        return resultdic
    
    
    
       ## 获取多个k的权重  freqLis表示频率矩阵 时间复杂度O((kend-kstart)*n*m^2)
       ## 无norm
    def getMulWeight_nonorm(self,sequences,kstart,kend):
        #获取第一个频率合集合字典 
        freqLis=self.getMulFreq_nonorm(sequences,kstart,kend)
        resultdic=self.getWeight(freqLis)
        return resultdic
    
       ## 获取多个k的权重  freqLis表示频率矩阵 时间复杂度O((kend-kstart)*n*m^2)
       ## 无norm
    def getMulWeight_suf(self,sequences,kstart,kend):
        #获取第一个频率合集合字典 
        freqLis=self.getMulFreq_suf(sequences,kstart,kend)
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