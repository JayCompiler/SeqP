# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 10:09:25 2018

@author: Yzi
"""
import Sequence
class Markov:
    # first-order markov
    dictOneOrder = {"A":0,"T":1,"C":2,"G":3}
    # second-order markov
    dictTwoOrder = {"AA": 0,"AT":1,"AC": 2,"AG":3,"TA":4,"TT":5,"TC":6,
                  "TG":7,"CA":8,"CT":9,"CC":10,"CG":11,"GA":12,"GT":13,
                  "GC":14,"GG":15}
    
    
    ## 单条序列的初识概率 序列长度n：tc=O(n)
    def init_Single_pro(self,sequence): 
        lis=["A","T","C","G"]
        prodic=dict.fromkeys(lis,0)
        for i in range(len(sequence)):
            prodic[sequence[i]]=prodic[sequence[i]]+1
        length=len(sequence)
        for key in prodic.keys():
            prodic[key]=prodic[key]/length
        return prodic
    
    ## 多条序列的初始矩阵 序列长度 n,m条序列，tc=O(mn)
    def init_MUl_pro(self,sequences):
        lis=["A","T","C","G"]
        prodic=dict.fromkeys(lis,0)
        for sequence in sequences:
            for i in range(len(sequence)):
                prodic[sequence[i]]=prodic[sequence[i]]+1
        length=0
        for sequence in sequences:
                length=len(sequence)+length
        for key in prodic.keys():
            prodic[key]=prodic[key]/length
        return prodic
        
    ## 单条序列的状态转移矩阵概率 序列长度为n tc=O(n),r>0
    def trans_SingleSeq_Matrix(self,sequence,sequences,r):
        if r<=0:
            print("r<=0，无转移矩阵，请重新输入")
            return
        lis=[]
        lis.append(sequence)
        # 获取前缀个数
        Sq=Sequence.Sequence()
        kmerset,dic=Sq.getSeqKerSet(sequences,r)
        # 统计个数
        dic0,count0=Sq.getSeqCount(lis,r,dic)
        # 去掉最后一位
        dic0[0][sequence[-r:]]=dic0[0][sequence[-r:]]-1
        # 获取后缀个数
        kmerset1,dic1=Sq.getSeqKerSet(sequences,r+1)
        # 统计个数
        dic2,count2=Sq.getSeqCount(lis,r+1,dic1)
        
        resultdic=dict.copy(dic2[0])
        for key in dic2[0].keys():
            if dic0[0][key[0:-1]]==0:
                resultdic[key]=0
            else:
                resultdic[key]=dic2[0][key]/dic0[0][key[0:-1]]
        return resultdic
    
    
    ## 多条序列的状态转移矩阵概率 序列长度为n,序列个数为m tc=O(nm)
    def trans_MulSeq_Matrix(self,sequences,r):
        if r<=0:
            print("r<=0，无转移矩阵，请重新输入")
            return
        Sq=Sequence.Sequence()
        # trim 化序列 
        trimSeq=Sq.trimSequence(sequences)
        # 获取前缀：
        kmerset,dic=Sq.getSeqKerSet(trimSeq,r)
        predic =Sq.getkmerCount(trimSeq,dic,r)
        #获取后缀
        kmerset1,dic1=Sq.getSeqKerSet(sequences,r+1)
        suffdic =Sq.getkmerCount(sequences,dic1,r+1)
        dic={}
        for key in suffdic.keys():
            dic[key]=suffdic[key]/predic[key[0:-1]]
        return dic
    
    ##0，1，2阶马尔柯夫模型的k-mer概率 ，以单条序列初始状态与状态转移矩阵估计概率 
    ## 序列长度n，kmer 长度k，tc=O(nk)
    def get_Single_kmer_Pro(self,sequence,sequences,k,r):
        if r>=k:
            r=0
        Sq=Sequence.Sequence()
        # 单条序列的kmerset,dic
        kmerSet,dic=Sq.getSeqKerSet(sequences,k)
        resultdic=dict.copy(dic)

        # 初始概率：
        initProdic=self.init_Single_pro(sequence)
        if r==0:
            for kmer in dict.keys(resultdic):
                pro=1
                for i in range (len(kmer)):
                    pro=initProdic[kmer[i]]*pro
                resultdic[kmer]=pro 
        elif r<0:
            print("r的值设定有误，不能小于0")
        else:
            # 状态转移
            transdic=self.trans_SingleSeq_Matrix(sequence,sequences,r)
            for kmer in dict.keys(resultdic):
            # 初始概率
                pro=1
                for i in range(r):
                    pro=initProdic[kmer[i]]*pro
            # kmer概率
                for loc in range(len(kmer)-r):
                    pro=pro*transdic[kmer[loc:loc+r+1]]
                resultdic[kmer]=pro    
        return Sq.addfloat(resultdic)
    
    def get_MulK_Single_kmer_pro(self,sequence,sequences,kstart,kend,r):
        resultdic = self.get_Single_kmer_Pro(sequence,sequences,kstart,r)
        for k in range(kstart+1,kend+1):
            tmpresultdic=self.get_Single_kmer_Pro(sequence,sequences,k,r)
            resultdic=dict(resultdic,tmpresultdic)
        return resultdic
    
    
    
    ##0，1，2阶马尔柯夫模型的k-mer概率 ，以多条序列初始状态与状态转移矩阵估计概率 
    ## 序列长度n，kmer 长度k，tc=O(nk)
    def get_Mul_kmer_Pro(self,sequence,sequences,k,r):
        if r>=k:
            r=0
        Sq=Sequence.Sequence()
        # 单条序列的kmerset,dic
#        kmerSet,dic=Sq.getSingleSeqKerSet(sequence,k)
        kmers,dic1=Sq.getSeqKerSet(sequences,k)
        resultdic=dict.copy(dic1)
#        print("sad",resultdic)
        # 初始概率：
        initProdic=self.init_MUl_pro(sequences)
        if r==0:
            for kmer in dict.keys(dic1):
                pro=1
                for i in range (len(kmer)):
                    pro=initProdic[kmer[i]]*pro
                resultdic[kmer]=pro 
        elif r<0:
            print("r的值设定有误，不能小于0")
        else:
            # 状态转移
            transdic=self.trans_MulSeq_Matrix(sequences,r)
            for kmer in dict.keys(dic1):
            # 初始概率
                pro=1
                for i in range(r):
                    pro=initProdic[kmer[i]]*pro
            # kmer概率
                for loc in range(len(kmer)-r):
                    pro=pro*transdic[kmer[loc:loc+r+1]]
                resultdic[kmer]=pro 
        return Sq.addfloat(resultdic)
    
    
    def get_Mulk_Mul_kmer_Pro(self,sequence,sequences,kstart,kend,r):
        resultdic = self.get_Mul_kmer_Pro(sequence,sequences,kstart,r)
        for k in range(kstart+1,kend+1):
            tmpresultdic=self.get_Mul_kmer_Pro(sequence,sequences,k,r)
            resultdic=dict(resultdic,tmpresultdic)
        return resultdic
        


if __name__ =="__main__":
     sequence=["ATCGGCGC","TTT","TCGGA"]
     MA=Markov()
     dic =MA.init_Single_pro(sequence[0])
     print(dic)
     lis=[]
     lis.append(sequence[0])
     lis.append(sequence[1])
     dic1 =MA.trans_SingleSeq_Matrix(sequence[0],lis,2)
     print(dic1)
#     Sq=Sequence.Sequence()
#     kmer,dic=Sq.getSeqKerSet(sequence,3)
#     pp=MA.get_Single_kmer_Pro(sequence[0],3,2,dic)
#     print(pp)
#     print("----------------")
#     
#     dic =MA.init_MUl_pro(sequence)
#     print(dic)
#     dic1 =MA.trans_MulSeq_Matrix(sequence,2)
#     print(dic1)
#     pp=MA.get_Mul_kmer_Pro(sequence[0],sequence,3,2)
#     print(pp)