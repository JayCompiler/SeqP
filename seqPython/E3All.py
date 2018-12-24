# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 23:02:13 2018

@author: Yzi
"""
## 以全部序列求权重，全部为

import ReadData
import Sequence
import Similarity
import Distance
import time
import markov
if __name__=="__main__":
    rd=ReadData.ReadData()
    ## 读取数据
    name="primates" 
    datasets,sqName=rd.getData3(name)
    ## 参与比较的序列个数
    numSeq=len(sqName)
    kstart=2
    kend=10
    flag=True
    sq=Sequence.Sequence()
    Sim=Similarity.Similarity()
    dis=Distance.Distance()
    ### 
    start = time.process_time()
    weight = sq.getMulWeight_mid(datasets,kstart,kend)
    
     
    end = time.process_time()
    print("权重计算时间：",(end-start))
       ## 计算MulD2特征---------------------1--------------
#    D2Lis=sq.getMulCount(datasets,kstart,kend,datasets) 
#   
#    ## 求权重 尝试单条序列的权重  这里为了使用后面的philip软件，求的是距离
#    matrix=[]
#    print(kend)
#    start = time.process_time()
#    for i in range(len(D2Lis)):
#        tmpRes=[sqName[i]]
#        for j in range(len(D2Lis)):
#             ## 如果两个序列一样，他们的相似度
#            if(i==j):
#                tmpRes.append(round(0.0,6))
#            else:
#                dist=dis.getMulD2Weight2(D2Lis[i],D2Lis[j],weight)
#                tmpRes.append(round(dist,6))
#        matrix.append(tmpRes)
#    ## 写入文件
#    f=open("resultE3/D2All_14",'w')
#    f.write(str(numSeq))
#    f.write("\n")
#    for i in range(len(matrix)):
#        for j in range(len(matrix[0])):
#            f.write(str(matrix[i][j]))
#            f.write("\t")
#        f.write("\n")
#    f.close()
#    end = time.process_time()
#    print("程序运行时间：",(end-start))




    for r in range(1,2):
        ## 求得 kmer概率
        start=time.process_time()
        ma=markov.Markov()
        kmer_pro=ma.get_Mulk_Mul_kmer_Pro(datasets,kstart,kend,r)
        end=time.process_time()
        
        print("计算马尔可夫概率时间：",(end-start))
        ###--------------------------------2--------------------------
        print("--------------------------MulD2s-----------------------")
        print("kstart=",kstart," kend=",kend," r=",r)
        start = time.process_time()
        ## 计算MulD2特征
        start1=time.process_time()
        D2Lis=sq.getD2SMulCount_mid(datasets,kstart,kend,r,flag,kmer_pro) 
        end1=time.process_time()
        print("特征计算时间",(end1-start1))
        matrix=[]
        print(kend)
        for i in range(len(D2Lis)):
            tmpRes=[sqName[i]]
            for j in range(len(D2Lis)):
                 ## 如果两个序列一样，他们的相似度
                if(i==j):
                    tmpRes.append(round(0.0,6))
                else:
                    dist=dis.getMulD2SWeight2(D2Lis[i],D2Lis[j],weight)
                    tmpRes.append(round(dist,6))
            matrix.append(tmpRes)
         ## 写入文件
        f=open("resultE3/D2sAll_101",'w')
        f.write(str(numSeq))
        f.write("\n")
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                f.write(str(matrix[i][j]))
                f.write("\t")
            f.write("\n")
        f.close()        
        end = time.process_time()
        print("程序运行时间：",(end-start))

       ###---------------------------------------------3----------------------
        print("--------------------------MulD2star-----------------------")
        print("kstart=",kstart," kend=",kend," r=",r)
        start = time.process_time()
        ## 计算MulD2特征
        D2Lis=sq.getD2StarMulCount_mid(datasets,kstart,kend,r,flag,kmer_pro) 
        ## 构成基因对,并设置标识
        matrix=[]
        print(kend)
        for i in range(len(D2Lis)):
            tmpRes=[sqName[i]]
            for j in range(len(D2Lis)):
                 ## 如果两个序列一样，他们的相似度
                if(i==j):
                    tmpRes.append(round(0.0,6))
                else:
                    dist=dis.getMulD2StarWeight2(D2Lis[i],D2Lis[j],weight)
                    tmpRes.append(round(dist,6))
            matrix.append(tmpRes)
        ## 写入文件
        f=open("resultE3/D2starAll_101",'w')
        f.write(str(numSeq))
        f.write("\n")
        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                f.write(str(matrix[i][j]))
                f.write("\t")
            f.write("\n")
        f.close()          
        end = time.process_time()
        print("程序运行时间：",(end-start))