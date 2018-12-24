# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 12:39:17 2018

@author: Yzi
"""

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
    kend=8
    flag=True
    sq=Sequence.Sequence()
    Sim=Similarity.Similarity()
    dis=Distance.Distance()
    ## 求权重 尝试单条序列的权重  这里为了使用后面的philip软件，求的是距离
    ##--------------------------------------------------------1-----------------------
#    matrix=[]
#    print(kend)
#    start = time.process_time()
#    for i in range(len(datasets)):
#        tmpRes=[sqName[i]]
#        for j in range(len(datasets)):
#            ## 如果两个序列一样，他们的相似度
#            if(i==j):
#                tmpRes.append(round(0.0,6))
#            else:
#                pairs=[]
#                pairs.append(datasets[i])
#                pairs.append(datasets[j])
#                # 获得特征
#                D2Lis=sq.getMulCount(pairs,kstart,kend,pairs)
#                # 计算权重
#                weight = sq.getMulWeight_suf(pairs,kstart,kend) 
#                dist=dis.getMulD2Weight2(D2Lis[0],D2Lis[1],weight) 
#                tmpRes.append(round(dist,6))
#        matrix.append(tmpRes)
#    ## 写入文件
#    f=open("resultE3/D2_14",'w')
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
    
    
    ## 计算 概率
    start = time.process_time()
    weight = sq.getMulWeight_mid(datasets,kstart,kend)
    end = time.process_time()
    print("计算markov概率时间：",end-start)
    for r in range(0,1):
 
        ###--------------------------------2---------------------------------------
#        print("--------------------------MulD2s-----------------------")
#        print("kstart=",kstart," kend=",kend," r=",r)
#     
#        matrix=[]
#        print(kend)
#        start = time.process_time()
#
#        for i in range(len(datasets)):
#            tmpRes=[sqName[i]]
#            for j in range(len(datasets)):
#                 ## 如果两个序列一样，他们的相似度
#                if(i==j):
#                    tmpRes.append(round(0.0,6))
#                else:
#                    SeqLis=[]
#                    SeqLis.append(datasets[i])
#                    SeqLis.append(datasets[j])
#            
#                    # 计算概率
#
#                    ma=markov.Markov()
#                    kmer_pro=ma.get_Mulk_Mul_kmer_Pro(SeqLis,kstart,kend,r)
##               
#                    
#                    # 获得特征
#                    D2Lis=sq.getD2SMulCount_mid(SeqLis,kstart,kend,r,True,kmer_pro)
#                    ## 计算权重
##                    weight = sq.getMulWeight_mid(SeqLis,kstart,kend) 
#                    
#                    dist=dis.getMulD2SWeight2(D2Lis[0],D2Lis[1],weight)                               
#                    tmpRes.append(round(dist,6))
#            matrix.append(tmpRes)
#            
#         ## 写入文件
#        f=open("resultE3/D2sSingle_8_mid_globalweight",'w')
#        f.write(str(numSeq))
#        f.write("\n")
#        for i in range(len(matrix)):
#            for j in range(len(matrix[0])):
#                f.write(str(matrix[i][j]))
#                f.write("\t")
#            f.write("\n")
#        f.close()        
#        end = time.process_time()
#        print("程序运行时间：",(end-start))

       ###---------------------------------------------3----------------------
        print("--------------------------MulD2star-----------------------")
        print("kstart=",kstart," kend=",kend," r=",r)
        start = time.process_time()
        matrix=[]
        print(kend)
        for i in range(len(datasets)):
            tmpRes=[sqName[i]]
            for j in range(len(datasets)):
                 ## 如果两个序列一样，他们的相似度
                if(i==j):
                    tmpRes.append(round(0.0,6))
                else:
                    SeqLis=[]
                    SeqLis.append(datasets[i])
                    SeqLis.append(datasets[j])
            
                    # 计算概率
#                    start=time.process_time()
                    ma=markov.Markov()
                    kmer_pro=ma.get_Mulk_Mul_kmer_Pro(SeqLis,kstart,kend,r)
#                    end=time.process_time()
                    
                    # 获得特征
                    D2Lis=sq.getD2StarMulCount_suf(SeqLis,kstart,kend,r,True,kmer_pro)
                    ## 计算权重
#                    weight = sq.getMulWeight_suf(SeqLis,kstart,kend) 
                    
                    dist=dis.getMulD2StarWeight2(D2Lis[0],D2Lis[1],weight)                               
                    tmpRes.append(round(dist,6))
            matrix.append(tmpRes)
        ## 写入文件
        f=open("resultE3/D2starSingle_8_globalweight",'w')
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
