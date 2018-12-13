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
    kend=6
    sq=Sequence.Sequence()
    Sim=Similarity.Similarity()
    dis=Distance.Distance()
    ## 求权重 尝试单条序列的权重  这里为了使用后面的philip软件，求的是距离
    matrix=[]
    for i in range(len(datasets)):
        tmpRes=[sqName[i]]
        for j in range(len(datasets)):
            ## 如果两个序列一样，他们的相似度
            if(i==j):
                tmpRes.append(round(0.0,6))
            else:
                pairs=[]
                pairs.append(datasets[i])
                pairs.append(datasets[j])
                # 获得特征
                D2Lis=sq.getMulCount(pairs,kstart,kend,pairs)
                # 计算权重
                weight = sq.getMulWeight_suf(pairs,kstart,kend) 
                dist=dis.getMulD2Weight2(D2Lis[0],D2Lis[1],weight) 
                tmpRes.append(round(dist,6))
        matrix.append(tmpRes)
    ## 写入文件
    f=open("resultE3/D2",'w')
    f.write(str(numSeq))
    f.write("\n")
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            f.write(str(matrix[i][j]))
            f.write("\t")
        f.write("\n")
    f.close()
    print("end")
