# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 10:50:07 2018

@author: Yzi

"""
## 多个k值的 以全部数据集求权重和特征
import ReadData
import Sequence
import Similarity
import time
import markov

## 返回列表的第2个 元素，即相似度
def getSim(Lis):
    return Lis[1]

if __name__=="__main__":
    rd=ReadData.ReadData()
    datasets,pos,neg=rd.getData2("fly_blastoderm")
#    datasets,pos,neg=rd.getData2("fly_pns")
    Sim=Similarity.Similarity()
    kstart=2
    kend=6
    sizePo=len(pos)
    top=int(sizePo*(sizePo-1))
    flag=True
#    top=3321
    sq=Sequence.Sequence()
    
    start = time.process_time()
    posweight = sq.getMulWeight_mid(pos,kstart,kend)
    negweight = sq.getMulWeight_mid(neg,kstart,kend)
    end = time.process_time()
    print(len(posweight))
    print(len(negweight))
    print("权重计算时间：",(end-start))

    print("--------------------------MulD2-----------------------")
    print("kstart=",kstart," kend=",kend)
    start = time.process_time()
    ## 计算MulD2特征
    posD2Lis=sq.getMulCount(pos,kstart,kend,pos) 
    negD2Lis=sq.getMulCount(neg,kstart,kend,neg)
    
    ## 构成基因对,并设置标识
    posPair=[]
    negPair=[]
    for i in range(len(posD2Lis)):
        for j in range(i+1,len(posD2Lis)):
            posSim=Sim.getMulD2WeightSim2(posD2Lis[i],posD2Lis[j],posweight)
            negSim=Sim.getMulD2WeightSim2(negD2Lis[i],negD2Lis[j],negweight)
            posPair.append(["+",posSim])
            negPair.append(["-",negSim])
    Pair=posPair+negPair
    pSim=sorted(Pair,key=getSim,reverse=True)
    corrCnt=[]
    corrFre=[]
    count=0
    for i in range(top):
        if pSim[i][0]=="+":
            count=count+1
        if (i+1)%10==0 or i+1==top:
            corrCnt.append(count)
            corrFre.append(count/(i+1))
    print("预测准确个数：")
    print(corrCnt)
    print("预测准确率：")
    print(corrFre)
    end = time.process_time()
    print("程序运行时间：",(end-start))


    for r in range(0,3):
        start=time.process_time()
        ma=markov.Markov()
        kmer_propos=ma.get_Mulk_Mul_kmer_Pro(pos,kstart,kend,r)
        kmer_proneg=ma.get_Mulk_Mul_kmer_Pro(neg,kstart,kend,r)
        end=time.process_time()
        print("计算马尔可夫概率时间：",(end-start))
        print("--------------------------MulD2s-----------------------")
        print("kstart=",kstart," kend=",kend," r=",r)
        start = time.process_time()
        ## 计算MulD2特征
        start1=time.process_time()
        posD2Lis=sq.getD2SMulCount_mid(pos,kstart,kend,r,flag,kmer_propos) 
        negD2Lis=sq.getD2SMulCount_mid(neg,kstart,kend,r,flag,kmer_proneg)
        end1=time.process_time()
        print("特征计算时间",(end1-start1))
        ## 构成基因对,并设置标识
        posPair=[]
        negPair=[]
        for i in range(len(posD2Lis)):
            for j in range(i+1,len(posD2Lis)):
                posSim=Sim.getMulD2sWeightSim2(posD2Lis[i],posD2Lis[j],posweight)
                negSim=Sim.getMulD2sWeightSim2(negD2Lis[i],negD2Lis[j],negweight)
                posPair.append(["+",posSim])
                negPair.append(["-",negSim])
        Pair=posPair+negPair
        pSim=sorted(Pair,key=getSim,reverse=True)
        corrCnt=[]
        corrFre=[]
        count=0
        for i in range(top):
            if pSim[i][0]=="+":
                count=count+1
            if (i+1)%10==0 or i+1==top:
                corrCnt.append(count)
                corrFre.append(count/(i+1))
        print("预测准确个数：")
        print(corrCnt)
        print("预测准确率：")
        print(corrFre)
        end = time.process_time()
        print("程序运行时间：",(end-start))
        
        
        print("--------------------------MulD2star-----------------------")
        print("kstart=",kstart," kend=",kend," r=",r)
        start = time.process_time()
        ## 计算MulD2特征
        posD2Lis=sq.getD2StarMulCount_mid(pos,kstart,kend,r,flag,kmer_propos) 
        negD2Lis=sq.getD2StarMulCount_mid(neg,kstart,kend,r,flag,kmer_proneg)
        ## 构成基因对,并设置标识
        posPair=[]
        negPair=[]
        for i in range(len(posD2Lis)):
            for j in range(i+1,len(posD2Lis)):
                posSim=Sim.getMulD2starWeightSim2(posD2Lis[i],posD2Lis[j],posweight)
                negSim=Sim.getMulD2starWeightSim2(negD2Lis[i],negD2Lis[j],negweight)
                posPair.append(["+",posSim])
                negPair.append(["-",negSim])
        Pair=posPair+negPair
        pSim=sorted(Pair,key=getSim,reverse=True)
        corrCnt=[]
        corrFre=[]
        count=0
        for i in range(top):
            if pSim[i][0]=="+":
                count=count+1
            if (i+1)%10==0 or i+1==top:
                corrCnt.append(count)
                corrFre.append(count/(i+1))
        print("预测准确个数：")
        print(corrCnt)
        print("预测准确率：")
        print(corrFre)
        end = time.process_time()
        print("程序运行时间：",(end-start))
     
            
        