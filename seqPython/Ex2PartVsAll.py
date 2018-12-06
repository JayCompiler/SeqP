# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 13:32:05 2018

@author: Yzi
"""


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
    flag=True
    top=150
    sq=Sequence.Sequence()
    posHalfLen=int(len(pos)/2)
    poslen=len(pos)
    ### 强行划分 训练集和 测试集 因为这个根本就不叫训练
    trainpos = pos[0:posHalfLen]
    trainneg = neg[0:posHalfLen]
    testpos=pos[posHalfLen:poslen]
    testneg=neg[posHalfLen:poslen]
    
    start = time.process_time()
    ### 改了
    posweight = sq.getMulWeight_mid(testpos,kstart,kend)
    negweight = sq.getMulWeight_mid(testneg,kstart,kend)
    end = time.process_time()
    print(len(posweight))
    print(len(negweight))
    print("权重计算时间：",(end-start))

#    print("--------------------------MulD2-----------------------")
#    print("kstart=",kstart," kend=",kend)
#    start = time.process_time()
#    ## 计算MulD2特征
#    posD2Lis=sq.getMulCount(pos,kstart,kend,pos) 
#    negD2Lis=sq.getMulCount(neg,kstart,kend,neg)
#    
#    ## 构成基因对,并设置标识
#    posPair=[]
#    negPair=[]
#    for i in range(len(posD2Lis)):
#        for j in range(i+1,len(posD2Lis)):
#            posSim=Sim.getMulD2WeightSim2(posD2Lis[i],posD2Lis[j],posweight)
#            negSim=Sim.getMulD2WeightSim2(negD2Lis[i],negD2Lis[j],negweight)
#            posPair.append(["+",posSim])
#            negPair.append(["-",negSim])
#    Pair=posPair+negPair
#    pSim=sorted(Pair,key=getSim,reverse=True)
#    corrCnt=[]
#    corrFre=[]
#    count=0
#    for i in range(top):
#        if pSim[i][0]=="+":
#            count=count+1
#        if (i+1)%10==0 or i+1==top:
#            corrCnt.append(count)
#            corrFre.append(count/(i+1))
#    print("预测准确个数：")
#    print(corrCnt)
#    print("预测准确率：")
#    print(corrFre)
#    end = time.process_time()
#    print("程序运行时间：",(end-start))


    for r in range(0,3):
        start=time.process_time()
        ma=markov.Markov()
        kmer_propos=ma.get_Mulk_Mul_kmer_Pro(trainpos,kstart,kend,r)
        kmer_proneg=ma.get_Mulk_Mul_kmer_Pro(trainneg,kstart,kend,r)

        testkmerpos=sq.getMulSeqKerSet(testpos,kstart,kend)
        testkmerneg=sq.getMulSeqKerSet(testneg,kstart,kend)
        
        ## 训练集的kmer
        for key in testkmerpos:
            if key in kmer_propos.keys():
                testkmerpos[key]=kmer_propos[key]
        for key in testkmerneg:
            if key in kmer_proneg.keys():
                testkmerneg[key]=kmer_proneg[key]
        
        end=time.process_time()
        print("计算马尔可夫概率时间：",(end-start))
        print("--------------------------MulD2s-----------------------")
        print("kstart=",kstart," kend=",kend," r=",r)
        start = time.process_time()
        ## 计算MulD2特征
        start1=time.process_time()
        ## 使用测试集来 做
        posD2Lis=sq.getD2SMulCount_mid(testpos,kstart,kend,r,flag,testkmerpos) 
        negD2Lis=sq.getD2SMulCount_mid(testneg,kstart,kend,r,flag,testkmerneg)
        
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
        ## 计算MulD2特征 测试集合
        posD2Lis=sq.getD2StarMulCount_mid(testpos,kstart,kend,r,flag,testkmerpos) 
        negD2Lis=sq.getD2StarMulCount_mid(testneg,kstart,kend,r,flag,testkmerneg)
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
     
            
        