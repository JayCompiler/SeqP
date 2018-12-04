# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 16:17:42 2018

@author: Yzi

"""

# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 14:06:56 2018

@author: Yzi
单条 markov
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
    
    #top 参数 36
#    name="fly_tracheal_system"
    
    
     #top 参数 300
#    name="human_muscle"
     
     
     
        #top 参数 153
#    name="fly_pns" or "human_HBB"
        
        
        
                #top 参数 136
    name="fly_eye" 
#or "human_HBB" 
#    name="human_HBB"

                
                
                #top =36
#    name="human_liver" 
# or"fly_tracheal_system"
#    name="fly_tracheal_system"
#    name="fly_pns"
    #top =300
#    name="fly_blastoderm"

    
    print("-----------------------",name,"-----------------------------")
    datasets,pos,neg=rd.getData2(name)
#    print(len(datasets))
    Sim=Similarity.Similarity()
    kstart=2
    kend=8
    top=136
    flag=True
    sq=Sequence.Sequence()


    print("--------------------------MulD2-----------------------")
    print("kstart=",kstart," kend=",kend)
    start = time.process_time()
    ## 计算MulD2特征
#    posD2Lis=sq.getMulCount(pos,kstart,kend,pos) 
#    negD2Lis=sq.getMulCount(neg,kstart,kend,neg)
    ## 构成基因对,并设置标识
    posPair=[]
    negPair=[]
    for i in range(len(pos)):
        for j in range(i+1,len(pos)):
            posSeqLis=[]
            negSeqLis=[]
            posSeqLis.append(pos[i])
            posSeqLis.append(pos[j])
            negSeqLis.append(neg[i])
            negSeqLis.append(neg[j])
            # 获得特征
            posD2Lis=sq.getMulCount_nonorm(posSeqLis,kstart,kend,posSeqLis)
            negD2Lis=sq.getMulCount_nonorm(negSeqLis,kstart,kend,negSeqLis)
            ## 计算权重
            posweight = sq.getMulWeight_nonorm(posSeqLis,kstart,kend) 
            negweight = sq.getMulWeight_nonorm(negSeqLis,kstart,kend)
            # getMulD2WeightSim2 接受一切处理好的特征
            posSim=Sim.getMulD2WeightSim2(posD2Lis[0],posD2Lis[1],posweight)
            negSim=Sim.getMulD2WeightSim2(negD2Lis[0],negD2Lis[1],negweight)
            posPair.append(["+",posSim])
            negPair.append(["-",negSim])
    Pair=posPair+negPair
    pSim=sorted(Pair,key=getSim,reverse=True)
    corrCnt=[]
    corrFre=[]
    count=0
#    print("asd",len(pSim))
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
        ## 计算MulD2s特征
        posPair=[]
        negPair=[]
        for i in range(len(pos)):
            for j in range(i+1,len(pos)):        
                posSeqLis=[]
                negSeqLis=[]
                posSeqLis.append(pos[i])
                posSeqLis.append(pos[j])
                negSeqLis.append(neg[i])
                negSeqLis.append(neg[j])
                # 计算概率
                start=time.process_time()
                ma=markov.Markov()
                kmer_propos=ma.get_Mulk_Mul_kmer_Pro(posSeqLis,kstart,kend,r)
                kmer_proneg=ma.get_Mulk_Mul_kmer_Pro(negSeqLis,kstart,kend,r)
                end=time.process_time()
                # 获得特征
                posD2Lis=sq.getD2SMulCount_nonorm(posSeqLis,kstart,kend,r,True,kmer_propos)
                negD2Lis=sq.getD2SMulCount_nonorm(negSeqLis,kstart,kend,r,True,kmer_proneg)
                ## 计算权重
                posweight = sq.getMulWeight_nonorm(posSeqLis,kstart,kend) 
                negweight = sq.getMulWeight_nonorm(negSeqLis,kstart,kend)
                
                posSim=Sim.getMulD2sWeightSim2(posD2Lis[0],posD2Lis[1],posweight)
                negSim=Sim.getMulD2sWeightSim2(negD2Lis[0],negD2Lis[1],negweight)   #修改过
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
        ## 计算MulD2star特征
        ## 构成基因对,并设置标识
        posPair=[]
        negPair=[]
        for i in range(len(pos)):
            for j in range(i+1,len(pos)):
                
                posSeqLis=[]
                negSeqLis=[]
                posSeqLis.append(pos[i])
                posSeqLis.append(pos[j])
                negSeqLis.append(neg[i])
                negSeqLis.append(neg[j])
                # 计算概率
                start=time.process_time()
                ma=markov.Markov()
                kmer_propos=ma.get_Mulk_Mul_kmer_Pro(posSeqLis,kstart,kend,r)
                kmer_proneg=ma.get_Mulk_Mul_kmer_Pro(negSeqLis,kstart,kend,r)
                end=time.process_time()
                # 获得特征
                posD2Lis=sq.getD2StarMulCount_nonorm(posSeqLis,kstart,kend,r,True,kmer_propos)
                negD2Lis=sq.getD2StarMulCount_nonorm(negSeqLis,kstart,kend,r,True,kmer_proneg)
                ## 计算权重
                posweight = sq.getMulWeight_nonorm(posSeqLis,kstart,kend) 
                negweight = sq.getMulWeight_nonorm(negSeqLis,kstart,kend)
                
                posSim=Sim.getMulD2starWeightSim2(posD2Lis[0],posD2Lis[1],posweight)
                negSim=Sim.getMulD2starWeightSim2(negD2Lis[0],negD2Lis[1],negweight)
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

