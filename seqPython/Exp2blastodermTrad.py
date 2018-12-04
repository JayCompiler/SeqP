# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 17:28:26 2018

@author: Yzi
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 10:50:07 2018

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
#    name="fly_blastoderm"
    
    
        #top 参数 36
#    name="fly_tracheal_system"
    
    
     #top 参数 300
#    name="human_muscle"
     
     
     
        #top 参数 153
#    name="fly_pns" 
#or
#    name= "human_HBB"

        
        
        
                #top 参数 136
    name="fly_eye" 
#or "human_HBB" 
                
                
                #top =36
                #    name="human_liver" or"fly_tracheal_system"
#                  name="human_liver"   
    #top=300
#    name="fly_tracheal_system"
    
    datasets,pos,neg=rd.getData2(name)
    
    Sim=Similarity.Similarity()
    k=4
    flag=True
    top=136
    sq=Sequence.Sequence()
    ma=markov.Markov()

    poskmerset,poskmersetdic=sq.getSeqKerSet(pos,k)
    
    negkmerset,negkmersetdic=sq.getSeqKerSet(neg,k)
    
    posPair=[]
    negPair=[]
    print("--------------------------D2-----------------------")
    print("k=",k)
    start = time.process_time()
    for i in range(len(pos)):
        for j in range(i+1,len(pos)):
            posSim=Sim.getD2Sim(pos[i],pos[j],k)
            negSim=Sim.getD2Sim(neg[i],neg[j],k)
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
        posKmer_pro=ma.get_Mul_kmer_Pro(pos,k,r)
        negKmer_pro=ma.get_Mul_kmer_Pro(neg,k,r)
        print("--------------------------D2s-----------------------")
        print("k=",k,"r=",r)
        start = time.process_time()
        
        ## 计算MulD2特征
        ## 构成基因对,并设置标识
        posPair=[]
        negPair=[]
        for i in range(len(pos)):
            for j in range(i+1,len(pos)):
                posSim=Sim.getD2sSim(pos[i],pos[j],k,r,flag,poskmersetdic,posKmer_pro)
                negSim=Sim.getD2sSim(neg[i],neg[j],k,r,flag,negkmersetdic,negKmer_pro)
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
        
        
        
      
        print("--------------------------D2star-----------------------")
        print("k=",k,"r=",r)
        start = time.process_time()
        
        ## 计算MulD2特征
#        posD2Lis=sq.getD2SMulCount2(pos,kstart,kend,r,flag) 
#        negD2Lis=sq.getD2SMulCount2(neg,kstart,kend,r,flag)
        ## 构成基因对,并设置标识
        posPair=[]
        negPair=[]
        for i in range(len(pos)):
            for j in range(i+1,len(pos)):
                posSim=Sim.getD2starSim(pos[i],pos[j],k,r,flag,pos,poskmersetdic,posKmer_pro)
                negSim=Sim.getD2starSim(neg[i],neg[j],k,r,flag,neg,negkmersetdic,negKmer_pro)
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
        
        
       
     
            
        