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

## 返回列表的第三个 元素，即相似度
def getSim(Lis):
    return Lis[1]

if __name__=="__main__":
    rd=ReadData.ReadData()
    datasets,pos,neg=rd.getData2("fly_blastoderm")
    Sim=Similarity.Similarity()
    k=6
    flag=True
    sq=Sequence.Sequence()
    ma=markov.Markov()

    poskmerset,poskmersetdic=sq.getSeqKerSet(pos,k)
    
    negkmerset,negkmersetdic=sq.getSeqKerSet(neg,k)
    for r in range(0,3):
        posKmer_pro=ma.get_Mul_kmer_Pro(pos,k,r)
        negKmer_pro=ma.get_Mul_kmer_Pro(neg,k,r)
        print("--------------------------MulD2s-----------------------")
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
                posSim=Sim.getD2sSim(pos[i],pos[j],k,r,flag,poskmersetdic,posKmer_pro)
                negSim=Sim.getD2sSim(neg[i],neg[j],k,r,flag,negkmersetdic,negKmer_pro)
                posPair.append(["+",posSim])
                negPair.append(["-",negSim])
        Pair=posPair+negPair
        pSim=sorted(Pair,key=getSim,reverse=True)
        corrCnt=[]
        corrFre=[]
        count=0
        for i in range(300):
            if pSim[i][0]=="+":
                count=count+1
            if (i+1)%10==0:
                corrCnt.append(count)
                corrFre.append(count/(i+1))
        print("预测准确个数：")
        print(corrCnt)
        print("预测准确率：")
        print(corrFre)
        end = time.process_time()
        print("程序运行时间：",(end-start))
        
        
        
      
        print("--------------------------MulD2star-----------------------")
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
        for i in range(300):
            if pSim[i][0]=="+":
                count=count+1
            if (i+1)%10==0:
                corrCnt.append(count)
                corrFre.append(count/(i+1))
        print("预测准确个数：")
        print(corrCnt)
        print("预测准确率：")
        print(corrFre)
        end = time.process_time()
        print("程序运行时间：",(end-start))
        
        
       
     
            
        