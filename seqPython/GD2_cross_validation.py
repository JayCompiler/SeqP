# -*- coding: utf-8 -*-
"""
Created on Fri Jan 11 14:25:38 2019

@author: Yzi
"""

import random
from deap import base
from deap import creator
from deap import tools
import numpy as np
import ReadData
import Sequence
from sklearn.metrics import roc_auc_score
import Distance
import time


## 将lis划分为num块
def chunkIt(Lis, num):
  random.shuffle(Lis)
  avg = len(Lis) / float(num)
  out = []
  last = 0.0

  while last < len(Lis):
    out.append(Lis[int(last):int(last + avg)])
    last += avg

  return out


def evalOneMax(individual):
    weight = individual 
    su=sum(weight)
    for i in range(len(weight)):
        weight[i]=(weight[i]+np.spacing(1))/su+np.spacing(1)
    sim=[]
    ## 计算权重
  
    for i in range(len(trainLis)):
        tmp_sim=0
        for key in dict.keys(d2countLis[0]):
            tmp_sim=tmp_sim+d2countLis[trainLis[i][0]][key]*d2countLis[trainLis[i][1]][key]*weight2M[key]*weight[0]
        for key in dict.keys(d3countLis[0]):
            tmp_sim=tmp_sim+d3countLis[trainLis[i][0]][key]*d3countLis[trainLis[i][1]][key]*weight3M[key]*weight[1]
        for key in dict.keys(d4countLis[0]):
            tmp_sim=tmp_sim+d4countLis[trainLis[i][0]][key]*d4countLis[trainLis[i][1]][key]*weight4M[key]*weight[2]    
        for key in dict.keys(d5countLis[0]):
            tmp_sim=tmp_sim+d5countLis[trainLis[i][0]][key]*d5countLis[trainLis[i][1]][key]*weight5M[key]*weight[3]
        for key in dict.keys(d6countLis[0]):
            tmp_sim=tmp_sim+d6countLis[trainLis[i][0]][key]*d6countLis[trainLis[i][1]][key]*weight6M[key]*weight[4] 
        sim.append(tmp_sim)               
    auc=roc_auc_score(trainlabel, sim)
    return auc,


##------------------------注册遗传算法相关的包函数----------------------------------------------------- -----   


creator.create("FitnessMax", base.Fitness, weights=(1.0,))   
creator.create("Individual", list, fitness=creator.FitnessMax)    

toolbox = base.Toolbox()    


toolbox.register("attr_bool", random.random)   


toolbox.register("individual", tools.initRepeat, creator.Individual,    #tools.initRepeat是干嘛的？？？
    toolbox.attr_bool, 5)


toolbox.register("population", tools.initRepeat, list, toolbox.individual)
 


# register the crossover operator
toolbox.register("mate", tools.cxTwoPoint)
    
#----------
# Operator registration
#----------
# register the goal / fitness function
# 这里的toolbox register语句的理解：注册了一个函数evaluae依据的是后面的evalOneMax 理通了!!!
toolbox.register("evaluate", evalOneMax)

# register a mutation operator with a probability to
# flip each attribute/gene of 0.05
toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)

# operator for selecting individuals for breeding the next
# generation: each individual of the current generation
# is replaced by the 'fittest' (best) of three individuals
# drawn randomly from the current generation.
toolbox.register("select", tools.selTournament, tournsize=4)    #这里选择的tournsize又是什么意思呢


##------------------------------------------------遗传函数主体
def main():
    random.seed(64)
    # hash(64) is used
    # random.seed方法的作用是给随机数对象一个种子值，用于产生随机序列。
    # 对于同一个种子值的输入，之后产生的随机数序列也一样。
    # 通常是把时间秒数等变化值作为种子值，达到每次运行产生的随机系列都不一样
    # create an initial population of 300 individuals (where
    # each individual is a list of integers)
    pop = toolbox.population(n=100)    #定义了300个个体的种群！！！
    
#    print(pop[0][0])
    for i in range(len(pop)):
        su=sum(pop[i])
        for j in range(len(pop[i])):
            pop[i][j]=pop[i][j]/su

    # CXPB  is the probability with which two individuals
    #       are crossed
    #
    # MUTPB is the probability for mutating an individual
    #
    # NGEN  is the number of generations for which the
    #       evolution runs   进化运行的代数！
    CXPB, MUTPB, NGEN = 0.6, 0.3, 20
    
    print("Start of evolution")
    
    # Evaluate the entire population
    fitnesses = list(map(toolbox.evaluate, pop)) ###-------------------------------
    for ind, fit in zip(pop, fitnesses):
        ind.fitness.values = fit
    
    print("  Evaluated %i individuals" % len(pop))   
    
    # Begin the evolution      
    for g in range(NGEN):
        print("-- Generation %i --" % g)
        
        # Select the next generation individuals
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(map(toolbox.clone, offspring))
        
        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):

            # cross two individuals with probability CXPB
            if random.random() < CXPB:
                toolbox.mate(child1, child2)

                del child1.fitness.values
                del child2.fitness.values

        for mutant in offspring:

            # mutate an individual with probability MUTPB
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values
    
        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fitnesses = map(toolbox.evaluate, invalid_ind)  #------------------------------------2-----------
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        
        print("  Evaluated %i individuals" % len(invalid_ind))
        
        # The population is entirely replaced by the offspring
        pop[:] = offspring
        
        # Gather all the fitnesses in one list and print the stats
        fits = [ind.fitness.values[0] for ind in pop]
        
        length = len(pop)
        mean = sum(fits) / length
        sum2 = sum(x*x for x in fits)
        std = abs(sum2 / length - mean**2)**0.5
        
        print("  Min %s" % min(fits))
        print("  Max %s" % max(fits))
        print("  Avg %s" % mean)
        print("  Std %s" % std)
        if std<=10E-8:
            break
        
    print("-- End of (successful) evolution ---")
    
    best_ind = tools.selBest(pop, 1)[0]
#    print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))
    
    tp=sum(best_ind)
    for i in range(len(best_ind)):
        best_ind[i]=best_ind[i]/tp
    return best_ind






if __name__ == "__main__":   
    ###--------------------------数据的预处理--------------------------1-----------------
#    name="human_muscle"
    
#    name="fly_blastoderm"
    
#    name="human_HBB"
    name="pns"
    
    ## 获取数据集 整个数据集，正数据集，负数据集 都是序列，没有标签
    #datasets,pos,neg=rd.getData2("fly_blastoderm")
      ## 获取数据 kmer 从2-->6
    print("------GA_D2---",name,"------------")
    rd=ReadData.ReadData()
    
    datasets,pos,neg=rd.getData2(name)
    ## 正例数据集条数
    possize=len(pos)
    ## 合并两个list
    datasets=pos+neg
    
    sq=Sequence.Sequence()
    
    flag=True
        ## 获取整个数据集的kmer集合,dict
        
        ## 获取字典集合
    d2set,d2dic=sq.getSeqKerSet(datasets,2)
    d3set,d3dic=sq.getSeqKerSet(datasets,3)
    d4set,d4dic=sq.getSeqKerSet(datasets,4)
    d5set,d5dic=sq.getSeqKerSet(datasets,5)
    d6set,d6dic=sq.getSeqKerSet(datasets,6)
    
    
    ## d2 以整个数据集计算 计算频率 不需要标准化，用于计算单个k的其他距离定义方法
    ## 获得权重--------------------------------------------------------------------------
    d2freLis,d2arf=sq.getSeqfreq(datasets,2,d2dic)
    weight2M=sq.getWeight(d2freLis)
    d3freLis,d3arf=sq.getSeqfreq(datasets,3,d3dic)
    weight3M=sq.getWeight(d3freLis)
    d4freLis,d4arf=sq.getSeqfreq(datasets,4,d4dic)
    weight4M=sq.getWeight(d4freLis)
    d5freLis,d5arf=sq.getSeqfreq(datasets,5,d5dic)
    weight5M=sq.getWeight(d5freLis)
    d6freLis,d6arf=sq.getSeqfreq(datasets,6,d6dic)
    weight6M=sq.getWeight(d6freLis)
    
    
    size=len(d2freLis)
    
    
    
        ## d2 特征 count 所有数据集
    d2countLis,d2arc=sq.getSeqCount(datasets,2,d2dic)
    d3countLis,d3arc=sq.getSeqCount(datasets,3,d3dic)
    d4countLis,d4arc=sq.getSeqCount(datasets,4,d4dic)
    d5countLis,d5arc=sq.getSeqCount(datasets,5,d5dic)
    d6countLis,d6arc=sq.getSeqCount(datasets,6,d6dic)
    ## 标准化过程： 标准化 count
    for i in range(len(datasets)):
        d2countLis[i]=sq.normdata_max_min(d2countLis[i])
        d3countLis[i]=sq.normdata_max_min(d3countLis[i])
        d4countLis[i]=sq.normdata_max_min(d4countLis[i])
        d5countLis[i]=sq.normdata_max_min(d5countLis[i])
        d6countLis[i]=sq.normdata_max_min(d6countLis[i])
    
    # 构造对集合
    pairPoslist=[]
    pairNeglist=[]
    for i in range(possize):
        # 装填正例
        for j in range(i+1,possize):
            tpos=[i,j,1]
            pairPoslist.append(tpos)
    # 装填负例
    for i in range(possize):
        # 装填正例
        for j in range(i+1,possize):
            tneg=[i+possize,j+possize,0]
            pairNeglist.append(tneg)
    
    
    
    ###构造数据集，交叉验证,n折划分---------------------------------------2-----------------------------
    num=5
    posChunk=chunkIt(pairPoslist,num)
    negChunk=chunkIt(pairNeglist,num)
    start=time.process_time()
    aucLis=[]
    for  o in range(num):
        ## 构造训练集合
        testPosLis=posChunk[o]
        testNegLis=negChunk[o]
        trainPosLis=[]
        trainNegLis=[]
        j=o+1
        q=o-1
        while j<num:
            trainPosLis=trainPosLis+posChunk[j]
            trainNegLis=trainNegLis+negChunk[j]
            j=j+1
        while q>=0:
            trainPosLis=trainPosLis+posChunk[q]
            trainNegLis=trainNegLis+negChunk[q]
            q=q-1
            j=j+1
        testLis=testPosLis+testNegLis
        trainLis=trainPosLis+trainNegLis
        Tes=np.array(testLis)
        Tra=np.array(trainLis)
        testlabel=Tes[:,2]
        trainlabel=Tra[:,2]
        
        
        print("训练集大小：",len(trainlabel))
        print("测试集大小：",len(testlabel))
        ##---------------------------------------------------------3 遗传开始-------------------------
        dis =Distance.Distance()
       
        print("GA--------------第",o,"次验证---------------" )
        ## 归一化处理
        w=main()
    
        print("权重个数----",len(w))
        sim=[]
        ## 计算权重
      
        for i in range(len(testLis)):
            tmp_sim=0
            for key in dict.keys(d2countLis[0]):
                tmp_sim=tmp_sim+d2countLis[testLis[i][0]][key]*d2countLis[testLis[i][1]][key]*weight2M[key]*w[0]
            for key in dict.keys(d3countLis[0]):
                tmp_sim=tmp_sim+d3countLis[testLis[i][0]][key]*d3countLis[testLis[i][1]][key]*weight3M[key]*w[1]    
            for key in dict.keys(d4countLis[0]):
                tmp_sim=tmp_sim+d4countLis[testLis[i][0]][key]*d4countLis[testLis[i][1]][key]*weight4M[key]*w[2]                 
            for key in dict.keys(d5countLis[0]):
                tmp_sim=tmp_sim+d5countLis[testLis[i][0]][key]*d5countLis[testLis[i][1]][key]*weight5M[key]*w[3]                
            for key in dict.keys(d6countLis[0]):
                tmp_sim=tmp_sim+d6countLis[testLis[i][0]][key]*d6countLis[testLis[i][1]][key]*weight6M[key]*w[4]
            sim.append(tmp_sim)     
        auc=roc_auc_score(testlabel, sim)    
        aucLis.append(auc)
        print("GA_D2第",o,"次auc的值",auc)
    su=0
    for i in range(len(aucLis)):
        su=su+aucLis[i]
    avgAuc=su/num
    print("GA_D2平均auc",avgAuc)
    end=time.process_time()
    print("程序时间：",(end-start))
    print("特征个数",len(d2countLis[0])+len(d3countLis[0])+len(d4countLis[0])+len(d5countLis[0])+len(d6countLis[0]))
    
    
    
