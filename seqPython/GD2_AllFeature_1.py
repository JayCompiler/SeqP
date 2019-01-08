# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 14:16:53 2019

第一个实验对每个特征加权重

@author: YZi
"""


import random

from deap import base
from deap import creator
from deap import tools
#import numpy as np
import ReadData
import Sequence
from sklearn.metrics import roc_auc_score
import Distance
import time


#name="human_muscle"
#name="fly_blastoderm"
name="human_HBB"
#name="fly_tracheal_system"

## 获取数据集 整个数据集，正数据集，负数据集 都是序列，没有标签
#datasets,pos,neg=rd.getData2("fly_blastoderm")
  ## 获取数据 kmer 从2-->6
print("---------",name,"------------")
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
d2freLis,d2arf=sq.getSeqfreq(datasets,2,d2dic)
d3freLis,d3arf=sq.getSeqfreq(datasets,3,d3dic)
d4freLis,d4arf=sq.getSeqfreq(datasets,4,d4dic)
d5freLis,d5arf=sq.getSeqfreq(datasets,5,d5dic)
d6freLis,d6arf=sq.getSeqfreq(datasets,6,d6dic)


### 将2-6 合成一个字典,并计算权重

size=len(d2freLis)
## 标准化后合并
freqLis=[None]*size
for i in range(len(d2freLis)):
    freqLis[i]=dict(sq.normdata_max_min(d2freLis[i]),**(sq.normdata_max_min(d3freLis[i])))
    freqLis[i]=dict(freqLis[i],**(sq.normdata_max_min(d4freLis[i])))
    freqLis[i]=dict(freqLis[i],**(sq.normdata_max_min(d5freLis[i])))
    freqLis[i]=dict(freqLis[i],**(sq.normdata_max_min(d6freLis[i])))

# 计算权重
weightMax=sq.getWeight(freqLis)


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
    
## 标准化后合并
countLis=[None]*size
for i in range(len(d2freLis)):
    countLis[i]=dict(d2countLis[i],**(d3countLis[i]))
    countLis[i]=dict(countLis[i],**(d4countLis[i]))
    countLis[i]=dict(countLis[i],**(d5countLis[i]))
    countLis[i]=dict(countLis[i],**(d6countLis[i]))

weightSize=len(countLis[0])    

## 构造训练集合和测试集合
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


## 产生随机种子，并将数据集打乱
#randnum = random.randint(0,100)
#random.seed(randnum)
## 打乱正负序列集 
random.shuffle(pairPoslist)
random.shuffle(pairNeglist)

## 划分训练集合和测试集合  并设置标签
trainLis=[]
trainlabel=[]
testLis=[]
testlabel=[]
trainPosSize=int(0.8*len(pairPoslist))
for i in range(len(pairPoslist)):
    if i<=trainPosSize:
        trainLis.append(pairPoslist[i])
        trainlabel.append(1)
        trainLis.append(pairNeglist[i])
        trainlabel.append(0)
    else:
        testLis.append(pairPoslist[i])
        testlabel.append(1)
        testLis.append(pairNeglist[i])
        testlabel.append(0)

##-----------------------------------------------------------------------------    


creator.create("FitnessMax", base.Fitness, weights=(1.0,))   
creator.create("Individual", list, fitness=creator.FitnessMax)    

toolbox = base.Toolbox()    


toolbox.register("attr_bool", random.random)   

## 设置权重个数
toolbox.register("individual", tools.initRepeat, creator.Individual,    #tools.initRepeat是干
    toolbox.attr_bool, weightSize)


toolbox.register("population", tools.initRepeat, list, toolbox.individual)
 


# register the crossover operator
toolbox.register("mate", tools.cxTwoPoint)

def evalOneMax(individual):
    weight = individual 
#    su=sum(weight)
#    for i in range(len(weight)):
#        weight[i]=(weight[i]+np.spacing(1))/su+np.spacing(1)
    sim=[]
    ## 计算权重
  
    for i in range(len(trainLis)):
        tmp_sim=0
        count=0
        for key in sorted(countLis[0]):
             tmp_sim=tmp_sim+countLis[trainLis[i][0]][key]*countLis[trainLis[i][1]][key]*weight[count]
             count=count+1
        sim.append(tmp_sim) 
    auc=roc_auc_score(trainlabel, sim)
    return auc,


#    return float(np.sum(np.abs(np.array(sim)-label))),
    
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
toolbox.register("select", tools.selTournament, tournsize=4)    #这里选择的tournsize又是什么意思呢？

#----------

def main():
    random.seed(64)
    # hash(64) is used
    # random.seed方法的作用是给随机数对象一个种子值，用于产生随机序列。
    # 对于同一个种子值的输入，之后产生的随机数序列也一样。
    # 通常是把时间秒数等变化值作为种子值，达到每次运行产生的随机系列都不一样
    # create an initial population of 300 individuals (where
    # each individual is a list of integers)
    pop = toolbox.population(n=50)    #定义了300个个体的种群！！！
    
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
    CXPB, MUTPB, NGEN = 0.6, 0.3, 25
    
    print("Start of evolution")
    
    # Evaluate the entire population
    fitnesses = list(map(toolbox.evaluate, pop))
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
        fitnesses = map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit
        
        print("Evaluated %i individuals" % len(invalid_ind))
        
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
    
    print("-- End of (successful) evolution ---")
    
    best_ind = tools.selBest(pop, 1)[0]
#    print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))
    
    tp=sum(best_ind)
    for i in range(len(best_ind)):
        best_ind[i]=best_ind[i]/tp
    return best_ind
#    file= open('test.txt', 'w')
#    for w in best_ind:
#        file.write(str(w))
#        file.write('\n')
#    file.close()

if __name__ == "__main__":
    
    print(len(freqLis[0]))
    dis =Distance.Distance()
    print("GA---------------------------")
    ## 归一化处理
    start= time.process_time()
    w=main()
    end= time.process_time()
    print("求权重时间消耗：",(end-start))
    count=0
#    c1=len(d2countLis[0])
#    c2=len(d3countLis[0])
#    c3=len(d4countLis[0])
#    c4=len(d5countLis[0])
#    c5=len(d6countLis[0])
    for k in range(2,7):
        nam=str(k)+"-weight-"+name
        file= open(nam, 'w')
        feature="d"+str(k)+"freLis"
        for i in range(len(eval(feature)[0])):
            file.write(str(w[count]))
            count=count+1
            file.write('\n')
        file.close()
    ## 保存完毕
    print("----------------------保存完毕-------------")
    print(len(w))
    sim=[]
    ## 计算权重  结合遗传算法
    sim=[]
    inner2=[]
    inner3=[]
    inner4=[]
    inner5=[]
    inner6=[]
    
    sim=[]
    ## 计算权重
  
    for i in range(len(testLis)):
        tmp_sim=0
        count=0
        for key in sorted(countLis[0]):
            tmp_sim=tmp_sim+countLis[testLis[i][0]][key]*countLis[testLis[i][1]][key]*w[count]
            count=count+1
        sim.append(tmp_sim)               
    auc=roc_auc_score(testlabel, sim)
    
    
    print(auc)
    print(len(d2countLis[0])+len(d3countLis[0])+len(d4countLis[0])+len(d5countLis[0])+len(d6countLis[0]))

    for k in range(2,7):
        ## eu:
       
        feature="d"+str(k)+"freLis"
        print("----------------EU-----------------")
        print("k=",k)
        sim=[]
        for i in range(len(testLis)):
            tmp=1/(dis.EuD_feature(eval(feature)[testLis[i][0]],eval(feature)[testLis[i][1]]))
            sim.append(tmp)
        auc=roc_auc_score(testlabel, sim)  
        print(auc)
        
        print("----------------Ma-----------------")
        print("k=",k)
        sim=[]
        for i in range(len(testLis)):
            tmp=1/(dis.manhattan_feature(eval(feature)[testLis[i][0]],eval(feature)[testLis[i][1]]))
            sim.append(tmp)
        auc=roc_auc_score(testlabel, sim)  
        print(auc)
        
        
        
        print("----------------kld-----------------")
        print("k=",k)
        sim=[]
        for i in range(len(testLis)):
            tmp=1/(dis.KLD_feature(eval(feature)[testLis[i][0]],eval(feature)[testLis[i][1]]))
            sim.append(tmp)
        auc=roc_auc_score(testlabel, sim)  
        print(auc)
        
        
        print("----------------pcc-----------------")
        print("k=",k)
        sim=[]
        for i in range(len(testLis)):
            tmp=1/(dis.pcc_feature(eval(feature)[testLis[i][0]],eval(feature)[testLis[i][1]]))
            sim.append(tmp)
        auc=roc_auc_score(testlabel, sim)  
        print(auc)
        
    
        print("----------------cosine-----------------")
        print("k=",k)
        sim=[]
        for i in range(len(testLis)):
            tmp=1/(dis.cosine_feature(eval(feature)[testLis[i][0]],eval(feature)[testLis[i][1]]))
            sim.append(tmp)
        auc=roc_auc_score(testlabel, sim)  
        print(auc)
        
        
        print("----------------D2-----------------")
        print("k=",k)
        feature="d"+str(k)+"countLis"
        
        sim=[]
        for i in range(len(testLis)):
            tmp=1/(dis.getD2_feature(eval(feature)[testLis[i][0]],eval(feature)[testLis[i][1]]))
            sim.append(tmp)
        auc=roc_auc_score(testlabel, sim)  
        print(auc)
        
     
