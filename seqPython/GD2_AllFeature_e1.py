# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 17:08:54 2018

@author: 51164
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 14:47:49 2018

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

  ## 获取数据 kmer 从2-->6
rd=ReadData.ReadData()
query,dataset,label=rd.getData("dataset1")
    ## 产生随机种子，并将数据集打乱
randnum = random.randint(0,100)
random.seed(randnum)
    
random.shuffle(dataset)
random.seed(randnum)
random.shuffle(label)
    
     ### 划分 训练集和测试集
trainset=[]
trainlabel=[]
testset=[]
testlabel=[]
poscount=0
negcount=0

    # 按照3：1的比例划分训练集合和训练集合
for i in range(len(dataset)):
    if (label[i]==1 and poscount<16) or (label[i]==0 and negcount<16):
        trainset.append(dataset[i])
        trainlabel.append(label[i])
        if label[i]==1:
            poscount=poscount+1
        else:
            negcount=negcount+1
    else:
        testset.append(dataset[i])
        testlabel.append(label[i])
    ## 训练集合加上 query sequence
trainSets=list.copy(trainset)
trainSets.append(query[0])
     ## 测试集合加上 query sequence
testSets=list.copy(testset)
testSets.append(query[0])
    ## 整体数据集加上 query sequence
sq=Sequence.Sequence()
datasets=list.copy(dataset)
datasets.append(query[0])
flag=True
    ## 获取kmer集合,dict
d2set,d2dic=sq.getSeqKerSet(datasets,2)
d3set,d3dic=sq.getSeqKerSet(datasets,3)
d4set,d4dic=sq.getSeqKerSet(datasets,4)
d5set,d5dic=sq.getSeqKerSet(datasets,5)
d6set,d6dic=sq.getSeqKerSet(datasets,6)
    ## d2 频率 训练集
d2freLis,d2arf=sq.getSeqfreq(trainSets,2,d2dic)
d3freLis,d3arf=sq.getSeqfreq(trainSets,3,d3dic)
d4freLis,d4arf=sq.getSeqfreq(trainSets,4,d4dic)
d5freLis,d5arf=sq.getSeqfreq(trainSets,5,d5dic)
d6freLis,d6arf=sq.getSeqfreq(trainSets,6,d6dic)



    ## d2 频率 测试集
d2freLisTest,d2arfTest=sq.getSeqfreq(testSets,2,d2dic)
d3freLisTest,d3arfTest=sq.getSeqfreq(testSets,3,d3dic)
d4freLisTest,d4arfTest=sq.getSeqfreq(testSets,4,d4dic)
d5freLisTest,d5arfTest=sq.getSeqfreq(testSets,5,d5dic)
d6freLisTest,d6arfTest=sq.getSeqfreq(testSets,6,d6dic)
    ## d2 特征 与权重
d2countLis,d2arc=sq.getSeqCount(trainSets,2,d2dic)
w2dic=sq.getWeight(d2freLis)
d3countLis,d3arc=sq.getSeqCount(trainSets,3,d3dic)
w3dic=sq.getWeight(d3freLis)
d4countLis,d4arc=sq.getSeqCount(trainSets,4,d4dic)
w4dic=sq.getWeight(d4freLis)
d5countLis,d5arc=sq.getSeqCount(trainSets,5,d5dic)
w5dic=sq.getWeight(d5freLis)
d6countLis,d6arc=sq.getSeqCount(trainSets,6,d6dic)
w6dic=sq.getWeight(d6freLis)

    ## d2 特征 测试集
d2countLisTest,d2arcTest=sq.getSeqCount(testSets,2,d2dic)
w2dicTest=sq.getWeight(d2freLisTest)
d3countLisTest,d3arcTest=sq.getSeqCount(testSets,3,d3dic)
w3dicTest=sq.getWeight(d3freLisTest)
d4countLisTest,d4arcTest=sq.getSeqCount(testSets,4,d4dic)
w4dicTest=sq.getWeight(d4freLisTest)
d5countLisTest,d5arcTest=sq.getSeqCount(testSets,5,d5dic)
w5dicTest=sq.getWeight(d5freLisTest)
d6countLisTest,d6arcTest=sq.getSeqCount(testSets,6,d6dic)
w6dicTest=sq.getWeight(d6freLisTest)




creator.create("FitnessMax", base.Fitness, weights=(1.0,))   
creator.create("Individual", list, fitness=creator.FitnessMax)    

toolbox = base.Toolbox()    


toolbox.register("attr_bool", random.random)   


toolbox.register("individual", tools.initRepeat, creator.Individual,    #
    toolbox.attr_bool, 5456)


toolbox.register("population", tools.initRepeat, list, toolbox.individual)
 


# register the crossover operator
toolbox.register("mate", tools.cxTwoPoint)

def evalOneMax(individual):
    weight = individual 
    su=sum(weight)
    for i in range(len(weight)):
        weight[i]=(weight[i]+np.spacing(1))/su+np.spacing(1)
    sim=[]
    
    ## 计算权重
    for i in range(len(trainset)):
        count=0
        tmp_sim=0
        for key in sorted(d2countLis[0]):
            tmp_sim=tmp_sim+d2countLis[i][key]*d2countLis[len(d2countLis)-1][key]*w2dic[key]*weight[count]
            count=count+1
        for key in sorted(d3countLis[0]):
            tmp_sim=tmp_sim+d3countLis[i][key]*d3countLis[len(d2countLis)-1][key]*w3dic[key]*weight[count]
            count=count+1
        for key in sorted(d4countLis[0]):
            tmp_sim=tmp_sim+d4countLis[i][key]*d4countLis[len(d2countLis)-1][key]*w4dic[key]*weight[count]
            count=count+1
        for key in sorted(d5countLis[0]):
            tmp_sim=tmp_sim+d5countLis[i][key]*d5countLis[len(d2countLis)-1][key]*w5dic[key]*weight[count]
            count=count+1
        for key in sorted(d6countLis[0]):
            tmp_sim=tmp_sim+d6countLis[i][key]*d6countLis[len(d2countLis)-1][key]*w6dic[key]*weight[count]
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
    CXPB, MUTPB, NGEN = 0.6, 0.5, 70
    
    
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
    dis =Distance.Distance()
    ## 归一化处理
    w=main()
    su=sum(w)
    for i in range(len(w)):
        w[i]=w[i]/su
    print(w)
    print(len(w))
    sim=[]
    ## 统计内积 大小
    inner2=[]
    inner3=[]
    inner4=[]
    inner5=[]
    inner6=[]
    
    
    for i in range(len(testset)):
        tmp_sim=0
        count=0
        ti2=0
        ti3=0
        ti4=0
        ti5=0
        ti6=0
        
        for key in sorted(d2countLisTest[0]):
            tmp_sim=tmp_sim+d2countLisTest[i][key]*d2countLisTest[len(d2countLisTest)-1][key]*w2dicTest[key]*w[count]
            ti2=ti2+d2countLisTest[i][key]*d2countLisTest[len(d2countLisTest)-1][key]
            count=count+1
        for key in sorted(d3countLisTest[0]):
            tmp_sim=tmp_sim+d3countLisTest[i][key]*d3countLisTest[len(d3countLisTest)-1][key]*w3dicTest[key]*w[count]
            ti3=ti3+d3countLisTest[i][key]*d3countLisTest[len(d3countLisTest)-1][key]
            count=count+1
        for key in sorted(d4countLisTest[0]):
            tmp_sim=tmp_sim+d4countLisTest[i][key]*d4countLisTest[len(d4countLisTest)-1][key]*w4dicTest[key]*w[count]
            ti4=ti4+d4countLisTest[i][key]*d4countLisTest[len(d4countLisTest)-1][key]
            count=count+1
        for key in sorted(d5countLisTest[0]):
            tmp_sim=tmp_sim+d5countLisTest[i][key]*d5countLisTest[len(d5countLisTest)-1][key]*w5dicTest[key]*w[count]
            ti5=ti5+d5countLisTest[i][key]*d5countLisTest[len(d5countLisTest)-1][key]
            count=count+1
        for key in sorted(d6countLisTest[0]):
            tmp_sim=tmp_sim+d6countLisTest[i][key]*d6countLisTest[len(d6countLisTest)-1][key]*w6dicTest[key]*w[count]
            ti6=ti6+d6countLisTest[i][key]*d6countLisTest[len(d6countLisTest)-1][key]
            count=count+1
        sim.append(tmp_sim)
    auc=roc_auc_score(testlabel, sim)
    print(auc)
    print(len(d2countLisTest[0])+len(d3countLisTest[0])+len(d4countLisTest[0])+len(d5countLisTest[0])+len(d6countLisTest[0]))

    for k in range(2,7):
        ## eu:
        print("----------------EU-----------------")
        print("k=",k)
        feature="d"+str(k)+"freLisTest"
        sim=[]
        for i in range(len(testset)):
            si=1/(dis.EuD_feature(eval(feature)[i],eval(feature)[len(eval(feature))-1]))
            sim.append(si)  
#        print(sim)
        auc=roc_auc_score(testlabel, sim)
        print(auc)       
        
        print("----------------Ma-----------------")
        print("k=",k)
        feature="d"+str(k)+"freLisTest"
        sim=[]
        for i in range(len(testset)):
            si=1/(dis.manhattan_feature(eval(feature)[i],eval(feature)[len(eval(feature))-1]))
            sim.append(si)  
#        print(sim)
        auc=roc_auc_score(testlabel, sim)
        print(auc)    
        
        print("----------------kld-----------------")
        print("k=",k)
        feature="d"+str(k)+"freLisTest"
        sim=[]
        for i in range(len(testset)):
            si=1/(dis.KLD_feature(eval(feature)[i],eval(feature)[len(eval(feature))-1]))
            sim.append(si)  
#        print(sim)
        auc=roc_auc_score(testlabel, sim)
        print(auc)   
        
        print("----------------pcc-----------------")
        print("k=",k)
        feature="d"+str(k)+"freLisTest"
        sim=[]
        for i in range(len(testset)):
            si=1/(dis.pcc_feature(eval(feature)[i],eval(feature)[len(eval(feature))-1]))
            sim.append(si)  
#        print(sim)
        auc=roc_auc_score(testlabel, sim)
        print(auc)   
        print("----------------cosine-----------------")
        print("k=",k)
        feature="d"+str(k)+"freLisTest"
        sim=[]
        for i in range(len(testset)):
            si=1/(dis.cosine_feature(eval(feature)[i],eval(feature)[len(eval(feature))-1]))
            sim.append(si)  
#        print(sim)
        auc=roc_auc_score(testlabel, sim)
        print(auc)   
        print("----------------che-----------------")
        print("k=",k)
        feature="d"+str(k)+"freLisTest"
        sim=[]
        for i in range(len(testset)):
            si=1/(dis.chebyshev_feature(eval(feature)[i],eval(feature)[len(eval(feature))-1]))
            sim.append(si)  
#        print(sim)
            
        auc=roc_auc_score(testlabel, sim)
        print(auc)   
        print("----------------D2-----------------")
        print("k=",k)
        feature="d"+str(k)+"freLisTest"
        sim=[]
        for i in range(len(testset)):
            si=1/(dis.getD2_feature(eval(feature)[i],eval(feature)[len(eval(feature))-1]))
            sim.append(si)  
#        print(sim)
        auc=roc_auc_score(testlabel, sim)
        print(auc)   
                


