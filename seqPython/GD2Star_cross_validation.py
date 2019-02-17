# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 21:44:59 2019
D2star 遗传算法
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
import markov

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
        for key in dict.keys(d2scountLis[0]):
            tmp_sim=tmp_sim+d2scountLis[trainLis[i][0]][key]*d2scountLis[trainLis[i][1]][key]*weight2M[key]*weight[0]
        for key in dict.keys(d3scountLis[0]):
            tmp_sim=tmp_sim+d3scountLis[trainLis[i][0]][key]*d3scountLis[trainLis[i][1]][key]*weight3M[key]*weight[1]
        for key in dict.keys(d4scountLis[0]):
            tmp_sim=tmp_sim+d4scountLis[trainLis[i][0]][key]*d4scountLis[trainLis[i][1]][key]*weight4M[key]*weight[2]   
        for key in dict.keys(d5scountLis[0]):
            tmp_sim=tmp_sim+d5scountLis[trainLis[i][0]][key]*d5scountLis[trainLis[i][1]][key]*weight5M[key]*weight[3]
        for key in dict.keys(d6scountLis[0]):
            tmp_sim=tmp_sim+d6scountLis[trainLis[i][0]][key]*d6scountLis[trainLis[i][1]][key]*weight6M[key]*weight[4] 
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
    CXPB, MUTPB, NGEN = 0.6, 0.3, 50
    
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
    name="human_muscle"
#    name="fly_blastoderm"
#    name="human_HBB"
    
    ## 获取数据集 整个数据集，正数据集，负数据集 都是序列，没有标签
    #datasets,pos,neg=rd.getData2("fly_blastoderm")
      ## 获取数据 kmer 从2-->6
    print("---GA_D2star------",name,"------------")
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
    
   ## 获取kmer 概率集合
    mar=markov.Markov()
    ## 获得概率
    kmer_pro2=mar.get_Mulk_Mul_kmer_Pro(datasets,2,2,0)
    kmer_pro3=mar.get_Mulk_Mul_kmer_Pro(datasets,3,3,0)
    kmer_pro4=mar.get_Mulk_Mul_kmer_Pro(datasets,4,4,0)
    kmer_pro5=mar.get_Mulk_Mul_kmer_Pro(datasets,5,5,0)
    kmer_pro6=mar.get_Mulk_Mul_kmer_Pro(datasets,6,6,0)

    ####### d2star特征--------------------
    d2scountLis=sq.getD2Star_Mul_seq_Count(datasets,2,0,True,d2dic,kmer_pro2)
    d3scountLis=sq.getD2Star_Mul_seq_Count(datasets,3,0,True,d3dic,kmer_pro3)
    d4scountLis=sq.getD2Star_Mul_seq_Count(datasets,4,0,True,d4dic,kmer_pro4)
    d5scountLis=sq.getD2Star_Mul_seq_Count(datasets,5,0,True,d5dic,kmer_pro5)    
    d6scountLis=sq.getD2Star_Mul_seq_Count(datasets,6,0,True,d6dic,kmer_pro6)
    
    for i in range(len(datasets)):
        d2scountLis[i]=sq.normdata_max_min(d2scountLis[i])
        d3scountLis[i]=sq.normdata_max_min(d3scountLis[i])
        d4scountLis[i]=sq.normdata_max_min(d4scountLis[i])
        d5scountLis[i]=sq.normdata_max_min(d5scountLis[i])
        d6scountLis[i]=sq.normdata_max_min(d6scountLis[i])  

    
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
       
        print("GA---------------第",o,"次验证---------------" )
        ## 归一化处理
        w=main()
    
        print("权重个数----",len(w))
        sim=[]
        ## 计算权重
      
        for i in range(len(testLis)):
            tmp_sim=0
            for key in dict.keys(d2scountLis[0]):
                tmp_sim=tmp_sim+d2scountLis[testLis[i][0]][key]*d2scountLis[testLis[i][1]][key]*weight2M[key]*w[0]
            for key in dict.keys(d3scountLis[0]):
                tmp_sim=tmp_sim+d3scountLis[testLis[i][0]][key]*d3scountLis[testLis[i][1]][key]*weight3M[key]*w[1]
            for key in dict.keys(d4scountLis[0]):
                tmp_sim=tmp_sim+d4scountLis[testLis[i][0]][key]*d4scountLis[testLis[i][1]][key]*weight4M[key]*w[2]
            for key in dict.keys(d5scountLis[0]):
                tmp_sim=tmp_sim+d5scountLis[testLis[i][0]][key]*d5scountLis[testLis[i][1]][key]*weight5M[key]*w[3]
            for key in dict.keys(d6scountLis[0]):
                tmp_sim=tmp_sim+d6scountLis[testLis[i][0]][key]*d6scountLis[testLis[i][1]][key]*weight6M[key]*w[4]
            sim.append(tmp_sim)     
        auc=roc_auc_score(testlabel, sim)    
        aucLis.append(auc)
        print("GA_D2star第",o,"次auc的值",auc)
    su=0
    for i in range(len(aucLis)):
        su=su+aucLis[i]
    avgAuc=su/num
    print("GA_D2star平均auc",avgAuc)
    end=time.process_time()
    print("程序时间：",(end-start))
    print("特征个数",len(d2scountLis[0])+len(d3scountLis[0])+len(d4scountLis[0])+len(d5scountLis[0])+len(d6scountLis[0]))
    
    
    
