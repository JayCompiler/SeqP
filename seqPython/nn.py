# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 10:25:24 2019

@author: Yzi

使用nn 解决序列比对的分类问题


"""

## 多个k值的 以全部数据集求权重和特征
import ReadData
import Sequence
import Similarity
import time
import numpy as np
import tensorflow as tf
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
#用于数据划分
from sklearn.model_selection import train_test_split
# 用于特征选择
from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2

def add_layer(inputs,in_size,out_size,activation_function=None):
    Weights =tf.Variable(tf.random_normal([in_size,out_size]))
    biases=tf.Variable(tf.zeros([1,out_size]))+0.1
    Wx_plus_b=tf.matmul(inputs,Weights)+biases ## 矩阵乘法
    if activation_function==None:
        outputs=Wx_plus_b
    else:
        outputs=activation_function(Wx_plus_b)
    return outputs



if __name__=="__main__":
    rd=ReadData.ReadData()
#    name="muscle"
#    name="fly_blastoderm"
    name="pns"
    datasets,pos,neg=rd.getData2(name)
    print("----------------------数据集：",name,"-------------------")
    Sim=Similarity.Similarity()
    kstart=2
    kend=6
#    sizePo=len(pos)
#    top=int(sizePo*(sizePo-1))
    flag=True
    sq=Sequence.Sequence()
    ## 获得 count
    posD2DicLis=sq.getMulCount_nonorm(pos,kstart,kend,datasets) 
    negD2DicLis=sq.getMulCount_nonorm(neg,kstart,kend,datasets)
    ## 构造 序列对特征
    posD2Lis=[]
    negD2Lis=[]
    #####将特征与 理顺： 正权重
    ## 正特征
    for i in range(len(posD2DicLis)):
        tmp=[]
        for key in sorted(posD2DicLis[0]):
            tmp.append(posD2DicLis[i][key])
        posD2Lis.append(tmp)
    posArr=np.array(posD2Lis)  
    ## 负特征
    for i in range(len(negD2DicLis)):
        tmp=[]
        for key in sorted(posD2DicLis[0]):
            tmp.append(negD2DicLis[i][key])
        negD2Lis.append(tmp)
    negArr=np.array(negD2Lis)
    
    ## 构成基因对,并设置标识
    posPair=[]
    negPair=[]
    for i in range(len(posD2Lis)):
        for j in range(i+1,len(posD2Lis)):
            ## 添加绝对值
            tmpPosPair=np.abs(posArr[i]-posArr[j])
            tmpNegPair=np.abs(negArr[i]-negArr[j])
            posPair.append(tmpPosPair)
            negPair.append(tmpNegPair)
    posPair=np.array(posPair)
    negPair=np.array(negPair)
    pair=np.append(posPair,negPair,axis=0)
    averagePair=np.mean(pair,axis=0)
    for i in range(len(pair[0])):
        pair[:,i][pair[:,i]<averagePair[i]]=0
        pair[:,i][pair[:,i]>averagePair[i]]=1
    label=[]
    for i in range(len(posPair)):
        label.append(1)
    for i in range(len(negPair)):
        label.append(0)
    label=np.array(label).reshape(len(pair),1)
    print("-------------特征计算 succcessfully----")
    ## 特征选择
    #1)方差
    preNum=len(pair[0])
    sel = VarianceThreshold(threshold=(.8 * (1 - .8)))
    pair=sel.fit_transform(pair)
    nexNum=len(pair[0])
    print("方差特征选择降维：",preNum-nexNum)
    #2)单变量特征选择
    X_new = SelectKBest(chi2, k=200).fit_transform(pair, label)
    pair=X_new
    nex2Num=len(pair[0])
    print("单一变量特征选择降维：",nexNum-nex2Num)
    

    ## 划分训练集
    x_train,x_test,y_train,y_test = train_test_split(pair,label,test_size=0.2,random_state=0)

    ### 训练过程

    
    x_row=len(x_train)
    x_col=len(x_train[0])
    
    
    xs=tf.placeholder(tf.float32,[None,x_col])
    ys=tf.placeholder(tf.float32,[None,1])
    
    
    
    ## 隐含层
    l1=add_layer(xs,len(x_train[0]),300,activation_function=tf.nn.softplus)
    ## 输出层
    prediction =add_layer(l1,300,1,activation_function=tf.nn.sigmoid)
    
    loss =tf.reduce_mean(tf.reduce_sum(tf.square(ys-prediction),
                        reduction_indices=[1]))
    
    ## 计算auc
    auc,op= tf.metrics.auc(y_train,prediction)

    train_step= tf.train.GradientDescentOptimizer(0.5).minimize(loss)
    
    init = tf.group(tf.global_variables_initializer(), tf.local_variables_initializer())
#    sess.run(init)
#    init =tf.initialize_all_variables()
    
    sess=tf.Session()
    sess.run(init)
    
    
    for i in range(1000):
        sess.run(train_step,feed_dict={xs:x_train,ys:y_train})
        if i%20==0:
            prediction=sess.run(op,feed_dict={xs:x_train,ys:y_train})
            prediction_value=sess.run(auc,feed_dict={xs:x_train,ys:y_train})
            print("第",i,"次训练auc结果",prediction_value)

    print("-----------------测试结果auc=：",sess.run(auc,feed_dict={xs:x_test,ys:y_test}))
    
    

