# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 13:20:43 2019

@author: Yzi

使用 cnn进行 分类
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
import math

def add_layer(inputs,in_size,out_size,activation_function=None):
    Weights =tf.Variable(tf.random_normal([in_size,out_size]))
    biases=tf.Variable(tf.zeros([1,out_size]))+0.1
    Wx_plus_b=tf.matmul(inputs,Weights)+biases ## 矩阵乘法
    if activation_function==None:
        outputs=Wx_plus_b
    else:
        outputs=activation_function(Wx_plus_b)
    return outputs

def weight_variable(shape):
  initial = tf.truncated_normal(shape, stddev=0.1)       ## 产生随机变量
  return tf.Variable(initial)


def bias_variable(shape):
  initial = tf.constant(0.1, shape=shape)       ## 产生 bias 初始化
  return tf.Variable(initial)

def conv2d(x, W):
    #strides=[1,x_movement,y_movement,1]
  return tf.nn.conv2d(x, W, strides=[1, 1, 1, 1], padding='SAME')

def max_pool_2x2(x):
  return tf.nn.max_pool(x, ksize=[1, 2, 2, 1],
                        strides=[1, 2, 2, 1], padding='SAME')

## 输入为array数组，count表示扩展为几维，输出也为array  将数据扩展成图片格式
def tranToPic(data,dim):   
    size=len(data)
    length=len(data[0])
    ma=np.max(data,axis=0)
    mi=np.min(data,axis=0)
    width =(ma-mi)/dim
    XX=[]
    for i in range(len(data)):
        tmp=[]
        for j in range(len(data[0])):
            num=0
            x=[0]*dim
            value=math.floor((data[i,j]-mi[j])/width[j])
            while num<dim:
                if value==num:
                    x[num]=value/dim
#                    x[num]=1
                    tmp.append(x)
                    break
                if value==dim:
#                    x[num]=1
                    x[dim-1]=value/dim
                    tmp.append(x)
                    break                
                num=num+1
        XX.append(tmp)
    XX=np.array(XX).reshape([size,length,dim])
    YY=np.zeros([size,dim,length])
    for i in range(len(XX)):
        YY[i]=XX[i].T
    return YY


# 定义一个函数，按批次取数据

def minibatches(inputs=None, targets=None, batch_size=None, shuffle=False):
    assert len(inputs) == len(targets)
    if shuffle:
        indices = np.arange(len(inputs))
        np.random.shuffle(indices)
    for start_idx in range(0, len(inputs) - batch_size + 1, batch_size):
        if shuffle:
            excerpt = indices[start_idx:start_idx + batch_size]
        else:
            excerpt = slice(start_idx, start_idx + batch_size)
        yield inputs[excerpt], targets[excerpt]
   
    

if __name__=="__main__":
    rd=ReadData.ReadData()
#    name="muscle"
#    name="pns"
    name="fly_blastoderm"
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
    ## 归一化数据
#    pair = (pair-pair.min(axis=0)) / (pair.max(axis=0)-pair.min(axis=0))
   
        
    label=[]
    for i in range(len(posPair)):
        label.append(1)
    for i in range(len(negPair)):
        label.append(0)    
    label=np.array(label).reshape(len(pair),1)
    
    ## 特征选择
    #1)方差
    preNum=len(pair[0])
    sel = VarianceThreshold(threshold=(.8 * (1 - .8)))
    pair=sel.fit_transform(pair)
    nexNum=len(pair[0])
    print("方差特征选择降维：",preNum-nexNum)
    
        #2)单变量特征选择
    X_new = SelectKBest(chi2, k=784).fit_transform(pair, label)
    pair=X_new
    nex2Num=len(pair[0])
    print("单一变量特征选择降维：",nexNum-nex2Num)
    
    
    ##扩展成图像格式
#    dim=28
#    pair=tranToPic(pair,dim)
    print("-------------特征计算与转换 succcessfully----")
    
    feaNum=pair*nex2Num
    ## 划分训练集
    x_train,x_test,y_train,y_test = train_test_split(pair,label,test_size=0.2,random_state=0)
    print("训练集大小：",len(x_train))
    print("测试集大小",len(x_test))
   

    xs=tf.placeholder(tf.float32,[None,784]) #  图像特证数大小
    ys=tf.placeholder(tf.float32,[None,1])      #  二分类函数
    keep_prob=tf.placeholder(tf.float32)
    x_image=tf.reshape(xs,[-1,28,28,1])  ## -1表示不管样本数量，后面的1表示灰度图像

    ## 搭建cnn
    ##conv1 layer##
    W_conv1=weight_variable([5,5,1,32]) #patch 5x5,in size 1,out size 32
    b_conv1=bias_variable([32])
    
    h_conv1=tf.nn.relu(conv2d(x_image,W_conv1)+b_conv1) # output 28*28*32
    h_pool1=max_pool_2x2(h_conv1)   ## output size 14x14x32
    
    ##conv2 layer##
    W_conv2=weight_variable([5,5,32,64]) #patch 5x5,in size 32,out size 64
    b_conv2=bias_variable([64])
    h_conv2=tf.nn.relu(conv2d(h_pool1,W_conv2)+b_conv2) # output (dim/2)*(nexNum/2)*32
    h_pool2=max_pool_2x2(h_conv2)    ## output size (dim/4)*(nexNum/4)x64
    



    ##func1 layer##
    W_fc1=weight_variable([7*7*64,1024])
    b_fc1=bias_variable([1024])
    ## [n_sample,7,7,64]->>[n_sample,7*7*64]
    h_pool2_flat =tf.reshape(h_pool2,[-1,7*7*64])
    h_fc1=tf.nn.softplus(tf.matmul(h_pool2_flat,W_fc1)+b_fc1)
    h_fc1_drop=tf.nn.dropout(h_fc1,keep_prob)
    
    

    
    #func2 layer
    W_fc2=weight_variable([1024,1])
    b_fc2=bias_variable([1])
    prediction=tf.nn.sigmoid(tf.matmul(h_fc1_drop,W_fc2)+b_fc2) ## 仍然使用sigmoid函数
    
#    cross_entropy=tf.reduce_mean(-tf.reduce_sum(ys*tf.log(prediction),reduction_indices=[1]))
    loss =tf.reduce_mean(tf.reduce_sum(tf.square(ys-prediction),
                        reduction_indices=[1]))
    # 计算auc
    auc,op= tf.metrics.auc(ys, prediction)
    
#    train_step=tf.train.AdamOptimizer(1e-1).minimize(cross_entropy)
    train_step= tf.train.GradientDescentOptimizer(0.5).minimize(loss)
#    train_step=tf.train.AdagradDAOptimizer().minimize(loss)

    
#    train_step= tf.train.GradientDescentOptimizer(0.5).minimize(cross_entropy)
    init = tf.group(tf.global_variables_initializer(), tf.local_variables_initializer())

    sess=tf.Session()
    sess.run(init)
#    if int((tf.__version__).split('.')[1]) < 12 and int((tf.__version__).split('.')[0]) < 1:
#        init = tf.initialize_all_variables()
#    else:
#        init = tf.global_variables_initializer()
#    sess.run(init)    
    # 声明批量大小
    # 批量大小是指通过计算图一次传入多少训练数据
    batch_size = 50
    print("------训练开始-----")
    for i in range(1000): 
        rand_index = np.random.choice(len(x_train), size=batch_size)
        rand_x = x_train[rand_index]
        rand_y = y_train[rand_index]
        sess.run(train_step, feed_dict={xs: x_train, ys: y_train, keep_prob: 0.2})
        if i%50==0:
            print("损失:",sess.run(loss,feed_dict={xs:rand_x,ys:rand_y,keep_prob: 0.2}))
            prediction=sess.run(op,feed_dict={xs:rand_x,ys:rand_y,keep_prob: 0.2})
            prediction_value=sess.run(auc,feed_dict={xs:rand_x,ys:rand_y,keep_prob: 0.2})
            print("第",i,"次训练auc结果",prediction_value)
    print("------训练结束------")
    print("-----------------测试结果auc=：",sess.run(auc,feed_dict={xs:x_test,ys:y_test,keep_prob: 0.5}))

#
#    print("------训练开始-----")
#    for i in range(10000): 
#        sess.run(train_step, feed_dict={xs: x_train, ys: y_train, keep_prob: 0.5})
#        if i%50==0:
#            print("损失:",sess.run(loss,feed_dict={xs:x_train,ys:y_train,keep_prob: 0.5}))
#            prediction=sess.run(op,feed_dict={xs:x_train,ys:y_train,keep_prob: 0.5})
#            prediction_value=sess.run(auc,feed_dict={xs:x_train,ys:y_train,keep_prob: 0.5})
#            print("第",i,"次训练auc结果",prediction_value)
#    print("------训练结束------")
#    print("-----------------测试结果auc=：",sess.run(auc,feed_dict={xs:x_test,ys:y_test,keep_prob: 0.5}))
