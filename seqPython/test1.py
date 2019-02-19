# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 09:13:50 2018

@author: Yzi
"""

## 多个k值的 以全部数据集求权重和特征
import ReadData
import Sequence
import Similarity
import time
import markov
import numpy as np
import math
from sklearn.feature_selection import VarianceThreshold



#print(np.random.choice(5, 3))

pair=[[1,2],[3,4],[5,6],[7,8],[9,10],[11,12]]
label=[1,1,1,0,0,0]

pair=np.array(pair)
label=np.array(label)
print(a1[1:2,:])
x_train_pos,x_test_pos,y_train_pos,y_test = train_test_split(pair,label,test_size=0.2,random_state=0)

#a=[0,4,5]
#print(a1[a])

#x=[0]*4
#print(x)

#X=np.array(X)
#ma=np.max(X,axis=0)
#mi=np.min(X,axis=0)
#count=3
#width =(ma-mi)/count
#print(width)
#print(mi)
#XX=[]
#for i in range(len(X)):
#    tmp=[]
#    for j in range(len(X[0])):
#        num=0
#        x=[0]*count
#        while num<count:
##            print(math.floor((X[i,j]-mi[j])/width[j]))
#            if math.floor((X[i,j]-mi[j])/width[j])==num:
#                x[num]=1
#                tmp.append(x)
#                break
#            if math.floor((X[i,j]-mi[j])/width[j])==count:
#                x[count-1]=1
#                tmp.append(x)
#                break                
#            num=num+1
#    XX.append(tmp)
#XX=np.array(XX).reshape([len(X),len(X[0]),count])
#print(XX)
#print("----------------")
#YY=np.zeros([len(X),count,len(X[0])])
#for i in range(len(XX)):
##    print(XX[i].T)
#    YY[i]=XX[i].T
#print(YY)


#def tranToPic(data,dim):   
#    size=len(data)
#    length=len(data[0])
#    ma=np.max(data,axis=0)
#    mi=np.min(data,axis=0)
#    width =(ma-mi)/dim
#    XX=[]
#    for i in range(len(data)):
#        tmp=[]
#        for j in range(len(data[0])):
#            num=0
#            x=[0]*dim
#            value=math.floor((data[i,j]-mi[j])/width[j])
#            while num<dim:
#                if value==num:
#                    x[num]=value/dim
#                    tmp.append(x)
#                    break
#                if value==dim:
#                    x[dim-1]=value/dim
#                    tmp.append(x)
#                    break                
#                num=num+1
#        XX.append(tmp)
#    XX=np.array(XX).reshape([size,length,dim])
#    YY=np.zeros([size,dim,length])
#    for i in range(len(XX)):
#        YY[i]=XX[i].T
#    return YY
#
#X = [[1, 9, 1], [3,6,5], [4, 3, 7], [6, 5, 9], [10, 12, 5], [7, 1, 8]]
#print(tranToPic(np.array(X),3))

#X = [[1, 9, 1], [3,6,5], [4, 3, 7], [6, 5, 9], [10, 12, 5], [7, 1, 8]]
#X=np.array(X)
#ma=np.max(X,axis=0)
#mi=np.min(X,axis=0)
#width =(ma-mi)/2
#print(width)
#
#XX=[]
#for i in range(len(X)):
#    tmp=[]
#    for j in range(len(X[0])):
##        if math.floor((X-ma[j])/width[j])=0
#        x=[1,0,0,0]
#        tmp.append(x)
#    XX.append(tmp)
#XX=np.array(XX).reshape([6,3,4])
#print(XX)
#print("----------------")
#YY=np.zeros([6,4,3])
#for i in range(len(XX)):
##    print(XX[i].T)
#    YY[i]=XX[i].T
#print(YY)



#X=np.array(X)
#X[2,2]=np.array([0,1,3])
#print(X)
#sel = VarianceThreshold(threshold=(.8 * (1 - .8)))
#X=sel.fit_transform(X)
#print(X)

#if __name__=="__main__":
#    rd=ReadData.ReadData()
#    datasets,pos,neg=rd.getData2("fly_blastoderm")
##    datasets,pos,neg=rd.getData2("fly_pns")
#    Sim=Similarity.Similarity()
#    kstart=2
#    kend=6
#    sizePo=len(pos)
#    
#    flag=True
#    top=int(sizePo*(sizePo-1)/2)
#    print(top)
#    
##    a=[1,2]
##    b=[3,4]
##    aArr=np.array(a)
##    bArr=np.array(b)
##    print(aArr[1])
##    print(np.abs(aArr-bArr))
#    
#    a1=[1,2]
#    a2=[3,4]
#    bb=[[1,2],[3,4]]
#    
#    print(bb[:,0])
#    
#    print(a1+a2)
#    
#    a = np.array([[1, 2], [3, 4]])
#    me=np.mean(a, axis=0)
#    print(me)
##    print(np.mean(a, axis=1))
#    a[:,0][a[:,0]<me[0]]=0
#    a[:,0][a[:,0]>me[0]]=1
#    print(a)
#
#
#
#
#
#
#
#
#
