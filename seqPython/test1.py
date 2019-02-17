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




from sklearn.feature_selection import VarianceThreshold
X = [[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 1, 1], [0, 1, 0], [0, 1, 1]]
X=np.array(X)
sel = VarianceThreshold(threshold=(.8 * (1 - .8)))
X=sel.fit_transform(X)
print(X)

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
