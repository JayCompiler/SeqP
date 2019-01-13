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

## 返回列表的第2个 元素，即相似度
def getSim(Lis):
    return Lis[1]

if __name__=="__main__":
    rd=ReadData.ReadData()
    datasets,pos,neg=rd.getData2("fly_blastoderm")
#    datasets,pos,neg=rd.getData2("fly_pns")
    Sim=Similarity.Similarity()
    kstart=2
    kend=6
    sizePo=len(pos)
    
    flag=True
    top=int(sizePo*(sizePo-1)/2)
    print(top)

    a=[[1,2],[3,4]]
    b=np.array(a)
    print(b[:,1])
#    a=[1]
#    b=[2]
#    c=[6]
#    print(a+b+c)
#    print(a+b)
#    b=np.array(a)
#    print(b[:,0])






