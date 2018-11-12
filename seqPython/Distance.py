# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 16:04:51 2018

@author: Yzi
"""

import Sequence
import numpy as np
import math
class Distance:
    ## 欧式距离
    def EuD(self,seqA,seqB,k):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        kmerSet,dic =Sq.getSeqKerSet(seqLis,k)
        lisFea,freq=Sq.getSeqfreq(seqLis,k,dic)
        su=0.0
        for key in kmerSet:
            su=(lisFea[0][key]-lisFea[1][key])**2+su
#        print(lisFea,freq)
        return math.sqrt(su)
    ## 曼哈顿距离
    def manhattan(self,seqA,seqB,k):
        seqLis=[]
        # 变成list
        seqLis.append(seqA)
        seqLis.append(seqB)
        Sq =Sequence.Sequence()
        # 获取 关键字集合 字典dic
        kmerSet,dic =Sq.getSeqKerSet(seqLis,k)
        lisFea,freq=Sq.getSeqfreq(seqLis,k,dic)
        su=0.0
        for key in kmerSet:
            su=abs(lisFea[0][key]-lisFea[1][key])+su
#        print(lisFea,freq)
        return su
        
        
if __name__=="__main__":
    seqA="ATCG"
    seqB="TCAG"
    distance =Distance()
    su=distance.EuD(seqA,seqB,2)
    print(su)
    su=distance.manhattan(seqA,seqB,2)
    print(su)
        