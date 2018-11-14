# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 14:42:05 2018

@author: Yzi
"""
#
#a=set([1,2,3])
#b=set([4,2,3])
#dic={1:2,2:2}
#c=[{"A":2,"F":4},{2:3,3:2}]
#
#print(c[0]["F"])
import numpy as np
import os
from sklearn import metrics

def file_name(file_dir):   
    L=[]   
    fils=[]
    for dirpath, dirnames, filenames in os.walk(file_dir):  
        for file in filenames :  
            if os.path.splitext(file)[1] == '.fasta':  
                L.append(os.path.join(dirpath, file)) 
                fils.append(file)
    
    return L,fils


if __name__=="__main__":
    L,fils=file_name("dataset1")
    data=[]
    label=[1,0,1,1,0,0]
    pred=[0.3,0.2,0.25,0.1,0.11,0.1]
    fpr, tpr, thresholds = metrics.roc_curve(label, pred,pos_label=1)
    auc=metrics.auc(fpr, tpr)
    print(auc)

    a=["ATACTASSGDGASDA","das","sda"]
    lis=[[1,2],[3,4],[5,8]]
    kk=np.array(lis)
    print(kk[1][1])
#    print(a[1:2])
#    print(b)