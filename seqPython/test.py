# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 14:42:05 2018

@author: Yzi
"""
import ReadData
import Sequence
import Similarity
from sklearn import metrics
import numpy as np

def ate1():
    print("a--------")


#dic ={"c":1,"b":4,"asd":3}
#for key in sorted(dic):
#    print(key,dic[key])
p1=[]
p2=[]
su=0
ma=0
mi=0
su1=0
ma1=0
mi1=0
for i in range(len(inner2)):
    su=inner2[i][0]+su
    su1=inner2[i][1]+su1
    if ma<inner2[i][0]:
        ma=inner2[i][0]
    if mi>inner2[i][0]:
        mi=inner2[i][0]
    if ma1<inner2[i][1]:
        ma1=inner2[i][1]
    if mi1>inner2[i][1]:
        mi1=inner2[i][1]
tmp1=[su/len(inner2),ma,mi]
tmp2=[su1/len(inner2),ma1,mi1]
p1.append(tmp1)
p2.append(tmp2)
    
su=0
ma=0
mi=0
su1=0
ma1=0
mi1=0
for i in range(len(inner3)):
    su=inner3[i][0]+su
    su1=inner3[i][1]+su1
    if ma<inner3[i][0]:
        ma=inner3[i][0]
    if mi>inner3[i][0]:
        mi=inner3[i][0]
    if ma1<inner3[i][1]:
        ma1=inner3[i][1]
    if mi1>inner3[i][1]:
        mi1=inner3[i][1]
tmp1=[su/len(inner3),ma,mi]
tmp2=[su1/len(inner3),ma1,mi1]
p1.append(tmp1)
p2.append(tmp2)
    
    
su=0
ma=0
mi=0
su1=0
ma1=0
mi1=0
for i in range(len(inner4)):
    su=inner4[i][0]+su
    su1=inner4[i][1]+su1
    if ma<inner4[i][0]:
        ma=inner4[i][0]
    if mi>inner4[i][0]:
        mi=inner4[i][0]
    if ma1<inner4[i][1]:
        ma1=inner4[i][1]
    if mi1>inner4[i][1]:
        mi1=inner4[i][1]
tmp1=[su/len(inner4),ma,mi]
tmp2=[su1/len(inner4),ma1,mi1]
p1.append(tmp1)
p2.append(tmp2)   


su=0
ma=0
mi=0
su1=0
ma1=0
mi1=0
for i in range(len(inner5)):
    su=inner5[i][0]+su
    su1=inner5[i][1]+su1
    if ma<inner5[i][0]:
        ma=inner5[i][0]
    if mi>inner5[i][0]:
        mi=inner5[i][0]
    if ma1<inner5[i][1]:
        ma1=inner5[i][1]
    if mi1>inner5[i][1]:
        mi1=inner5[i][1]
tmp1=[su/len(inner5),ma,mi]
tmp2=[su1/len(inner5),ma1,mi1]
p1.append(tmp1)
p2.append(tmp2)     
    

su=0
ma=0
mi=0
su1=0
ma1=0
mi1=0
for i in range(len(inner6)):
    su=inner6[i][0]+su
    su1=inner6[i][1]+su1
    if ma<inner6[i][0]:
        ma=inner6[i][0]
    if mi>inner6[i][0]:
        mi=inner6[i][0]
    if ma1<inner6[i][1]:
        ma1=inner6[i][1]
    if mi1>inner6[i][1]:
        mi1=inner6[i][1]
tmp1=[su/len(inner6),ma,mi]
tmp2=[su1/len(inner6),ma1,mi1]
p1.append(tmp1)
p2.append(tmp2)   


print(p1)
print(p2)  
    
    