# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 14:49:55 2018

@author: Yzi
"""
import math

## 计算均值 时间复杂度为O(n)
def mean(dic):
    length=len(dic)
    su=0.0
    for key in dict.keys(dic):
        su=su+dic[key]
    mea=su/length
    return mea
## 计算方差，时间复杂度为O(n)
def var(dic):
    length=len(dic)-1
    su=0.0
    mea=mean(dic)
    for key in dict.keys(dic):
        su=su+(dic[key]-mea)**2
    return su/length

## 计算标准差 时间复杂度O(n)
def sd(dic):
    va=var(dic)
    return math.sqrt(va)