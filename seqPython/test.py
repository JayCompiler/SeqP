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
#    print(L)
#    print(len(L))
#    print(fils)
    data=[]
    for file in L:
        with open(file) as a:
            lis=a.readlines()
            data.append(lis[1])
    for da in data:
        print(da)
    # os.walk 返回当前路径，文件夹，文件
#    rootpath,dirs,files=os.walk(".",topdown=False)
#    #for root, dirs, files in os.walk(".", topdown=False):
#    #    print(root)
#    #    for name in files:
#    #        print(os.path.join(root, name))
#    #    for name in dirs:
#    #        print(os.path.join(root, name))  
#    #ss=os.path.join(root,dirs)
#    print(rootpath)
#    print("-------------------------")
#    print(dirs)
#    print("-------------------------")
#    print(files)