# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 10:29:12 2018

@author: Yzi
"""

import numpy as np
import os

class ReadData:
    ## 获取文件名列表 和路径列表
    def file_name(self,file_dir):   
        L=[]   
        fils=[]
        for dirpath, dirnames, filenames in os.walk(file_dir):  
            for file in filenames:  
                if os.path.splitext(file)[1] == '.fasta':  
                    L.append(os.path.join(dirpath, file)) 
                    fils.append(file)
        return L,fils
    
    ## 获取第一个数据集合
    def getData(self,file_dir):
        L,fils=self.file_name(file_dir)
        data=[]
        label=[]
        query=[]
        for file in L:
#            print(file.split('.')[0][-1])
            if file.split('.')[0][-1]!='-' and file.split('.')[0][-1]!='+':
                with open(file) as a:
                    lis=a.readlines()
                    # 去掉换行符
                    query.append(lis[1].strip('\n'))
            else:
                 with open(file) as a:
                    lis=a.readlines()
                    # 去掉换行符
                    data.append(lis[1].strip('\n'))
                    if file.split('.')[0][-1]=='-':
                        label.append(0)
                    else:
                        label.append(1)
        return self.normdata(query),self.normdata(data),label
    ## 去掉序列中的的单词
    def normdata(self,lis):
        newlis=[]
        start=0
        end=0
        for sequence in lis:
            tmp=[]
            for i in range(len(sequence)):
                if sequence[i]!="A" and sequence[i]!="T" and sequence[i]!="C"and \
                sequence[i]!="G":
                    end=i
                    tp=sequence[start:end]
                    tmp.append(tp)
                    start=end+1
                ## 处理末尾的情况
                if (sequence[i]=="A" or sequence[i]=="T" or sequence[i]=="C"or \
                sequence[i]=="G") and i==len(sequence)-1:
                    tmp.append(sequence[start:])
            newSeq=""
            for seq in tmp:
                newSeq=newSeq+seq
            newlis.append(newSeq)
            start=0
            end=0
        return newlis
        
                    

if __name__=="__main__":
    rd=ReadData()
    print(rd.normdata(["ATACTASSGDGASDA","TACGASDJAS"]))
#    query,data,label=rd.getData("dataset1")
#    print(query)
#    print(data)
#    print(label)
#    for da in data:
#        print(da)
#    