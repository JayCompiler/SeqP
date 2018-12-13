# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 10:29:12 2018

@author: Yzi
"""

#import numpy as np
import os

class ReadData:
    ## 获取文件名列表 和路径列表
    def file_name(self,file_dir):   
        L=[]   
        fils=[]
        for dirpath, dirnames, filenames in os.walk(file_dir):  
            for file in filenames:  
                if os.path.splitext(file)[1] == '.fasta'or os.path.splitext(file)[1] == '.fa' or os.path.splitext(file)[1] == '.txt':  
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
    
    
    def getData2(self,file_dir):
        ## 获取fRS目录下的文件
        L,files=self.file_name("fRS")
        dataset=[]
        pos=[]
        neg=[]
        data=[]
        for file in L:
            if file_dir in file:
                with open(file) as a:
                    data=a.readlines()
        count=0
        while count <len(data):
            tmp=""
            ## 提取数据集，当开头不为>，或者count==len(data)-1
            while data[count][0]!=">" or count==len(data)-1:
                tmp=tmp+data[count].strip('\n')
                count=count+1
                if count==len(data):
                    break
            if tmp!="":
#                print(tmp)
                dataset.append(tmp)
            else:
                count=count+1
        ## 去除不必要的字符
        dataset=self.normdata(dataset)
        length=len(dataset)
        for i in range(int(length/2)):
            pos.append(dataset[i])
        for i in range(int(length/2),length):
            neg.append(dataset[i])
        return dataset,pos,neg
    
    
    ## 读取第三个文件的数据集和文件名
    def getData3(self,file_dir):
        ## 获取fRS目录下的文件
        L,files=self.file_name("e3")
        dataset=[]
#        pos=[]
#        neg=[]
        data=[]
        sequenceName=[]
        for file in L:
            if file_dir in file:
                with open(file) as a:
                    data=a.readlines()
        count=0
        while count <len(data):
            tmp=""
            ## 提取数据集，当开头不为>，或者count==len(data)-1
            if data[count][0]==">":
                sequenceName.append(data[count][1:].strip('\n'))
            while data[count][0]!=">" or count==len(data)-1:
                tmp=tmp+data[count].strip('\n')
                count=count+1
                if count==len(data):
                    break
            if tmp!="":
#                print(tmp)
                dataset.append(tmp)
            else:
                count=count+1
        ## 去除不必要的字符
        dataset=self.normdata(dataset)
        for i in range(len(sequenceName)):
            if len(sequenceName[i])>=10:
                sequenceName[i]=sequenceName[i][0:10]
            else:
                c=10-len(sequenceName[i])
                p=""
                for j in range(c):
                    p=p+" "
                sequenceName[i]=sequenceName[i]+p
        return dataset,sequenceName  
    
    
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
    data,pos,neg=rd.getData2("fly_blastoderm")
    print(pos)
#    print(len(L))
    print(len(pos))
    print(len(neg))
#    print(len(L))
#    L,files=rd.file_name("fRS")
#    print(L)
#    print(files)
#    print("ss1" in "sss")
#    print(rd.normdata(["ATCG","TACGASDJAS"]))
