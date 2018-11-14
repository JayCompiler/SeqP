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
            for file in filenames :  
                if os.path.splitext(file)[1] == '.fasta':  
                    L.append(os.path.join(dirpath, file)) 
                    fils.append(file)
        return L,fils
    
    ## 获取第一个数据集合
    def getData(self,file_dir):
        L,fils=self.file_name(file_dir)
        data=[]
        for file in L:
            with open(file) as a:
                lis=a.readlines()
                # 去掉换行符
                data.append(lis[1].strip('\n'))
        return data

if __name__=="__main__":
    rd=ReadData()
    data=rd.getData("dataset1")
    for da in data:
        print(da)
    