# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 09:11:11 2018

@author: 51164
"""
import Distance
class Similarity:
    def getEuSim(self,seqA,seqB,k):
        dis=Distance.Distance()
        sim=dis.EuD(seqA,seqB,k)
        return 1/sim
    def getmanhattanSim(self,seqA,seqB,k):
        dis=Distance.Distance()
        sim=dis.manhattan(seqA,seqB,k)
        return 1/sim
    def getKLDSim(self,seqA,seqB,k):
        dis=Distance.Distance()
        sim=dis.KLD(seqA,seqB,k)
        return 1/sim
    def getpccSim(self,seqA,seqB,k):
        dis=Distance.Distance()
        sim=dis.pcc(seqA,seqB,k)
        return 1/sim
    def getcosineSim(self,seqA,seqB,k):
        dis=Distance.Distance()
        sim=dis.cosine(seqA,seqB,k)
        return 1/sim
    def getchebyshevSim(self,seqA,seqB,k):
        dis=Distance.Distance()
        sim=dis.chebyshev(seqA,seqB,k)
        return 1/sim
    def getD2Sim(self,seqA,seqB,k):
        dis=Distance.Distance()
        sim=dis.getD2(seqA,seqB,k)
        return 1/sim
    
    
    def getMulD2WeightSim(self,seqA,seqB,kstart,kend,sequences,weight):
        dis=Distance.Distance()
        sim=dis.getMulD2Weight(seqA,seqB,kstart,kend,sequences,weight)
        return 1/sim
    def getMulD2WeightSim2(self,feaA,feaB,weight):
        dis=Distance.Distance()
        sim=dis.getMulD2Weight2(feaA,feaB,weight)
        return 1/sim
    
    
    
    def getD2WeightSim(self,seqA,seqB,k,sequences,weight):
        dis=Distance.Distance()
        sim=dis.getD2Weight(seqA,seqB,k,sequences,weight)
        return 1/sim
    
    
    def getD2sSim(self,seqA,seqB,k,r,flag,kmersetdic,kmer_pro):
        dis=Distance.Distance()
        sim=dis.getD2S(seqA,seqB,k,r,flag,kmersetdic,kmer_pro)
        return 1/sim
    
    def getD2sWeightSim(self,seqA,seqB,k,r,flag,kmersetdic,weight,kmer_pro):
        dis=Distance.Distance()
        sim=dis.getD2SWeight(seqA,seqB,k,r,flag,kmersetdic,weight,kmer_pro)
        return 1/sim
    
    def getMulD2sWeightSim(self,seqA,seqB,kstart,kend,r,flag,sequences,weight,kmer_pro):
        dis=Distance.Distance()
        sim=dis.getMulD2SWeight(seqA,seqB,kstart,kend,r,flag,sequences,weight,kmer_pro)
        return 1/sim
    
    def getMulD2sWeightSim2(self,feaA,feaB,weight):
        dis=Distance.Distance()
        sim=dis.getMulD2SWeight2(feaA,feaB,weight)
        return 1/sim
    
    
    def getD2starSim(self,seqA,seqB,k,r,flag,sequences,kmersetdic,kmer_pro):
        dis=Distance.Distance()
        sim=dis.getD2Star(seqA,seqB,k,r,flag,sequences,kmersetdic,kmer_pro)
        return 1/sim
    def getMulD2starWeightSim(self,seqA,seqB,kstart,kend,r,flag,sequences,weight,kmer_pro):
        dis=Distance.Distance()
        sim=dis.getMulD2StarWeight(seqA,seqB,kstart,kend,r,flag,sequences,weight,kmer_pro)
        return 1/sim
    def getMulD2starWeightSim2(self,feaA,feaB,weight):
        dis=Distance.Distance()
        sim=dis.getMulD2StarWeight2(feaA,feaB,weight)
        return 1/sim
    
    def getD2starWeightSim(self,seqA,seqB,k,r,flag,sequences,kmersetdic,weight,kmer_pro):
        dis=Distance.Distance()
        sim=dis.getD2StarWeight(seqA,seqB,k,r,flag,sequences,kmersetdic,weight,kmer_pro)
        return 1/sim
    

