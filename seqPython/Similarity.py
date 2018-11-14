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
    def getD2sSim(self,seqA,seqB,k,r,flag,kmersetdic):
        dis=Distance.Distance()
        sim=dis.getD2S(seqA,seqB,k,r,flag,kmersetdic)
        return 1/sim
    def getD2starSim(self,seqA,seqB,k,r,flag,sequences,kmersetdic):
        dis=Distance.Distance()
        sim=dis.getD2Star(seqA,seqB,k,r,flag,sequences,kmersetdic)
        return 1/sim
