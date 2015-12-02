# -*- coding: utf-8 -*-
"""
Created on Wed Dec 02 00:22:43 2015

@author: Stephen
"""

from random import *
from math import *

import matplotlib.pyplot as plt

J = 1
H = 0
        
def initSys(size):
    sys = [[0 for k1 in xrange(size)] for k2 in xrange(size)]
    for i in xrange(size):
        for j in xrange(size):
            sys[i][j] = choice((-1,1))
            
    return sys
    
def randSwap(sys,size,beta):
    (sx, sy) = (randrange(size), randrange(size))
     
    dE = 0
     
    if sx > 0:
        dE += 2*J*sys[sx][sy]*sys[sx-1][sy]
    if sx < size-1:
        dE += 2*J*sys[sx][sy]*sys[sx+1][sy]
    if sy > 0:
        dE += 2*J*sys[sx][sy]*sys[sx][sy-1]
    if sy < size-1:
        dE += 2*J*sys[sx][sy]*sys[sx][sy+1]
         
    if random() < exp(-beta*dE):
        sys[sx][sy] = -sys[sx][sy]
         
    return sys

def calcM(sys,size):
    M = 0
    for i in xrange(size):
        for j in xrange(size):
            M += sys[i][j]
            
    return float(M)/(size**2)
    
def mcmcSweep(sys,size,beta):
    for steps in xrange(2500):
        sys = randSwap(sys,size,beta)
    
    return sys
    
            
if __name__ == '__main__':
    beta_c = log(1+sqrt(2))/(2*J)
    print beta_c
    
    size = 16
    beta = beta_c/2
    sys = initSys(size)
    M_arr = []
    for times in xrange(1000):
        sys = mcmcSweep(sys,size,beta)
        M_arr += [calcM(sys,size)]
    
    plt.plot(M_arr)
        
    