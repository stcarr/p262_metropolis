# -*- coding: utf-8 -*-
"""
Created on Wed Dec 02 00:22:43 2015

@author: Stephen
"""

from random import *
from math import *

import matplotlib.pyplot as plt

import numpy
import time

J = 1
        
def initSys(size):
    sys = [[0 for k1 in xrange(size)] for k2 in xrange(size)]
    for i in xrange(size):
        for j in xrange(size):
            sys[i][j] = choice((-1,1))
            
    return sys
    
def randSwap(sys,size,beta,h):
    (sx, sy) = (randrange(size), randrange(size))
     
    dE = h*sys[sx][sy]
    
    if sx > 0:
        dE += 2*J*sys[sx][sy]*sys[sx-1][sy]
    else:
        dE += 2*J*sys[sx][sy]*sys[-1][sy]
        
    if sx < size-1:
        dE += 2*J*sys[sx][sy]*sys[sx+1][sy]
    else:
        dE += 2*J*sys[sx][sy]*sys[0][sy]
        
    if sy > 0:
        dE += 2*J*sys[sx][sy]*sys[sx][sy-1]
    else:
        dE += 2*J*sys[sx][sy]*sys[sx][-1]
        
    if sy < size-1:
        dE += 2*J*sys[sx][sy]*sys[sx][sy+1]
    else:
        dE += 2*J*sys[sx][sy]*sys[sx][0]
        
    ''' Non-Periodic
    if sx > 0:
        dE += 2*J*sys[sx][sy]*sys[sx-1][sy]
    if sx < size-1:
        dE += 2*J*sys[sx][sy]*sys[sx+1][sy]
    if sy > 0:
        dE += 2*J*sys[sx][sy]*sys[sx][sy-1]
    if sy < size-1:
        dE += 2*J*sys[sx][sy]*sys[sx][sy+1]
    '''
         
    if random() < exp(-beta*dE):
        sys[sx][sy] = -sys[sx][sy]
         
    return sys

def calcM(sys,size):
    M = 0
    for i in xrange(size):
        for j in xrange(size):
            M += sys[i][j]
            
    return float(M)/(size**2)
    
def mcmcSweep(sys,size,beta,h):
    for steps in xrange(size**2):
        sys = randSwap(sys,size,beta,h)
    
    return sys
    
def calcE(sys,size,beta,h):

    E = 0    
    
    for i in range(size):
        for j in range(size):

            E += -h*sys[i][j]            
            
            if i-1 > -1:
                E += -J * sys[i-1][j] * sys[i][j]
            else:
                E += -J * sys[-1][j] * sys[i][j]
                
            if i+1 < size:
                E += -J * sys[i+1][j] * sys[i][j]
            else:
                E += -J * sys[0][j] * sys[i][j]
                
            if j-1 > -1:
                E += -J * sys[i][j-1] * sys[i][j]
            else:
                E += -J * sys[i][-1] * sys[i][j]
            
            if j+1 < size:
                E += -J * sys[i][j+1] * sys[i][j]
            else:
                E += -J * sys[i][0] * sys[i][j]
                
    
    return E

def plotMCMC(E_arr,M_arr,size,beta):
    
    plt.figure(1)
    if size == 4:
        plt.subplot(411)
        plt.title('Initial Instability of Ising MCMC')
        plt.plot(M_arr,linewidth=3.0,label=str(size) + 'x' + str(size) + ' , beta = ' + str(beta))
        plt.ylabel('M')
        axes1 = plt.gca()
        axes1.set_ylim([-1.2,1.2])
        plt.legend()
    if size == 8:
        plt.subplot(412)
        plt.plot(M_arr,linewidth=3.0,label=str(size) + 'x' + str(size) + ' , beta = ' + str(beta))
        plt.ylabel('M')
        axes2 = plt.gca()
        axes2.set_ylim([-1.2,1.2])
        plt.legend()
    if size == 16:
        plt.subplot(413)
        plt.plot(M_arr,linewidth=3.0,label=str(size) + 'x' + str(size) + ' , beta = ' + str(beta))
        plt.ylabel('M')
        axes3 = plt.gca()
        axes3.set_ylim([-1.2,1.2])
        plt.legend()
    if size == 32:
        plt.subplot(414)
        plt.plot(M_arr,linewidth=3.0,label=str(size) + 'x' + str(size) + ' , beta = ' + str(beta))
        plt.ylabel('M')
        plt.xlabel('MCMC time')
        axes4 = plt.gca()
        axes4.set_ylim([-1.2,1.2])
        plt.legend()
    '''
    plt.subplot(212)
    plt.plot(M_arr,label=str(size) + 'x' + str(size) + ' , beta = ' + str(beta))
    plt.legend()
    '''
    plt.show()
    

def equilibTest():
    
    beta_c = log(1+sqrt(2))/(2*J)
    
    stop = 1000
    
    size_arr = [4,8,16,32]
    beta_arr = [beta_c/2,beta_c,beta_c*2]
    
    
    h = 0
    
    for size in size_arr:
        for beta in beta_arr:
            sys = initSys(size)
            
            M_arr = []
            E_arr = []
            
            for times in xrange(stop):
                sys = mcmcSweep(sys,size,beta,h)
                M_arr += [calcM(sys,size)]
                E = calcE(sys,size,beta,h)
                E_arr += [E]
        
            plotMCMC(E_arr,M_arr,size,beta)
    
    
def criticalSearch():
    beta_c = log(1+sqrt(2))/(2*J)
    h = 0
    
    start = 2000
    stop = 10000
    
    size_arr = [4,8,16,32]
    beta_arr = [beta_c/5,beta_c/3,beta_c/1.5,beta_c,beta_c*1.5,beta_c*2,beta_c*3,beta_c*4,beta_c*5,beta_c*10]

    E_avg = [[]]    
    M_avg = [[]]
    C_avg = [[]]
    
    size_index = 0
    
    for size in size_arr:
        
        
        for beta in beta_arr:
            
            sys = initSys(size)
            
            M_arr = []
            E_arr = []    
            E2_arr = []
            
            for times in xrange(stop):
                sys = mcmcSweep(sys,size,beta,h)
                M_arr += [calcM(sys,size)]
                E = calcE(sys,size,beta,h)
                E_arr += [E]
                E2_arr += [E**2]
                
            E_avg[size_index] += [numpy.mean(E_arr[start:stop])/(size**2)]
            M_avg[size_index] += [numpy.mean(M_arr[start:stop])]
            C_avg[size_index] += [beta**2*(numpy.mean(E2_arr[start:stop]) - (numpy.mean(E_arr[start:stop]))**2)/(size**2)]
            
            print "Beta = " + str(beta) + " done!"
        
        E_avg += [[]]
        M_avg += [[]]
        C_avg += [[]]
        
        size_index += 1
        
        print "----- System Size = " + str(size) + " done! ------" 
        
    t_arr = []
    for i in xrange(len(beta_arr)):
        t_arr += [1 - beta_c/beta_arr[i]]

    plt.figure(1)
    
    plt.subplot(311)
    plt.title('Critical behavior at different system sizes')
    for i in range(len(size_arr)):
        plt.plot(t_arr,E_avg[i],label=size_arr[i])
    plt.legend(loc=3)
    plt.ylabel('Avg E')
    axes1 = plt.gca()
    axes1.set_xlim([-1.5,1])

    plt.subplot(312)
    for i in range(len(size_arr)):
        plt.plot(t_arr,M_avg[i],label=size_arr[i])
    plt.legend(loc=3)
    plt.ylabel('Avg M')
    axes2 = plt.gca()
    axes2.set_xlim([-1.5,1])  
    axes2.set_ylim([-1.2,1.2])
    
    plt.subplot(313)
    for i in range(len(size_arr)):
        plt.plot(t_arr,C_avg[i],label=size_arr[i])
    plt.legend(loc=2)
    axes3 = plt.gca()
    axes3.set_xlim([-1.5,1])
    plt.ylabel('C')
    plt.xlabel('t')
    
    plt.show()           
    

def criticalExps():
    
    tic = time.time()
    
    size = 32
    start = 2000
    stop = 5000
    
    beta_c = log(1+sqrt(2))/(2*J)
    
    b_min = beta_c*1.05
    b_max = beta_c*1.25
    b_sample = 16
    
    h_min = 0.05
    h_max = 0.3
    h_sample = 16
    
    beta_arr = numpy.linspace(b_min,b_max,b_sample)
    H_arr = numpy.linspace(h_min,h_max,h_sample)
    
    M_avg_b = []
    M_avg_h = []
    
    for beta in beta_arr:
        
        h = 0 
            
        sys = initSys(size)
            
        M_arr = []
        
        for times in xrange(stop):
            sys = mcmcSweep(sys,size,beta,h)
            M_arr += [calcM(sys,size)]
            
        M_avg_b += [abs(numpy.mean(M_arr[start:stop]))]
        
        print "Beta = " + str(beta) + " done!"
    
    for h in H_arr:
        
        sys = initSys(size)
        
        beta = beta_c
       
        M_arr = []
        
        for times in xrange(stop):
            sys = mcmcSweep(sys,size,beta,h)
            M_arr += [calcM(sys,size)]
            
        M_avg_h += [numpy.mean(M_arr[start:stop])]
        
        print "H = " + str(h) + " done!"
    
    t_arr = []
    t_log = []
    M_b_log = []
    for i in xrange(len(beta_arr)):
        t_arr += [1 - beta_c/beta_arr[i]]
        t_log += [log(1 - beta_c/beta_arr[i])]
        M_b_log += [log(abs(M_avg_b[i]))]
    
    H_log = []
    M_h_log = []    
    for i in xrange(len(H_arr)):
        H_log += [log(abs(H_arr[i]))]
        M_h_log += [log(abs(M_avg_h[i]))]
        
    
    p1 = np.polyfit(t_log,M_b_log,1)
    p2 = np.polyfit(H_log,M_h_log,1)
    
    
    def fit1(x):
        return np.power(x,p1[0])
        
    def fit2(x):
        return np.power(x,p2[0])
    
    print "beta: " + str(p1[0])
    print "delta: " + str(1/p2[0])  
    
    plt.figure(1)
    plt.subplot(211)
    plt.plot(t_arr,M_avg_b, label="Simulation")
    plt.plot(t_arr,fit1(t_arr),label="Fit") 
    plt.title("M vs t, beta = " + str(p1[0]))
    plt.legend(loc=4)
    
    plt.subplot(212)
    plt.plot(H_arr,M_avg_h,label="Simulation")
    plt.plot(H_arr,fit2(H_arr),label="Fit")
    plt.title("M vs H, delta = " + str(1/p2[0]))
    plt.legend(loc=4)
    
    plt.show()
    
    return (M_avg_b,M_avg_h)
    
    print "Time Elapsed: " + str(time.time() - tic) + "s."
    

def finiteSizeScaling():
    beta_c = log(1+sqrt(2))/(2*J)
    h = 0
    
    start = 2000
    stop = 10000
    
    size_arr = [4,8,12]
    beta = beta_c

    E_avg = [] 
    M_avg = []
    C_avg = []
    X_avg = []
    
    for size in size_arr:
        
        sys = initSys(size)
        
        M_arr = []
        E_arr = []    
        E2_arr = []
        M2_arr = []
        
        for times in xrange(stop):
            sys = mcmcSweep(sys,size,beta,h)
            M = calcM(sys,size)*(size**2)
            E = calcE(sys,size,beta,h)
            E_arr += [E]
            E2_arr += [E**2]
            M_arr += [M]
            M2_arr += [M**2]
            
        E_avg += [numpy.mean(E_arr[start:stop])/(size**2)]
        M_avg += [numpy.mean(M_arr[start:stop])]
        C_avg += [beta**2*(numpy.mean(E2_arr[start:stop]) - (numpy.mean(E_arr[start:stop]))**2)/(size**2)]
        X_avg += [beta*(numpy.mean(M2_arr[start:stop]) - (numpy.mean(M_arr[start:stop]))**2)/(size**2)] 
        
        print "----- System Size = " + str(size) + " done! ------" 
        
    s_log = []
    C_log = []
    X_log = []   
        
    for i in xrange(len(size_arr)):
        s_log += [log(size_arr[i])]
        C_log += [log(C_avg[i])]
        X_log += [log(X_avg[i])]
    
    p1 = np.polyfit(s_log,C_log,1)
    p2 = np.polyfit(s_log,X_log,1)
    
    
    
    print "alpha/nu = " + str(p1[0])
    print "gamma/nu = " + str(p2[0])

    def fit1(x):
        return np.power(x,p1[0])*exp(p1[1])
        
    def fit2(x):
        return np.power(x,p2[0])*exp(p2[1])
    
    plt.figure(1)
    
    plt.subplot(211)
    plt.title('Finite Size Scaling at t = 0')
    plt.plot(size_arr,C_avg,label="Simulation")
    plt.plot(size_arr,fit1(size_arr),label="Fit (" + str(p1[0]) + ")") 
    plt.legend(loc=4)
    plt.ylabel('C/N')

    plt.subplot(212)
    plt.plot(size_arr,X_avg,label="Simulation")
    plt.plot(size_arr,fit2(size_arr),label="Fit (" + str(p2[0]) + ")") 
    plt.legend(loc=4)
    plt.ylabel('X/N')
    plt.xlabel('size')
    
    plt.show()
    
if __name__ == '__main__':
    
    #equilibTest()
    #criticalSearch()
    #(M_b_new,M_h_new) = criticalExps()
    #finiteSizeScaling()

    M_b = M_b_new
    M_h = M_h_new
    
    size = 32
    start = 2000
    stop = 5000
    
    beta_c = log(1+sqrt(2))/(2*J)
    
    b_min = beta_c*1.05
    b_max = beta_c*1.25
    b_sample = 16
    
    h_min = 0.05
    h_max = 0.3
    h_sample = 16
    
    beta_arr = numpy.linspace(b_min,b_max,b_sample)
    beta_arr = numpy.concatenate([beta_arr[:]])
    H_arr = numpy.linspace(h_min,h_max,h_sample)
    
    M_avg_b =numpy.concatenate([M_b[:]])
    M_avg_h = M_h    
    
    t_arr = []
    t_log = []
    M_b_log = []
    for i in xrange(len(beta_arr)):
        t_arr += [1 - beta_c/beta_arr[i]]
        t_log += [log(1 - beta_c/beta_arr[i])]
        M_b_log += [log(abs(M_avg_b[i]))]
    
    H_log = []
    M_h_log = []    
    for i in xrange(len(H_arr)):
        H_log += [log(abs(H_arr[i]))]
        M_h_log += [log(abs(M_avg_h[i]))]
        
    
    p1 = np.polyfit(t_log,M_b_log,1)
    p2 = np.polyfit(H_log,M_h_log,1)
    
    
    def fit1(x):
        return np.power(x,p1[0])*exp(p1[1])
        
    def fit2(x):
        return np.power(x,p2[0])*exp(p2[1])
    
    print "beta: " + str(p1[0])
    print "delta: " + str(1/p2[0])  
    
    plt.figure(1)
    plt.subplot(211)
    plt.plot(t_arr,M_avg_b, label="Simulation")
    plt.plot(t_arr,fit1(t_arr),label="Fit") 
    plt.title("M vs t, beta = " + str(p1[0]))
    plt.legend(loc=4)
    
    plt.subplot(212)
    plt.plot(H_arr,M_avg_h,label="Simulation")
    plt.plot(H_arr,fit2(H_arr),label="Fit")
    plt.title("M vs H, delta = " + str(1/p2[0]))
    plt.legend(loc=4)
    
    plt.show()    
