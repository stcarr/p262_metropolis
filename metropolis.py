# -*- coding: utf-8 -*-
"""
Created on Tue Dec 01 20:16:30 2015

@author: Stephen
"""

from random import *
from math import *
import numpy
import matplotlib.pyplot as plt


if __name__ == '__main__':
    
    h=0
    j=1
    
    s_arr = [4,8,16]    
    
    beta_c = math.log(1+math.sqrt(2))/(2*j)
    print beta_c
    beta_arr = [beta_c/20,beta_c/15,beta_c/10,beta_c/5,beta_c/3,beta_c,1.5*beta_c,5*beta_c]

    
    E_sweep_arr = []
    mag_sweep_arr = []
    E_avg = [[]]
    M_avg = [[]]
    C_avg = [[]]
    
    start = 2000
    stop = 5000
    
    size_index = 0
    
    for s in s_arr:
            
        n=s**2
    
        for beta in beta_arr:
        
            mag = []
            E_arr = []
            E2_arr = []
            mag_arr = []
            mat = [[randrange(-1,2,2) for k1 in range(s)] for k2 in range(s)]
            
            for x in range(stop):
                for i in range(n):
                    
                    (sx, sy) = (randrange(s), randrange(s))
                    de = h*mat[sx][sy]
                    
                    if sy-1 > -1: 
                        de += 2*(j * (mat[sx][sy]) * mat[sx][sy-1])
                    if sy+1 < s: 
                        de += 2*(j * (mat[sx][sy]) * mat[sx][sy+1])
                    if sx-1 > -1: 
                        de += 2*(j * (mat[sx][sy]) * mat[sx-1][sy])
                    if sx+1 < s: 
                        de += 2*(j * (mat[sx][sy]) * mat[sx+1][sy])
                    if de < 0 or random() < exp(-beta * de):  
                        mat[sx][sy] = -mat[sx][sy]
                
                mag_here = 0
                E_here = 0
                
                for i in range(s):
                    for j in range(s):
                        
                        mag_here += mat[i][j]
                        
                        temp_E = 0
                        if i-1 > -1:
                            temp_E += -j * mat[i-1][j] * mat[i][j]
                        if i+1 < s:
                            temp_E += -j * mat[i+1][j] * mat[i][j]
                        if j-1 > -1:
                            temp_E += -j * mat[i][j-1] * mat[i][j]
                        if j+1 < s:
                            temp_E += -j * mat[i][j+1] * mat[i][j]
                        E_here += temp_E/2 - h*mat[i][j]
                
                E_arr += [float(E_here)]
                E2_arr += [float(E_here**2)]
                mag_arr += [float(mag_here)]
                
                
            E_sweep_arr += [E_arr]
            mag_sweep_arr += [mag_arr]
    
            E_avg[size_index] += [numpy.mean(E_arr[start:stop])/(s**2)]
            M_avg[size_index] += [numpy.mean(mag_arr[start:stop])/(s**2)]
            C_avg[size_index] += [beta**2*(numpy.mean(E2_arr[start:stop]) - (numpy.mean(E_arr[start:stop]))**2)/(s**2)]
            
            print "Beta = " + str(beta) + " done!"
        
        E_avg += [[]]
        M_avg += [[]]
        C_avg += [[]]
        size_index += 1
        print "----- System Size = " + str(s) + " done! ------"
    
    
    plt.figure(1)
    plt.subplot(311)
    for i in range(len(s_arr)):
        plt.plot(beta_arr,E_avg[i],label=s_arr[i])
    plt.legend()
    plt.subplot(312)
    for i in range(len(s_arr)):
        plt.plot(beta_arr,M_avg[i],label=s_arr[i])
    plt.legend()
    plt.subplot(313)
    for i in range(len(s_arr)):
        plt.plot(beta_arr,C_avg[i],label=s_arr[i])
    plt.legend()
    plt.show()
    
    '''
    plt.figure(1)
    plt.subplot(211)
    for i in range(len(E_sweep_arr)):
        plt.plot(E_sweep_arr[i], label=beta_arr[i])
        plt.legend()
    plt.subplot(212)
    for i in range(len(mag_sweep_arr)):
        plt.plot(mag_sweep_arr[i], label=beta_arr[i])
        plt.legend()
    plt.show()
    '''