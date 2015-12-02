# -*- coding: utf-8 -*-
"""
Created on Tue Dec 01 20:16:30 2015

@author: Stephen
"""

from random import *
from math import *
import matplotlib.pyplot as plt

h=0
j=-1
beta = 1000

s=16
n=s**2
mag = []
E_arr = []
mag_arr = []
mat = [[randrange(-1,2,2) for i in range(s)] for j in range(s)]

for x in range(100):
    for i in range(n):
        (sx, sy) = (randrange(s), randrange(s))
        de = (-h * -mat[sx][sy]) - (-h * mat[sx][sy])
        
        if sy-1 > -1: de += (-j * (-mat[sx][sy]) * mat[sx][sy-1]) - \
                            (-j * (mat[sx][sy]) * mat[sx][sy-1])
        if sy+1 < s: de += (-j * (-mat[sx][sy]) * mat[sx][sy+1]) - \
                            (-j * (mat[sx][sy]) * mat[sx][sy+1])
        if sx-1 > -1: de += (-j * (-mat[sx][sy]) * mat[sx-1][sy]) - \
                            (-j * (mat[sx][sy]) * mat[sx-1][sy])
        if sx+1 < s: de += (-j * (-mat[sx][sy]) * mat[sx+1][sy]) - \
                            (-j * (mat[sx][sy]) * mat[sx+1][sy])
        if de < 0 or random() < exp(-beta * de):  mat[sx][sy] = -mat[sx][sy]
    
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
            E_here += temp_E
    
    E_arr += [float(E_here)/(s**2)]
    mag_arr += [float(mag_here)/(s**2)]
    
    
    print "step " + str(x) + " done!"

plt.plot(E_arr)
plt.plot(mag_arr)   

plt.figure(1)
plt.subplot(211)
plt.plot(E_arr, 'b')

plt.subplot(212)
plt.plot(mag_arr, 'r')
plt.show()