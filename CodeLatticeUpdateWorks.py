# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math
import time

def periodic(i, N):
    return (i % N + N) % N

def change_lattice(L_new, L, c_L, N, Dx, Dy, lx, ly, dt):
    for i in range (N): 
        for j in range (N):
            """Set periodic boundary conditions"""
            if i == 0: #first row
                top = (N-1,j) #sum using last row
            else:
                top = (i-1,j)

            if j == N-1: #last column
                right = (i,0) #sum using first column
            else:
                right = (i,j+1)

            if i == N-1: #last row
                bottom = (0,j) #sum using first row
            else:
                bottom = (i+1,j)

            if j == 0: #first column
                left = (i, N-1) #sum using last column
            else:
                left = (i, j-1)
    
            #do the x and y diff
            top  = L[i,j] - L[top]
            bottom = L[i,j] - L[bottom]
            right = L[i,j] - L[right]        
            left = L[i,j] - L[left]
            
            diffX = Dx*np.sum(np.sin(right) + np.sin(left))
            diffY = Dy*np.sum(np.sin(top) + np.sin(bottom))
            nonlinX = lx*(np.sum(np.cos(right) + np.cos(left)) - 1)
            nonlinY = ly*(np.sum(np.cos(top) + np.cos(bottom)) - 1)
            noise = 2*math.pi*c_L*np.random.uniform(-0.5,0.5)
        
            L_step = dt*(diffX + diffY + nonlinX + nonlinY + noise)
            L_new[i,j] = L[i,j] - L_step
            
    L_new = abs(L_new)%(2*math.pi)                
    return L_new

def energy(L, N, Dx, Dy):
    E = 0
    for i in range (N): 
        for j in range (N):
            """Set periodic boundary conditions"""
            if i == 0: #first row
                top = (N-1,j) #sum using last row
            else:
                top = (i-1,j)

            if j == N-1: #last column
                right = (i,0) #sum using first column
            else:
                right = (i,j+1)

            if i == N-1: #last row
                bottom = (0,j) #sum using first row
            else:
                bottom = (i+1,j)

            if j == 0: #first column
                left = (i, N-1) #sum using last column
            else:
                left = (i, j-1)

            E += -Dx*np.sum(np.cos(L[i,j] - L[right]) + np.cos(L[i,j] - L[left])) - Dy*np.sum(np.cos(L[i,j] - L[top]) + np.cos(L[i,j] - L[bottom]))
      
    return E

#Fastest way to do this? Vectorise?
def KPZ(N, Dx, Dy, lx, ly, c_L, dt, max_n):
    #initial configuration of NxN lattice, each site with initial phase between 0 and 2pi chosen using random number generator
    L = np.random.uniform(0, 2*math.pi, (N,N)) 
    #L = np.full((N, N), 0.5, dtype=float)
    convergence = 0
    n = 1 #number of iterations
    old_result = 0
    L_new = np.zeros((N,N))
    
    iters = []
    energy_density = []
    while(convergence == 0):
        n = n + 1
        #print n-1
        
        L = change_lattice(L_new, L, c_L, N, Dx, Dy, lx, ly, dt)
        #print L
        result = energy(L, N, Dx, Dy)
        #print result
        
        iters.append(n-1)
        energy_density.append(result/N**2)
        
        change = (old_result - result).any()
        if change < 0:
            convergence = 1
        if (n > max_n):
            convergence = 1
            
        old_result = result 
    
    return iters, energy_density

#start timing
time_start = time.clock()

#print KPZ(2, 1, 1, 1, 1, 0, 0.05, 10)


#print x, y
x,y = (KPZ(64, 1, 1, 0, 0, 0.5, 0.05, 100))
plt.plot(x,y)

#stop time
time_elapsed = (time.clock() - time_start)
print 'time: ', time_elapsed