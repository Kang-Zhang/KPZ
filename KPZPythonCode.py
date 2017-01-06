# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import math
import time


def periodic(i, N):
    """Set periodic boundary conditions"""
    return (i % N + N) % N

def change_lattice(L_new, L, c_L, N, Dx, Dy, lx, ly, dt):
    """Returns new lattice configuration for each iteration"""
    for i in range (N): 
        for j in range (N):
            top = (periodic(i-1,N),j)
            right = (i,periodic(j+1,N))
            bottom = (periodic(i+1,N),j)
            left = (i, periodic(j-1,N))
            
            diffX = Dx*(np.sin(L[i,j] - L[right]) + np.sin(L[i,j] - L[left]))
            diffY = Dy*(np.sin(L[i,j] - L[top]) + np.sin(L[i,j] - L[bottom]))
            nonlinX = lx*((np.cos(L[i,j] - L[right] ) + np.cos(L[i,j] - L[left])) - 1)
            nonlinY = ly*((np.cos(L[i,j] - L[top]) + np.cos(L[i,j] - L[bottom])) - 1)
            noise = 2*math.pi*c_L*np.random.uniform(-0.5,0.5)
        
            L_step = dt*(diffX + diffY + nonlinX + nonlinY + noise)
            L_new[i,j] = L[i,j] - L_step
            
    L_new = L_new%(2*math.pi)   
    #L_new = abs(L_new)%(2*math.pi)             
    return L_new

def energy(L, N, Dx, Dy):
    """Returns the energy of each updated lattice"""
    E=0
    for i in range (N): 
        for j in range (N):
            top = (periodic(i-1,N),j)
            right = (i,periodic(j+1,N))
            bottom = (periodic(i+1,N),j)
            left = (i, periodic(j-1,N))

            E += -Dx*(np.cos(L[i,j] - L[right]) + np.cos(L[i,j] - L[left])) - Dy*(np.cos(L[i,j] - L[top]) + np.cos(L[i,j] - L[bottom]))
      
    return E
  
def vortices(L, N):
    """Returns the number of vortices of each updated lattice"""
    num = 0
    for i in range (N): 
        for j in range (N):
        #anticlockwise
            top = (periodic(i-1,N),j)
            top_left = (periodic(i-1,N), periodic(j-1,N))
            left = (i, periodic(j-1,N))
            
            d1 = L[top] - L[i,j]
            d2 = L[top_left] - L[top]
            d3 = L[left] - L[top_left]
            d4 = L[i,j] - L[left]
    
            num += (d1+d2+d3+d4)/2*math.pi
    return num

#Fastest way to do this? Vectorise?
def KPZ(N, Dx, Dy, lx, ly, c_L, dt, max_n):
    #initial configuration of NxN lattice, each site with initial phase between 0 and 2pi chosen using random number generator
    L = np.random.uniform(0, 2*math.pi, (N,N)) 
    #L = np.full((N, N), 0.5, dtype=float)
    convergence = 0
    n = 1 #number of iterations
    old_result = 0
    L_new = np.zeros((N,N))
    
    #iters = []
    #energy_density = []
    energy_density = np.zeros(max_n)
    while(convergence == 0):
        n = n + 1
        #print n-1
        
        L = change_lattice(L_new, L, c_L, N, Dx, Dy, lx, ly, dt) #calculate new lattice for current iteration
        #print L
        
        result = energy(L, N, Dx, Dy) #use new lattice to calculate new energy
        #print result
        
        energy_density[i] = result/N**2 #keep track of energy density values for each iteration
        
        #iters.append(n-1)
        #energy_density.append(result/N**2)
        
        #vor = vortices(L,N)
        #print vor
        
        #test for convergence
        change = (old_result - result).any()
        if change < 0:
            convergence = 1
        if (n > max_n):
            convergence = 1
            
        old_result = result 
        
    return energy_density
    
def avg_energy(N, Dx, Dy, lx, ly, c_L, dt, max_n, R):
    """Returns average energy density over a number of realisations R"""
    steady_stateE = np.zeros(R) #keep track of the last energy value for each realisation R
    x = np.arange(1,max_n+1)
    for i in range (R):
        y = KPZ(N, Dx, Dy, lx, ly, c_L, dt, max_n)
        #print x, y
        #plt.plot(x,y)
        #plt.title(r'Energy density over time for $c_L =  %s$'%(c_L) + '\n' + r'and lattice size $%s^2$ with %s realisations' %(N, R))
        #plt.xlabel('time [T]')
        #plt.ylabel(r'$E/N^2$')
        
        #save values for average energy density vs noise parameter plot
        steady_stateE[i] = y[-1]
    return np.mean(steady_stateE)


#start timing
time_start = time.clock()
###This takes too long - 34 minutes; plot in Python notebook instead
N = 64 #set lattice size
c_L = 0.5
R = 5
mean = avg_energy(N,1,1,0,0,c_L,0.05,100,R)

#stop time
time_elapsed = (time.clock() - time_start)
print 'time: ', time_elapsed
"""
#now 123 seconds using np.arrays
#134 seconds compared to 263 seconds doing each separately for c_L = 0.5 and 5 realisations
#261 seconds for c_L = 0.5 and 10 realisations
#so half an hour sounds about right for the whole thing

average of the last element of energy density from each realisation for each c_L
c_L = np.arange(0,5.5,0.5)
avg_steadyE = []
for i in c_L:
    print i
    plt.figure()
    if i <= 2.5:
        mean = avg_energy(N,1,1,0,0,i,0.05,100,5)
    else:
        mean = avg_energy(N,1,1,0,0,i,0.05,100,10)
    avg_steadyE.append(mean)
    

plt.figure()
plt.plot(c_L,avg_steadyE)
plt.title(r'Average energy density against noise parameter for lattice size $%s^2$' %(N))
plt.xlabel(r'$c_L$')
plt.ylabel(r'$<E>/N^2$')
"""

#plot derivative of average energy density against noise parameter (kind of like a Dirac delta)

