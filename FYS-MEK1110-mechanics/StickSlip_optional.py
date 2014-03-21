# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 23:56:11 2014

@author: goranbs

**************************************************
     FYS-MEK1110 - Newtonian Mechanics 
     Project 9.1 - Mechanics - A.M.Sørensen
     Oblig #4 FYS-MEK1110 - spring 2014
     
            Stick-slip friction
    Model of a box being draged along a horizontal plane
    by a spring with constant velocity. The box has 
    stationary friction coefficient mu_s and dynamic
    friction coeff. mu_d. We asume that the spring can
    be modeled by Hook's law. And that the friction force
    is linked to the normal force by R = mu*N

"""
from scitools.std import *
from random import *
import numpy as np

m = 0.1     # [kg] weight of box
#g = 9.81    # [m/s^2] acceleration of gravity
b = 0.1     # [m] equilibrium length F = 0
v0 = 0.0    # [m/s] initial velocity
u = 0.1     # [m/s] constant velocity of spring attachment
k = 0.01     # [N/m] spring coefficient
K = 0.1     # [N/m] spring coefficient
my_s = 0.6  # static friction
my_d = 0.5  # dynamic friction
N = 1.0     # [N] Normal force
M = 10      # number of blocks
u = 0.1     # velocity of crust.
t0 = 0      # [s] start time
tmax = 1000.0  # [s] stop time of simulation
n = 10000    # Number of timesteps
dt = (tmax-t0)/float(n) # timestep

print 'Timestep dt=%0.4f' %dt

# matrices that holds the acceleration, velocity, position
# and forces in the coloumns.
a = np.zeros((n,M))
v = np.zeros((n,M))    # velocity of blocks
vb = np.zeros((n,M))   # velocity of attachmentpoints are constant
x = np.zeros((n,M))    # positions of blocks
xb = np.zeros((n,M))   # positions of attachment points
F = np.zeros((n,M))
R = np.zeros((n,M))
t = zeros(n)

vb[0,:] = u

# initial positions:
for i in range(M):
    x[0,i] = b*i + b*0.2*(1-2*random())
    xb[0,i] = x[i,0]
    
    
##############################################
# timeloop:

for j in range(n-1):         # timeloop

    absv = abs(v[j,i])

    for i in range(1,M-1):     # blockloop for blocks [1:M-1]
        F[j,i] = k*(xb[j,i] - x[j,i]) + K*(x[j,i+1] - x[j,i] -b) - K*(x[j,i] -x[j,i-1] -b)
        if absv == 0:
            R[j,i] = -F[j,i]
            if (abs(R[j,i]) >= my_s*N):
                R[j,i] = -my_d*N*np.sign(F[j,i])
        else:
            R[j,i] = -my_d*N*v[j,i]/absv
            

    # For the first block we only have a contribution to the force from the springs
    #  which are to the right. For the rightmost block we have only contributions 
    # from the springs which are to the left:
    F[j,0] = k*(xb[j,0] - x[j,0]) + K*(x[j,i+1] - x[j,i] - b)
    F[j,M-1] = k*(xb[j,0] - x[j,0]) - K*(x[j,i] - x[j,i-1] - b)
    
    if absv == 0:
        R[j,0] = -F[j,0]
        R[j,M-1] = -F[j,M-1]        
        if (abs(R[j,M-1]) >= my_s*N):
            R[j,M-1] = -my_d*N*np.sign(F[j,M-1])
        if (abs(R[j,0]) >= my_s*N):
            R[j,0] = -my_d*N*np.sign(F[j,0])
    else:
        R[j,M-1] = -my_d*N*v[j,M-1]/absv
        R[j,0] = -my_d*N*v[j,0]/absv


    # Euler-Chromer loop:
    
    for i in range(M):
        a[j,i] = (F[j,i] + R[j,i])/m
        v[j+1,i] = v[j,i] + a[j,i]*dt
        if (v[j+1,i] > 0 and v[j,i] < 0) or (v[j+1,i] < 0 and v[j,i] > 0):
            v[j+1,i] = 0

        x[j+1,i] = x[j,i] + v[j+1,i]*dt
        vb[j+1,i] = vb[j,i]
        xb[j+1,i] = xb[j,i] + vb[j+1,i]*dt
        t[j+1] = t[j] + dt

a[-1,:] = a[-2,:]
F[-1,:] = F[-2,:]
R[-1,:] = R[-2,:]

#########################################################
# Plotting

import matplotlib.pyplot as plt

plt.figure()
for i in range(M):
    plt.plot(t[:],x[:,i],'b-')
    plt.plot(t[:],xb[:,i],'r-')
    plt.hold(True)


plt.figure()
for i in range(M):
    plt.plot(t[:],F[:,i],'b-')
    plt.plot(t[:],R[:,i],'r-')
    plt.hold(True)

plt.show(True)
        
   
