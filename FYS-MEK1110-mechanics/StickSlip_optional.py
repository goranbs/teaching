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

random_number = random() # random floating point number in range (0,1)
print random_number

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

t0 = 0      # [s] start time
tmax = 2.0  # [s] stop time of simulation
n = 10000    # Number of timesteps
dt = (tmax-t0)/float(n) # timestep
print 'Timestep dt=%0.4f' %dt
x0 = 0


# matrices that holds the acceleration, velocity, position
# and forces in the coloumns.
a = np.zeros((M,n))
v = np.zeros((M,n))
x = np.zeros((M,n))
xb = np.zeros((M,n))
F = np.zeros((M,n))
R = np.zeros((M,n))
t = zeros(n)

# initial positions:
for i in range(M):
    x(i,0) = b*i + b*0.2*(1-2*random()))

##############################################
# function and for-loop:

def friction(F):
    if F >= my_s*N:
        #print 'my_d'
        return my_d*N
    else:
        #print 'my_s'
        return my_s*N

for t in range(n):
    for i in range(M):
        F[i,t] = k*(b*i -x[i,t]) - K*(x[i+1,t] - x[i,t] -b) - K*(x[i,t] - x[i-1,t] -b)
        R[i,t] = friction(F[i,t])
        
        a[i,t] = (F[i,t] -R[i,t]/m
        v[i+1,t] = v[i,t] + a[i,t]*dt
        x[i+1,t] = x[i,t] + v[i,t]*dt

        t(i+1) = t(i) + dt


for i in range(M):
    a[i,-1] = a[i,-2]

###################################################
#Plotting

