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

m = 0.1     # [kg] weight of box
g = 9.81    # [m/s^2] acceleration of gravity
b = 0.0     # [m] equilibrium length F = 0
v0 = 0.1    # [m/s] initial velocity
u = 0.01     # [m/s] constant velocity of spring attachment
k = 200     # [N/m] spring coefficient
my_s = 0.7  # static friction
my_d = 0.3  # dynamic friction
N = m*g     # [N] Normal force

t0 = 0      # [s] start time
tmax = 2.0  # [s] stop time of simulation
n = 1000    # Number of timesteps
dt = (tmax-t0)/float(n) # timestep
print 'Timestep dt=%0.4f' %dt
x0 = 0

a = zeros(n)
v = zeros(n)
x = zeros(n)
xb = zeros(n)
F = zeros(n)
R = zeros(n)
t = zeros(n)

v[0] = v0
x[0] = x0

def friction(F):
    if F >= my_s*N:
        print 'my_d'
        return my_d*N
    else:
        print 'my_d'
        return my_s*N

def spring(x,t):
    return (k/m)*(Xb(t) - x - b)

def Xb(t):
    return u*t
        
for i in range(n-1):
    F[i] = spring(x[i],t[i])
    R[i] = friction(F[i])
    a[i] = (F[i] - R[i])/m
    v[i+1] = v[i] + a[i]*dt
    xb[i] = Xb(t[i])
    x[i+1] = x[i] + v[i]*dt # xb[i]
    t[i+1] = t[i] + dt


a[-1] = a[-2]
xb[-1] = xb[-2]


###################################################
#Plotting
import matplotlib.pyplot as plt

plt.figure()
plt.plot(t[:],F[:],'r-o')
hold('on')
plt.plot(t[:],R[:],'b-*')
plt.legend('F(t)', 'R(t)')
plt.title('Forces')
plt.xlabel('time [s]')
plt.ylabel('Force [N]')
hold('off')

plt.figure()
plt.subplot(3,1,1)
plt.plot(t[:],a[:])
plt.title('acceleration, velocity and position of box')
plt.ylabel('acceleration [m/s^2]')
plt.legend('a(t)')
plt.subplot(3,1,2)
plt.plot(t[:],v[:])
plt.ylabel('Velocity [m/s]')
plt.legend('v(t)')
plt.subplot(3,1,3)
plt.plot(t[:],x[:])
plt.ylabel('position [m]')
plt.legend('x(t)')
plt.xlabel('time [s]')


plt.figure()
plt.plot(t,xb)
plt.ylabel('position [m]')
plt.xlabel('time [s]')
plt.legend('xb(t)')
plt.title('position of spring')

plt.show(True)
