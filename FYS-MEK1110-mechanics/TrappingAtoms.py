# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 21:23:34 2014

@author: goranbs
"""

from scitools.std import *
# Initialize parameters:
U0 = 150
m = 23
x0 = 2
alpha = 39.48
t = linspace(0,10,501)
dt = t[1] - t[0]
n = size(t)
x = zeros(n)
v = zeros(n)

# Initial Conditions
#x[0], v[0] = -5.0, 8.0
#x[0], v[0] = 0, 8.0
x[0], v[0] = 0, sqrt(4*U0/m)

for i in xrange(n-1):
    # Magnetic force M
    # Photon force F
    if abs(x[i]) >= x0:
        M,F = 0, 0
        
    else:
        M = -U0/x0*cmp(x[i],0)
        F = -alpha*v[i]
    
    a = (M+F)/m
    v[i+1] = v[i] + a*dt
    x[i+1] = x[i] + v[i+1]*dt
    
figure()
plot(t,x)
title('position of particle')
xlabel('time')
ylabel('position')

figure()
plot(t,v)
title('velocity of particle')
xlabel('time')
ylabel('velocity')