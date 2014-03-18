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
#x[0], v[0] = -5.0, 10.0
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


###########################################################
# Find largest initial velocity for the atom not to escape:

noEscape = []
x1 = zeros(n)
v1 = zeros(n)
x1[0] = -5
# Need to have a closer look at this part:
for v0 in linspace(8,10,10001):
    v1[0] = v0
    for i in xrange(n-1):
        # Magnetic force M
        # Photon force F
        if abs(x1[i]) >= x0:
            M,F = 0, 0
            
        else:
            M = -U0/x0*cmp(x1[i],0)
            F = -alpha*v1[i]
                
        a = (M+F)/m
        v1[i+1] = v1[i] + a*dt
        x1[i+1] = x1[i] + v1[i+1]*dt
                
        if x1[i+1] > x0:
            break
        
    noEscape.append(v0)
     
        
print 'For initial position x[0]=%.2f . \n Largest initial velocity that the atom can have \n and still not escape is v0=%.3f ' % (x1[0],noEscape[-1])

