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
b = 0.1     # [m] equilibrium length F = 0
v0 = 0.0    # [m/s] initial velocity
u = 0.1     # [m/s] constant velocity of spring attachment
k = 100     # [N/m] spring coefficient
my_s = 0.0  # static friction
my_d = 0.0  # dynamic friction
N = m*g     # [N] Normal force
omega = sqrt(k/m) # [1/s] freqency

t0 = 0      # [s] start time
tmax = 2.0  # [s] stop time of simulation
n = 10000    # Number of timesteps
dt = (tmax-t0)/float(n) # timestep
print 'Timestep dt=%0.4f' %dt
x0 = 0.0

a = zeros(n)
v = zeros(n)
x = zeros(n)
xb = zeros(n)
F = zeros(n)
R = zeros(n)
t = zeros(n)

v[0] = v0
x[0] = x0

################################################
# the for-loop        

for i in range(n-1):
    xb[i] = b + u*t[i]
    F[i] = k*(xb[i]-x[i]-b)
    if (v[i] == 0):
        R[i] = -F[i]
        if (abs(R[i]) >= my_s*N):
            R[i] = -my_d*N
    else:
        R[i] = -my_d*N

    a[i] = (F[i] + R[i])/m
    v[i+1] = v[i] + a[i]*dt
    if ((v[i+1] < 0 and v[i] > 0) or (v[i+1] > 0 and v[i] < 0)):
        print "if statement is true, t=%.3f" % t[i]         
        v[i+1] = 0
    
    x[i+1] = x[i] + v[i+1]*dt
    t[i+1] = t[i] + dt


a[-1] = a[-2]
xb[-1] = xb[-2]

x_exact = u*t - (u/omega)*sin(omega*t)

###################################################
#Plotting
import matplotlib.pyplot as plt

plt.figure()
plt.plot(t[:],F[:],'r-o')
plt.hold(True)
plt.plot(t[:],-R[:],'b-*')
plt.legend(('F(t)', 'R(t)'))
plt.title('Forces. dt=%.4f' % dt )
plt.xlabel('time [s]')
plt.ylabel('Force [N]')
plt.hold(False)

plt.figure()
plt.subplot(3,1,1)
plt.plot(t[:],a[:])
plt.title('acceleration, velocity and position of box. dt=%.5f' % dt)
plt.ylabel('acceleration [m/s^2]')
plt.legend('a(t)')
plt.subplot(3,1,2)
plt.plot(t[:],v[:])
plt.ylabel('Velocity [m/s]')
plt.legend('v(t)')
plt.subplot(3,1,3)
plt.plot(t[:],x[:], 'b-')
plt.hold(True)
plt.plot(t,x_exact, 'r-')
plt.ylabel('position [m]')
plt.legend(('numerical','exact'))
plt.xlabel('time [s]')
plt.hold(False)

plt.figure()
plt.plot(t,x_exact)
plt.hold(True)
plt.plot(t,x)
plt.ylabel('position [m]')
plt.xlabel('time [s]')
plt.legend(('exact', 'numerical'))
plt.title('position of spring. u=%.2f' % u)
plt.hold(False)

plt.show(True)
