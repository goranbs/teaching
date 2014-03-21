# Newtons Cradle simulates a Newtons cradle.

import numpy as np
from scitools.std import *

def force(dx,d,k,q):
    if dx<d:
        F = k*abs(dx-d)**q
    else:
        F = 0.0
    return F

N = 3      # Number of balls
m = 1.0     # mass of balls
k = 100.0    # [N/m] 
#q = 4.0      # exponent
#q = 3.0/4   # exponent 
q = 1.0     # exponent
d = 0.1     # [m] distance between balls
v0 = 1.1    # [m/s] initial velocity of ball A
time = 1.0    # [s] simulation time
dt = 0.001  # [s] timestep
n = int(round(time/dt)) # number of timesteps
x = zeros((n,N),float) # position of every ball (0,1,...,N) at time ti = 0,1,2,...,n
v = x.copy()
t = zeros(n,float)
##########################################################
# Initial conditions

v[0,0] = v0       # initial velocity of ball 0

for i in range(N):
    x[0,1] = d*i # intial position of ball 0,1,...,N

for i in range(n-1):                         # Timeloop
    # Find force in vector F
    F = zeros(N,float)                       # update forces
    for j in range(1,N):                     # loop over balls 1:N = 1,2...,N
        #print 'loop range(1,N), j=%g' %j
        dx = x[i,j] - x[i,j-1]
        F[j] = F[j] + force(dx,d,k,q)
    for j in range(N-1):                     # loop over 
        dx = x[i,j+1] - x[i,j]
        F[j] = F[j] - force(dx,d,k,q)

    # Euler-Cromer vectorized step:

    a = F/m
    v[i+1] = v[i] + a*dt
    x[i+1] = x[i] + v[i+1]*dt
    t[i+1] = t[i] + dt
    

#########################################################
# plotting:

legends = []
figure()
for j in range(N):
    plot(t,v[:,j])
    temp = str(j)
    legends.append('Ball '+ temp)
    if j==0:
        hold('on')
    if j==N-1:
        hold('off')

title('Velocity. Number of balls N=%g' % N)
xlabel('time [s]')
ylabel('velocity [m/s]')
legend((legends), loc='center right')
#hardcopy('cradleOfNewton_ex_h',fontsize=25) # sets the fontsize in PostScript points only

legends = []
figure()
for j in range(N):
    plot(t,x[:,j])
    temp = str(j)
    legends.append('Ball '+ temp)
    if j==0:
        hold('on')
    if j==N-1:
        hold('off')

title('Position of balls. Number of balls N=%g' % N)
xlabel('time [s]')
ylabel('position [m]')
legend((legends), loc='center right')
#hardcopy('cradleOfNewton_ex_h',fontsize=25) # sets the fontsize in PostScript points only


print 'v/v0 = ' , v[n-1,:]/v0
