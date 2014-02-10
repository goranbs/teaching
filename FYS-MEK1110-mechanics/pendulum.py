#            Project 7.1 - Ball in a spring - FYS-MEK1110 - 09.02.12
#
# Aim of project is to model a pendulum consisting of a massless wire
# and attached to a solid ceiling and a ball in the other end. We 
# assume that the streching of the wire can be modeled from Hook's 
# law. We neglect airresistance and look at small angles theta.
#######################################################################

from scitools.std import *
import numpy as np

# Constants
m = 0.1       # [kg] mass of ball
k = 200       # [N/kg]
g = 9.81      # [m/s^2] acceleration of gravity
L0 = 1.0      # [m] lenght of wire, no force applied
N = 10000       # number of timesteps

# initial conditions
v0 = 6.0        # [m/s] initial speed 
theta0 = 30    # initial angle of wire, relative to vertical.
t0 = 0        
tmax = 10
dt = (float(tmax-t0)/float(N)) # timestep
print ' timestep used=%.4f' % dt
theta = np.zeros((N,1))
time = np.zeros((N,1))
r = np.zeros((N,2))
a = np.zeros((N,2))
v = np.zeros((N,2))
S = np.zeros((N,2))
G = np.zeros((N,2))

theta[0] = theta0*(180/pi)
r[0,0] = -L0*sin(theta[0])   # initial position
r[0,1] = -L0*cos(theta[0])
v[0,0] = v0*sin(theta[0])   # initial velocity
v[0,1] = v0*cos(theta[0])
#r[1,0] = 1
#r[1,1] = 0.5


########################################################################

def wire(radius):
    normR = linalg.norm(radius)
    if normR < L0:
        print 'rope tension is zero!'
        return 0
    else:
        return -k*(normR - L0)*radius

for i in range(N-1):
    G[i,1] = -m*g
    S[i,:] = wire(r[i,:])
    #print S[i,:]
    a[i,:] = (G[i,:] + S[i,:])/m
    v[i+1,:] = v[i,:] + a[i,:]*dt
    r[i+1,:] = r[i,:] + v[i+1,:]*dt

    R = linalg.norm(r[i+1,:])

    theta[i+1] = np.arcsin(r[i+1,0]/R)
    time[i+1] = time[i] + dt



#########################################################################
# plotting

import matplotlib.pyplot as plt


figure()
plt.subplot(3,1,1)
plt.plot(a[:,0],a[:,1])
plt.title('acceleration')
plt.xlabel('a_x [m/s^2]')
plt.ylabel('a_y [m/s^2]')
plt.legend('a(theta)')

plt.subplot(3,1,2)
plt.plot(v[:,0],v[:,1])
plt.title('velocity')
plt.xlabel('v_x [m/s]')
plt.ylabel('v_y [m/s]')
plt.legend('v(theta)')

plt.subplot(3,1,3)
plt.plot(r[:,0],r[:,1])
plt.title('position')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.legend('r(theta)')
plt.show(True)

