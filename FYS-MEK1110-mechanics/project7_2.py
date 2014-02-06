# Weather Balloon
#
# project 7.2 - FYS-MEK1110

from scitools.std import *
import numpy as np

##################################################################
# Constants

m = 0.1     # [kg]
g = 9.81    # [m/s^2]
Bouyancy = 11.3  # [m/s^2]
D = 0.1     # [kg/m]
A = 0.2     # m^2
eps = 0.001 # Default threshold
vt_x = 0    # terminal velocity in x direction
vt_z = 0    # terminal in z

##################################################################
# initial conditions
t0 = 0               # [s]
tmax = 10            # [s]
r0 = np.array([0,0]) # initial position
v0 = np.array([0,0]) # initial velocity
w0 = np.array([1,0]) # wind
d = 10               # some hight
N = 1000             # timesteps

dt = (tmax-t0)/float(N) # timestep
t =  zeros(N)
r = np.zeros((N,2))
v = np.zeros((N,2))
a = np.zeros((N,2))
w = np.zeros((N,2))
Fd = np.zeros((N,2))
G = np.zeros((N,2))
B = np.zeros((N,2))

#################################################################
#functions

def Drag(v,w):
    return D*sqrt((v-w)**2)*v

def wind(z):
    return w0[1]*(1 - exp(-z/d))

##################################################################
# main calculation

for i in range(N-1):
    w[i,0] = wind(r[i,1])
    G[i,1] = m*g
    B[i,1] = Bouyancy
    
    Fd[i,:] = Drag(v[i,:] ,w[i,:])
    a[i,:] = (B[i,:] - G[i,:] - Fd[i,:])/m
    v[i+1,:] = v[i,:] + a[i,:]*dt
    r[i+1,:] = r[i,:] + v[i+1]*dt
    #print Fd[i,:]

    if v[i,0] + eps > v[i+1,0]:
        vt_x = v[i+1,1];
    if v[i+1,1] < v[i,1] + eps:
        vt_z = v[i+1,1]

    t[i+1] = t[i] + dt
    


##################################################################
# plotting

import matplotlib.pyplot as plt

figure()                                     # figure #1
plt.subplot(3,1,1)
plt.plot(t,a[:,1])
plt.title('a,v and z in vertical direction')
plt.xlabel('time [s]')
plt.ylabel('acceleration [m/s^2]')
plt.legend('a(t)')

plt.subplot(3,1,2)
plt.plot(t,v[:,1])
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.legend('v(t)')

plt.subplot(3,1,3)
plt.plot(t,r[:,1])
plt.xlabel('time [s]')
plt.ylabel('position [m]')
plt.legend('z(t)')
plt.show(True)

figure()                                      # figure #2
plt.subplot(3,1,1)
plt.plot(t,a[:,0])
plt.title('a,v and x in in horizontal direction')
plt.xlabel('time [s]')
plt.ylabel('acceleration [m/s^2]')
plt.legend('a(t)')

plt.subplot(3,1,2)
plt.plot(t,v[:,0])
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.legend('v(t)')

plt.subplot(3,1,3)
plt.plot(t,r[:,0])
plt.xlabel('time [s]')
plt.ylabel('position [m]')
plt.legend('x(t)')
plt.show(True)

figure()                         # forceplot
plt.subplot(2,1,1)
plt.plot(t,Fd[:,0])
plt.title('Drag in horizontal direction')
#plt.xlabel('time [s]')
plt.ylabel('Force [N]')
plt.legend('Fd_x(t)')

plt.subplot(2,1,2)
plt.plot(t,Fd[:,1])
plt.hold('on')
plt.plot(t,B[:,1])
plt.plot(t,G[:,1])
plt.title('Forces in vertical direction')
plt.xlabel('time [s]')
plt.ylabel('Force [N]')
plt.legend('Fd_z(t)','G','B')
plt.show(True)

    



