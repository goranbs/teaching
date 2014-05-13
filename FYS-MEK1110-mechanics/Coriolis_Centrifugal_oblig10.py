# -*- coding: utf-8 -*-
"""
Created on Mon May  5 12:57:27 2014

@author: goranbs


FYS-MEK1110 - 2014 - Oblig 10 - Coriolis & Centrifugal force

           - vertical gunshot from earth's crust -
"""
from math import cos,sin,pi,sqrt
import numpy as np


m = 0.002                          # mass bullet
Rb = 0.0028                        # radius of bullet
A = pi*Rb**2                       # cross-section of bullet
rho = 1.225                        # [kg/m**3] density of air
Cv = 0.3                           # assumed drag coeff of bullet
D = 0.5*rho*A*Cv                   # total drag coeff for bullet

g = 9.83                           # acceleration of gravity measured at the north pole
phi = 54*pi/180                    # [rad] North angular position on earth's crust - latitude
theta = pi/2                       # [rad] - longitude
omega = 2*pi/((23*60+56)*60)       # angular velocity of earth
v0 = 200                           # initial velocity. v0 = vz
Re = 6380*10**3                    # radius of earth


n = 2000                           # number of timesteps
t0 = 0
tmax = 82.0
dt = float(tmax-t0)/n
time = np.linspace(t0,tmax,n)
v = np.zeros((n,3))
r = np.zeros((n,3))                # vector from origin on earth
R = np.zeros((n,3))                # vector from center of earth
a = np.zeros((n,3))
Fc = np.zeros((n,3))               # Coriolis
Fs = np.zeros((n,3))               # Centrifugal force
Fd = np.zeros((n,3))               # Drag force

v[0,:] = np.array([0,0,v0])        # initial velocity
R0 = Re*np.array([cos(phi)*cos(theta), cos(phi)*sin(theta), sin(phi)])

Time = 0
max_it = 0

for t in range(n-1):
   #                                  (vz*cos(phi)   - vy*sin(phi)) i  +  vx*sin(phi) j  +  vx*cos(phi) k 
    Fc[t,:] = -m*2*omega*np.array([(v[t,2]*cos(phi) - v[t,1]*sin(phi)) , v[t,0]*sin(phi) , v[t,0]*cos(phi) ])    
    Fs[t,:] = -m*(omega**2)*np.array([(-R[t,0]*cos(phi)**2 - sin(phi)*((R[t,2])*cos(phi) - R[t,1]*sin(phi)) ) , (sin(phi)*((R[t,2])*cos(phi) - r[t,1]*sin(phi))) , cos(phi)*(R[t,1]*sin(phi) - (R[t,2])*cos(phi)) ])    
    #Fs[t,:] = -(m*omega**2*Re)*np.array([0,sin(phi)*cos(phi), cos(phi)**2])
    vel = abs(sqrt(v[t,0]*v[t,0] + v[t,1]*v[t,1] + v[t,2]*v[t,2]))
    Fd[t,:] = -m*D*(vel**2)*(v[t,:]/vel)
    #Fd[t,:] = -m*D*(vel**2)*np.sign(v[t,:]) 
    G = np.array([0,0,-m*g])
    
    a[t,:] = (Fc[t,:] + G + Fs[t,:] + Fd[t,:])/m # including centrifugal & drag
    #a[t,:] = (Fc[t,:] + G + Fs[t,:])/m # including centrifugal 
    #a[t,:] = (Fc[t,:] + G)/m # without centrifugal 
        
    v[t+1,:] = v[t,:] + a[t,:]*dt
    r[t+1,:] = r[t,:] + v[t+1,:]*dt
    R[t+1,:] = R0 + r[t+1,:]
    
    # when does the bullet hit the ground?
    if r[t+1,2] <= 0 and r[t,2] > 0 :
        Time = time[t]
        max_it = t
    

a[-1, :] = a[-2,:]
Fc[-1,:] = Fc[-2,:]


print 'Time before the bullet hits the ground, t=%.2f' % Time

# plotting
import matplotlib.pyplot as plt

plt.figure()
plt.subplot(1,2,1)
plt.plot(r[0:max_it,0], r[0:max_it,2])
plt.xlabel('x [m]')
plt.ylabel('z [m]')
plt.title('Displacement x: %.3f' % r[max_it,0])

plt.subplot(1,2,2)
plt.plot(r[0:max_it,0], r[0:max_it,1])
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.title('Displacement y: %.3f' % r[max_it,1])

plt.show(True)


#print v[max_it/6,:]
