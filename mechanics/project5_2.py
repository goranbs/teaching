# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 12:58:24 2013

@author: goranbs

Introduction to Mehanics - FYS-MEK1110

Project 5.2
            Realistic model for a 100m race.
"""
from scitools.std import *

# b) assume constant applied force F = 400N

F = 400       # Applied force [N]
m = 80        # mass of sprinter [kg]
rho = 1.293   # [kg/m^3] density of air
Cd = 1.2      # Drag coefficient 
A = 0.45      # [m^2] front area of runner
w = 0         # [m/s] wind
fv = 25.8     # [Ns/m] Driving force coefficient 
tc = 0.67     # [s] characteristic time
fc = 488      # [N] initial driving force in crouched position

def D(v,t):
    #print Area(t)
    return 0.5*rho*Cd*Area(t)*(v-w)**2

def Fv(v):
    return fv*v

def FC(t):
    return fc*exp(-(t/tc)**2)
    
def Area(t):
    return A*(1-0.25*exp(-(t/tc)**2))

t0 = 0                                  # start time
t_max = 10                              # stop time of simulation
timesteps = 100                         # number of timesteps
dt = (t_max - t0)/float(timesteps)      # timestep

t = zeros(timesteps)
v = zeros(timesteps)                # velocity array
x = zeros(timesteps)                # position array
Drag = zeros(timesteps)             # Dragforce
Acc = zeros(timesteps)              # Acceleration
Fc = zeros(timesteps)               # initial driving force
Fd = zeros(timesteps)               # Driving force
F_v = zeros(timesteps)              # Decreasing driving force 

eps = 0.006                         # error
goal = 10000                        # finsh time
vt = 0                              # terminal velocity
# initial conditions:
v[0] = 0                         
x[0] = 0
Fc[0] = fc
Fd[0] = F + Fc[0]- fv*v[0]

for n in range(len(t)-1):
    Acc[n] = (Fd[n] - Drag[n])/m     # Acceleration 
    v[n+1] = v[n] + dt*Acc[n]        # Forward Euler
    x[n+1] = x[n] + dt*v[n+1]        # Backward Euler
    t[n+1] = t[n] + dt
    Drag[n+1] = D(v[n+1],t[n+1])     # Drag force
    Fc[n+1] = FC(t[n+1])             # initial driving force
    F_v[n+1] = Fv(v[n+1])            #
    Fd[n+1] = F + Fc[n+1] - F_v[n+1] # Total driving force   
    if v[n+1] <= v[n] + eps and v[n+1] >= v[n] - eps:
        vt = v[n+1]
    if x[n+1] < 101 and x[n+1] > 99:
        goal = t[n+1]
        
Acc[-1] = (Fd[-1] - Drag[-1])/m

figure(1)
plot(t,v,'r-*')
hold('on')
plot(goal*ones(shape(t)),v,'b-')
title('velocity as a function of time. vt = %.2f m/s' % vt)
legend('v(t)','Goal')
xlabel('time [s]')
ylabel('[m/s]')

figure(2)
plot(t,x,'r-*')
hold('on')
plot(t,100*ones(shape(x)),'b-')
plot(goal*ones(shape(t)),x,'b-')
title('position as a function of time. finish time = %.2f s' % goal)
xlabel('time [s]')
ylabel('[m]')
hold('off')

figure(3)
plot(t,Drag,'b-d')
hold('on')
plot(t,Fd,'r-*')
plot(t,Fc,'g-o')
plot(t,F_v,'bk-x')
plot(t,F*ones(shape(Fc)),'m-')
title('Drag and driving forces as function of time')
legend('Drag(t)','Fd(t)','Fc(t)','Fv(t)','F(t)')
xlabel('time [s]')
ylabel('[N]')

figure(4)
plot(t,Acc,'r-*')
title('Acceleration as function of time. dt = %.3f s' % dt)
legend('a(t)')
xlabel('time [s]')
ylabel('[m/s^2]')





    