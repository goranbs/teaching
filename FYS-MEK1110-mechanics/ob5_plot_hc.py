########################################################
# Oblig 5 - Midtveis/Hjemmeeksamen - FYS-MEK1110 - v2014
#
# exercise e) plot horizontal component of spring force

from scitools.std import linspace, sqrt, zeros

y = linspace(-2,2,100)
d= 0.5
k = 100
m = 1.0
g = 9.81

Fx = k*d*(1-d/(sqrt(y*y + d*d)))
Fy = -k*y*(1-d/(sqrt(y*y + d*d))) 
import matplotlib.pyplot as plt

plt.figure()
plt.subplot(2,1,1)
plt.plot(y,Fx)
plt.title('Horizontal force compontent')
plt.legend(('Fx(y)'))
plt.xlabel('[N]')
plt.ylabel('[m]')

plt.subplot(2,1,2)
plt.plot(y,Fy)
plt.title('Vertical force compontent')
plt.legend(('Fy(y)'))
plt.xlabel('[N]')
plt.ylabel('[m]')

#plt.show(True)
plt.show(False)

#############################################################################
# Find the equilibrium position of the cylinder, using 
# fsolve from scipy.optimize

from scipy.optimize import fsolve

def sum_Forces_y(y):
    return k*y*(1-d/sqrt(d*d + y*y)) - m*g

roots = fsolve(sum_Forces_y,0.0)
print roots
'''
goranbs@kracken:~/goran/teaching/FYS-MEK1110-mechanics$ python ob5_plot_hc.py 
[ 0.41944805]

'''

# we may test the value returned by fsolve by inserting it into 
# the function sum_Forces_y:

is_zero = sum_Forces_y(0.41944804)
print is_zero # -6.44469443145e-07 -which is very close to zero!! :-)

##############################################################################
# Finding the verical position of the cylinder as a function of time

y0 = -1.0  # [m] - initial vertical position of cylinder

dt = 0.001 # timestep
t0 = 0.0   # [s]
tmax = 10.0 # [s]
n = int((tmax-t0)/dt) # number of timesteps

ay = zeros(n)
vy = zeros(n)
y_pos = zeros(n)
F_y = zeros(n)
G = zeros(n)
t = zeros(n)

for i in range(n-1):
    F_y[i] = -k*y_pos[i]*(1 - d/sqrt(d*d + y_pos[i]*y_pos[i]))
    G[i] = -m*g
    ay[i] = (F_y[i] + G[i])/m
    vy[i+1] = vy[i] + ay[i]*dt
    y_pos[i+1] = y_pos[i] + vy[i+1]*dt
    t[i+1] = t[i] + dt

ay[-1] = ay[-2]
F_y[-1] = F_y[-2]
G[-1] = G[-2]


############################################################################
# Plotting
plt.figure()                                     # figure #1
plt.subplot(3,1,1)
plt.plot(t,ay)
plt.title('Acceleration, Velocity and position')
plt.xlabel('time [s]')
plt.ylabel('acceleration [m/s^2]')
plt.legend(('a(t)'))
plt.subplot(3,1,2)
plt.plot(t,vy)
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.legend(('v(t)'))
plt.subplot(3,1,3)
plt.plot(t,y_pos)
plt.xlabel('time [s]')
plt.ylabel('position [m]')
plt.legend(('z(t)'))

#plt.show(True)

plt.figure()                     # forceplot
plt.plot(t,F_y)
plt.hold(True)
plt.plot(t,G)
plt.title('Forces in vertical direction')
plt.xlabel('time [s]')
plt.ylabel('Force [N]')
plt.legend(('Fy(t)','G(t)'))

##############################################################
# Plot the kinetic energy as a function of time

ek = 0.5*m*vy*vy

plt.figure()
plt.plot(t,ek,'r-*')
plt.title('Kinetic energy as function f time')
plt.xlabel('time [s]')
plt.ylabel('Kinetic energy [Nm]')
plt.legend(('Ek(t)'))

plt.show(True)
                       
