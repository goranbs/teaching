########################################################
# Oblig 5 - Midtveis/Hjemmeeksamen - FYS-MEK1110 - v2014
#
# exercise e) plot horizontal component of spring force

from scitools.std import *

y = linspace(-2,2,100)
d= 0.5
k = 100
m = 1.0
g = 9.81

Fx = k*d*(1-d/(sqrt(y*y + d*d)))
Fy = -k*y*(1-d/(sqrt(y*y + d*d))) 
import matplotlib.pyplot as plt

# plotting
plt.figure()                                         #figure #1
plt.subplot(2,1,1)
plt.plot(y,Fx)
plt.title('d) Horizontal force compontent, springforce')
plt.legend(('Fx(y)'))
plt.ylabel('[N]')
plt.xlabel('[m]')

plt.subplot(2,1,2)              
plt.plot(y,Fy)
plt.title('e) Vertical force compontent, springforce')
plt.legend(('Fy(y)'))
plt.ylabel('[N]')
plt.xlabel('[m]')

#plt.show(True)
plt.show(False)
print 'd) see plot'
print 'e) see plot'

#############################################################################
# Find the equilibrium position of the cylinder, using 
# fsolve from scipy.optimize                            exercise f)

from scipy.optimize import fsolve

def sum_Forces_y(y):
    return -k*y*(1-d/sqrt(d*d + y*y)) - m*g

roots = fsolve(sum_Forces_y,0.0)
#print roots
'''
goranbs@kracken:~/goran/teaching/FYS-MEK1110-mechanics$ python ob5_plot_hc.py 
[ 0.41944805]

'''

# we may test the value returned by fsolve by inserting it into 
# the function sum_Forces_y:

is_zero = sum_Forces_y(roots[0])
print 'f) equilibrium position y=%.3f m, f(y)=%.8f' % (roots[0],is_zero) # -6.44469443145e-07 -which is very close to zero!! :-)

##############################################################################
# Finding the verical position of the cylinder as a function of time
#                                                          exercise g)

y0 = -1.0  # [m] - initial vertical position of cylinder

dt = 0.001 # timestep
t0 = 0.0   # [s]
tmax = 2.0 # [s]
n = int((tmax-t0)/dt) # number of timesteps

ay = zeros(n)
vy = zeros(n)
y_pos = zeros(n)
F_y = zeros(n)
G = zeros(n)
t = zeros(n)

y_pos[0] = y0

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

print 'g) highest position of cylinder=%.4f m' % max(y_pos)

############################################################################
# Plotting
plt.figure()                                     # figure #2
plt.subplot(3,1,1)
plt.plot(t,ay)
plt.title('g) Acceleration, Velocity and position')
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
plt.legend(('y(t)'))

print 'g) see plots #2 and #3'
#plt.show(False)

plt.figure()                                      # figure #3
plt.plot(t,F_y,'g-')
plt.hold(True)
plt.plot(t,G,'r-')
plt.plot(t,F_y[:] + G[:],'g-*')
plt.title('g) Forces in vertical direction')
plt.xlabel('time [s]')
plt.ylabel('Force [N]')
plt.legend(('Fy(t)','G(t)','F_tot(t)'))
plt.hold(False)

#plt.show(False)

############################################################################
# Plot the kinetic energy as a function of time. exercise h)

ek = 0.5*m*vy*vy

plt.figure()                                       # figure #4
plt.plot(t,ek,'r-*')
plt.title('h) Kinetic energy as function f time')
plt.xlabel('time [s]')
plt.ylabel('Kinetic energy [Nm]')
plt.legend(('Ek(t)'))

plt.show(False)

#############################################################################
# Find the work needed to bring the cylinder from its equilibrium position to
# y = -1.
# Work = F*d, where F is the "working force" and d is the displacement of the 
# object that the force is working on.
#                                                  exercise i)



def integrate(f,x0,x,dx):
    #noe er galt med denne integratoren.. finn ut av det, eller bruk quad.
    n = int((x - x0)/dx) # integration steps
    n = abs(n)
    print 'in integral; n=%f' % n
    integral = 0
    for i in range(n-1):
        a = x0*(1+i*dx)
        b = x0*(2+i*dx)
        integral += (f(a) + f(b))

    integral = 0.5*dx*integral
    return integral

dy = 0.001
work = integrate(sum_Forces_y,-0.415,-1,dy)
print 'i) Work=%.7f Nm integral from y= %.4f to y=-1, with dy=%.4f' % (-work,roots,dy)
from scipy.integrate import quad
work2 = quad(sum_Forces_y,-0.414,-1)
print work2
#############################################################################
# Modify the program to include the effect of friction.
# my_d = 0.1                                       exercise j)

my_d = 0.1   # dynamic friction

def friction(N_force,vel):
    # the friction is proportioal to the Normal force N.
    # N is dependent of the horizontal spring force component.
    # and is directed in the opposit direction of the velocity
    if vel == 0:
        return 0
    else:
        return -my_d*abs(N_force)*(vel/abs(vel))

a = zeros(n)
v = zeros(n)
y = zeros(n)
F_y = zeros(n)
F_x = zeros(n)
G = zeros(n)
N = zeros(n)
R = zeros(n)
t = zeros(n)

y[0] = y0

for i in range(n-1):

    F_y[i] = -k*y_pos[i]*(1 - d/sqrt(d*d + y[i]*y[i]))
    F_x[i] = k*d*(1 -d/(sqrt(d*d + y[i]*y[i])))
    N[i] = -F_x[i]
    R[i] = friction(N[i],v[i])
    G[i] = -m*g
   
    #print y_[i]    
    #print 'Normal force=%.4f friction= %.4f' % (N[i],R[i])

    a[i] = (F_y[i] + G[i] + R[i])/m
    v[i+1] = v[i] + a[i]*dt
    y[i+1] = y[i] + v[i+1]*dt
    t[i+1] = t[i] + dt

a[-1] = a[-2]
F_y[-1] = F_y[-2]
G[-1] = G[-2]
R[-1] = R[-2]
N[-1] = N[-2]
F_x[-1] = F_x[-2]

print 'k) highest position of cylinder=%.4f m' % max(y)

# Plotting
plt.figure()                                     # figure #5
plt.subplot(3,1,1)
plt.plot(t,a)
plt.title('j) Acceleration, Velocity and position')
plt.xlabel('time [s]')
plt.ylabel('acceleration [m/s^2]')
plt.legend(('a(t)'))
plt.subplot(3,1,2)
plt.plot(t,v)
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.legend(('v(t)'))
plt.subplot(3,1,3)
plt.plot(t,y)
plt.xlabel('time [s]')
plt.ylabel('position [m]')
plt.legend(('z(t)'))


plt.figure()                                     # figure #6
plt.plot(t,F_y,'b-')
plt.hold(True)
plt.plot(t,G,'r-')
plt.plot(t,R,'y-')
plt.plot(t,F_y[:] + G[:] + R[:],'b-*')
plt.title('j) Forces in vertical direction')
plt.xlabel('time [s]')
plt.ylabel('Force [N]')
plt.legend(('Fy(t)','G(t)','R(t)','F_tot(t)'))
plt.hold(False)

print 'j) see plots #5 and #5'
plt.show(False)
############################################################################
# Plot the kinetic energy as a function of time. exercise l)

ek = 0.5*m*v*v

plt.figure()                                     # figure #7
plt.plot(y,ek,'r-*')
plt.title('l) Kinetic energy as function of height')
#plt.title('l) Kinetic energy as function f time')
plt.xlabel('y [m]')
#plt.xlabel('time [s]')
plt.ylabel('Kinetic energy [Nm]')
plt.legend(('Ek(y)'))
#plt.legend(('Ek(t)'))
print 'l) see plot #7'
############################################################################
# two springs 

import numpy as np

k1 = 100
k2 = 99
m = 1.0

r = np.zeros((n,2))
a = np.zeros((n,2))
v = np.zeros((n,2))
Fk1 = np.zeros((n,2))
Fk2 = np.zeros((n,2))
G = np.zeros((n,2))
t = np.zeros(n)

lo1 = np.array([0, 0.5])
lo2 = np.array([0, -0.5])

#print lo1
#print lo2

r[0,0] = -0.001     # initial x-position
r[0,1] = 0.0        # initial y-position
v[0,0] = 0.0        # initial velocity x
v[0,1] = 0.0        # initial velocity y

for i in range(n-1):
    Fk1[i,:] = -k1*(linalg.norm(r[i,:]+lo1[:]) - linalg.norm(lo1[:]))*(r[i,:]+lo1[:])/linalg.norm(r[i,:]+lo1)
    Fk2[i,:] = -k2*(linalg.norm(r[i,:]+lo2[:]) - linalg.norm(lo2[:]))*(r[i,:]+lo2[:])/linalg.norm(r[i,:]+lo2)
    G[i,1] = -m*g
    #print 'FK1_y=%.4f FK2_y=%.4f' % (Fk1[i,1],Fk2[i,1])
    
    a[i,:] = (Fk1[i,:] + Fk2[i,:] + G[i,:])/m
    
    v[i+1,:] = v[i,:] + a[i,:]*dt
    r[i+1,:] = r[i,:] + v[i+1,:]*dt

    t[i+1] = t[i] + dt

Fk1[-1,:] = Fk1[-2,:]
Fk2[-1,:] = Fk2[-2,:]
G[-1,:] = G[-2,:]
a[-1,:] = a[-2,:]

# plotting
# Plotting
plt.figure()                                     # figure #8
plt.subplot(3,1,1)
plt.plot(t,a[:,1])
plt.title('y-direction, r0=[%.3f,%.3f], v0=[%.3f,%.3f]'% (r[0,0],r[0,1],v[0,0],v[0,1]))
plt.xlabel('time [s]')
plt.ylabel('acceleration [m/s^2]')
plt.legend(('ay(t)'))
plt.subplot(3,1,2)
plt.plot(t,v[:,1])
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.legend(('vy(t)'))
plt.subplot(3,1,3)
plt.plot(t,r[:,1])
plt.xlabel('time [s]')
plt.ylabel('position [m]')
plt.legend(('y(t)'))

plt.figure()                                     # figure #9
plt.subplot(3,1,1)
plt.plot(t,a[:,0])
plt.title('x-direction, r0=[%.3f,%.3f], v0=[%.3f,%.3f]'% (r[0,0],r[0,1],v[0,0],v[0,1]))
plt.xlabel('time [s]')
plt.ylabel('acceleration [m/s^2]')
plt.legend(('ax(t)'))
plt.subplot(3,1,2)
plt.plot(t,v[:,0])
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.legend(('vx(t)'))
plt.subplot(3,1,3)
plt.plot(t,r[:,0])
plt.xlabel('time [s]')
plt.ylabel('position [m]')
plt.legend(('x(t)'))


plt.figure()                                     # figure #10
plt.plot(t,Fk1[:,1],'g-')
plt.hold(True)
plt.plot(t,Fk2[:,1],'r-')
plt.plot(t,G[:,1],'m-')
plt.plot(t,Fk1[:,1] + G[:,1] + Fk2[:,1],'b-*')
plt.title('two springs, forces in vertical direction')
plt.xlabel('time [s]')
plt.ylabel('Force [N]')
plt.legend(('Fk1_y(t)','Fk2_y(t)','G(t)','F_y(t)'))
plt.hold(False)


plt.figure()                                     # figure #11
plt.plot(r[:,0], r[:,1], 'b-')
plt.title('Two springs. Position r')
plt.xlabel('x [m]')
plt.ylabel('y [m]')
plt.legend('r')

print '----------------------------------------------------'
#plt.show(True)
plt.show(False)
