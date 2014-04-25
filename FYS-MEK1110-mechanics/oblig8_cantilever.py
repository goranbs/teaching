'''
Mathematical mandatory exercise, Obli8 - FYS-MEK1110

                  Micro-Electromechanical system

Keywords: Cantilver, rotation, moment of inertia.
'''

from math import pi,sin
from scitools.std import *
import matplotlib.pyplot as plt

M = 0.0001
g = 9.81
kappa = 0.1
L = 0.0001
I = 22.0/3*M*L**2
theta_max = M*g*L*5/(2*kappa)
theta = linspace(0,theta_max,10)
Ug = -5.0/4*L*M*g*sin(theta)
Uh = 0.5*kappa*(theta[:]*theta[:])

# i) Describe the motion of the cantilever

omega = (Ug[:]+Uh[:])*(2/I)

plt.figure()
plt.plot(theta,omega)
plt.title('omega squared')
plt.legend(['omega(theta)'])
plt.xlabel('theta')
plt.ylabel('omega^2')


# k) potential energy of the cantilever

plt.figure()
plt.plot(theta,Uh,'r-')
plt.hold(True)
plt.plot(theta,Ug,'b-')
plt.plot(theta,Uh+Ug,'g--')
plt.plot(theta,zeros(len(theta)),'y-')
plt.title('Potential energy')
plt.legend(['Uh','Ug','Uh+Ug','zero'])
plt.xlabel('theta')
plt.ylabel('potential energy')
plt.hold(False)


plt.show(True)


