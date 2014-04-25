'''
Mathematical mandatory exercise, Obli8 - FYS-MEK1110

                  Micro-Electromechanical system

Keywords: Cantilver, rotation, moment of inertia.
'''

from math import pi,sin
from scitools.std import *
import numpy as np
import matplotlib.pyplot as plt

n = 10
LMg = 0.1
g = 9.81
k = [0.05, 0.1, 0.2]
theta = np.zeros((len(k),n))
Ug = np.zeros((len(k),n))
Uh = np.zeros((len(k),n))
omega = np.zeros((len(k),n))

for i in range(len(k)):
    kappa = k[i]
    I = 22.0/3*(LMg/g)**2
    theta_max = LMg*5/(2*kappa)
    theta[i,:] = linspace(0,theta_max,10)
    Ug[i,:] = -5.0*LMg*sin(theta[i,:])
    Uh[i,:] = 0.5*kappa*(theta[i,:]*theta[i,:])
    omega[i,:] = (Ug[i,:]+Uh[i,:])*(2/I)
# i) Describe the motion of the cantilever



plt.figure()
plt.plot(theta[0,:],omega[0,:],'r-')
plt.hold(True)
plt.plot(theta[1,:],omega[1,:],'b-')
plt.plot(theta[2,:],omega[1,:],'g-')
plt.title('omega squared')
plt.legend(['omega(theta)'])
plt.xlabel('theta')
plt.ylabel('omega^2')
plt.hold(False)

# k) potential energy of the cantilever


plt.figure()
plt.plot(theta[2,:],Uh[2,:],'r-')
plt.hold(True)
plt.plot(theta[2,:],Ug[2,:],'b-')
plt.plot(theta[2,:],Uh[2,:]+Ug[2,:],'g--')
plt.plot(theta[2,:],zeros(len(theta[2,:])),'y-')
plt.title('Potential energy')
plt.legend(['Uh','Ug','Uh+Ug','zero'])
plt.xlabel('theta')
plt.ylabel('potential energy')
plt.hold(False)

plt.show(True)


