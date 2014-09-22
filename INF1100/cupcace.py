# -*- coding: utf-8 -*-
"""
Created on Sun Sep 21 19:05:28 2014

@author: Yang Trang

"""

import matplotlib.pyplot as mplt
import math
import numpy
'''
The program is for plotting data from fall-experiments
with cupcake cups. The code is rewritten from
Arnt Inge Vistnes’ code in Matlab suitable for execution
in Python.
'''
# (1)
def avg(numbers):
    n = float(len(numbers))
    return sum(numbers)/n
z1 = [0.0, 5.11, 9.48, 13.86]
z2 = [0,0, 5.11, 9.45, 13.87]
# east wing
# west wing
# (2)
# fall time from 1.floor
t11 = avg([3.62, 3.83, 3.41, 3.60])
t12 = avg([2.19, 1.98, 2.40, 2.17])
t14 = avg([1.88, 2.02, 1.74, 1.93])
# fall time from 2.floor
t21 = avg([6.63, 6.69, 6.48, 6.85])
t22 = avg([4.31, 3.97, 4.47, 5.00])
t24 = avg([3.43, 2.84, 3.53, 3.78])
# fall time from 3.floor
t31 = avg([10.0, 9.21, 10.80, 10.18])
t32 = avg([5.79, 5.55, 6.02, 5.84])
t34 = avg([4.53, 4.96, 5.00, 4.81])
# (3)
T1 = [0, t11, t21, t31]
T2 = [0, t12, t22, t32]
T4 = [0, t14, t24, t34]
# (5)
mplt.figure()
mplt.plot(z1, T1, '.-r')
mplt.plot(z1, T2, '.-b')
mplt.plot(z1, T4, '.-g')
mplt.xlabel('Drop [m]')
mplt.ylabel('Fall time [s]')
mplt.legend(['Rel. mass 1', 'Rel.mass 2', 'Rel. mass 4'])
mplt.title('Fall time as function of fall height for three relative masses')
10# comments about the result of the plot
print '\nThe plot shows three fairly straight lines which'
print 'indicates that the velocity is constaint for all heights.'
print 'The cupcake cups acheive terminal velocity pretty fast after'
print 'initialization. The terminal velocity differ with the mass.\n'
# (6)
v1 = avg([z1[i]/T1[i] for i in range(1,4)])
v2 = avg([z1[i]/T2[i] for i in range(1,4)])
v4 = avg([z1[i]/T4[i] for i in range(1,4)])
print 'Average speed for relativ mass 1 is %.2f m/s' % v1
print 'Average speed for relativ mass 2 is %.2f m/s' % v2
print 'Average speed for relativ mass 4 is %.2f m/s\n' % v4
print 'The results show that the terminal velocity is proportional'
print 'with the square root of relative mass. And they are NOT'
print 'valid with Galilei’s observations in 1589.\n'
# quick indication if the terminal velocity is proportional with sqrt(m)
mplt.figure()
mplt.plot(z1, [t/math.sqrt(2) for t in T1], '.-k')
mplt.plot(z1, [t/2 for t in T1], '.-c')
mplt.xlabel('time [s]')
mplt.ylabel('height [m]')
mplt.title('terminal velocity prop. sqrt(m)')
mplt.legend(['T1/sqrt(2)', 'T1/2'])
mplt.show()


v_term_cupcake = 10

rho = 1.2041
r = (6.5/2)*0.01
g = 9.81
m = 0.32*0.001
C = (2.0*m*g)/(v_term_cupcake**2*rho**math.pi*r**2)
print 'Factor of friction C is %.2f ' % C

N = 100
g = 9.81
h = numpy.linspace(0, 14, N)
t = numpy.sqrt(2*h/float(g))
mplt.figure()
mplt.plot(h, t, '-b')
mplt.xlabel('Drop height [m]')
mplt.ylabel('Fall time [s]')
mplt.title('Fall time vs fall height no air resistance')
mplt.show()