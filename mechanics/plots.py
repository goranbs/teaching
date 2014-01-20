# Ukesoppgaver, plotte posisjon,hastighet og aksellerasjon

from scitools.std import *

alpha = 3.0             # m/s**3
t0 = 0
tn = 10
x0 = 0
t = linspace(t0,tn,100) # timearray

def x(t):
    return x0 + alpha*(t**3 - t0**3)/6

def v(t):
    return (alpha*t**2)/2

def a(t):
    return alpha*t

figure()
plot(t, x(t),'r-o')
hold('on')
plot(t, v(t),'b-*')
plot(t, a(t),'y-d')
legend('x(t)','v(t)','a(t)')
xlabel('time [s]')
hold('off')
