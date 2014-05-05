from numpy import *
from scitools.easyviz.gnuplot_ import *
time = 10.0
dt = 0.001
m = 23.0
U0 = 150.0
x0 = 2.0
alpha = 39.48
n = int(round(time/dt))
x = zeros(n, float)
v = zeros(n, float)
t = zeros(n, float)
x[0] = -5.0
v[0] = 8.725

for i in range(n-1):
    
    if (abs(x[i])) < x0:                # if atom is in area -x0<=x<=x0
        Fx = -((U0/x0)*(x[i]/abs(x[i]))) # Magnetic force
        Fs = -alpha*v[i]                # Electrostatic force
    else:
        Fs = 0
        Fx = 0
    
    F = Fx + Fs

    a = F/m
    v[i+1] = v[i] + a*dt
    x[i+1] = x[i] + v[i+1]*dt
    t[i+1] = t[i] + dt

plot(t, x, xlabel='t', ylabel='x', legend='Atoms movement', title='Trapping atoms', hardcopy='trapping_atoms2_fysmek_08.png')
