
dt = 0.001;
n = floor(10./dt);
d = 0.5;
k = 100.0;
g = 9.81;
m = 1.0;
a = zeros(n,1);
v = zeros(n,1);
y = zeros(n,1);
t = zeros(n,1);
Fk1 = zeros(n,1);
G = zeros(n,1);
y(1) = -1;

for i=1:n-1
    Fk1(i) = -k*y(i)*(1 - d/(sqrt(y(i)*y(i) + d*d)));
    G(i) = -m*g;
    a(i) = (Fk1(i)+G(i))/m;
    v(i+1) = v(i) + a(i)*dt;
    y(i+1) = y(i) + v(i+1)*dt;
    t(i+1) = t(i) + dt;
end

figure()
plot(t,y)
