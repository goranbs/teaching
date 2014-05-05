
k = 20.0; % spring coeffesient
L0 = 1;  % Length of rope/spring
m = 0.1; % Mass of ball
theta = pi/3; % starting angle
%r0 = [L0*sin(theta) L0*cos(theta)]; % Starting position
v0 = [20 0]; % Starting velocity
r0 = [0 1];

time = 10;
dt = 0.001;
n = round(time/dt);
r = zeros(n,2);
v = zeros(n,2);
t = zeros(n,1);
r(1,:) = r0;
v(1,:) = v0;
g = 9.81;
km = k/m;

for i = 1:n-1
    rr = norm(r(i,:));
    
    if km*(rr - L0) > 0
        a = - g - km*(rr-L0)*r(i,:)/rr;
    else
        a = -g;
    end
        
    v(i+1,:) = v(i,:) + a.*dt;
    r(i+1,:) = r(i,:) + v(i+1,:).*dt;
    t(i+1) = t(i) + dt;
end

plot(r(:,1), r(:,2))
xlabel('x [m]')
ylabel('y [m]')
hold on