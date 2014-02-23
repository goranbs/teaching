function main
%%
%   Project 7.2 - FYS-MEK1110 - mechanics - oblig #2
%
%                Weather Balloon
%
% Aim is to model the motion of a weather Balloon in two dimentons
% as it is ascending up into the air. 
%%
% constants

radius = 0.3; % radius of balloon
V = 4/3*pi*radius^3;
A = 4*pi*radius^2;
eps = 0.001
vt_x = 0;
vt_z = 0;
m = 0.1;
w0 = 1.0;
h = 10;
g = 9.81;
rho = 1.255;     % 15degrees Celcius and at sealevel.
rho2 = 1.293;    % 20degrees Celcius and at sealevel.
Cd = 1.0;        % shapefactor
D = 0.5*rho2*Cd*A;
Bouyancy = rho2*V*g;
t0 = 0;
tmax = 10;
N = 1000;
dt = (tmax-t0)/N;

A
V
D


a = zeros(N,2);
v = zeros(N,2);
w = zeros(N,2);
r = zeros(N,2);
Fd = zeros(N,2);
G = zeros(N,2);
B = zeros(N,2);
t = zeros(N,1);

%%

    function Drg = Drag(v,w)
        velocity = v - w;
        Drg = D*norm(velocity).*velocity;
    end

    function wnd = wind(z)
        wnd = w0*(1-exp(-z/h));
    end

for i=1:N-1
    w(i,1) = wind(r(i,2));
    G(i,2) = m*g;
    B(i,2) = Bouyancy;
    
    %Fd(i,2) = D*abs(v(i,2))*(v(i,2));
    Fd(i,:) = Drag(v(i,:),w(i,:));
    
    a(i,:) = (B(i,:) - G(i,:) - Fd(i,:))./m;
    v(i+1,:) = v(i,:) + a(i,:).*dt;
    r(i+1,:) = r(i,:) + v(i+1,:).*dt;
    
    if v(i,1) + eps > v(i+1,1)
        vt_x = v(i+1,1);
    end
    
    if v(i+1,2) < v(i,2) + eps 
        vt_z = v(i+1,2);
    end
    t(i+1) = t(i) + dt;

end

a(end,:) = a(end-1,:);

%%
% plotting

figure()

subplot(3,1,1)
plot(t(:),a(:,2))
title('acceleration in vertical direction')
xlabel('time [s]')
ylabel('acceleration [m/s^2]')
legend('a(t)')

subplot(3,1,2)
plot(t(:),v(:,2))
Title1 = ['velocity in the vertical direction, vt=' num2str(vt_z,'%.3f')];
title(Title1)
xlabel('time [s]')
ylabel('velocity [m/s]')
legend('v(t)')

subplot(3,1,3)
plot(t(:),r(:,2))
title('postition int the vertical direction')
xlabel('time [s]')
ylabel('position [m]')
legend('x(t)')

figure()

subplot(3,1,1)
plot(t(:),a(:,1))
title('acceleration in the horizontal direction')
xlabel('time [s]')
ylabel('acceleration [m/s^2]')
legend('a(t)')

subplot(3,1,2)
plot(t(:),v(:,1))
Title2 = ['velocity in the horizontal direction, vt=' num2str(vt_x,'%.3f')];
title(Title2)
xlabel('time [s]')
ylabel('velocity [m/s]')
legend('v(t)')

subplot(3,1,3)
plot(t(:),r(:,1))
title('postition in the horizontal direction')
xlabel('time [s]')
ylabel('position [m]')
legend('x(t)')

figure()
subplot(2,1,1)
plot(t(:),Fd(:,1))
title('Drag force in horizontal direction')
legend('Fd_x(t)')
ylabel('[N]')
xlabel('time [s]')
subplot(2,1,2)
plot(t(:), Fd(:,2))
hold on
plot(t(:),B(:,2))
plot(t(:),G(:,2))
title('Forces in vertical direction')
ylabel('[N]')
xlabel('time [s]')
legend('Fd','B','G')
hold off

vt_x
vt_z
    
end