function NewtonsCradle()
% Newtons Cradle - Mandatory assignment nr.7 FYS-MEK1110

N = 3;
m = 1.0;
k =1000.0;
q = 1.0;
d = 0.1;
v0 = 1.0;
time = 1.0;
dt = 0.001;

n = ceil(time/dt);
x = zeros(n,N);
v = zeros(n,N);
F = zeros(n,N);
t = zeros(n,1);

for j=1:N
    x(1,j) = d*(j-1);
end
v(1,1)=v0;
for i=1:n-1
    % Find the forces oon each block,j.
    % First, force from block to the left.
    for j=2:N
        dx = x(i,j)-x(i,j-1);
        F(i,j)=F(i,j)+force(dx,d,k,q);
    end
    % Second, force from block to the right
    for j=1:N-1
        dx = x(i,j+1)-x(i,j);
        F(i,j) = F(i,j) - force(dx,d,k,q);
    end
    % Euler-Chromer:
    
    a = F(i,:)/m;
    v(i+1,:) = v(i,:) + a*dt;
    x(i+1,:) = x(i,:) + v(i+1,:)*dt;
    t(i+1,1) = t(i,1) + dt;
end



    function F = force(dx,d,k,q)
        if dx<d
            F = k*abs(dx-d).^q;
        else
            F = 0.0;
        end
    end

% Plotting the results:

fontsize=18;
legends = ['ball 1', 'ball 2', 'ball 3'];
f = figure()
for j=1:N
    plot(t,v(:,j));
    tmp = ['Ball ' num2str(j,'%g')];
    if j==1
        hold on;
    end
    if j==N
        hold off
    end
end

set(gca,'FontSize',fontsize)
%get(0,'DefaultAxesColorOrder')
name = ['Newtons Cradle. k=' num2str(k,'%.1f') ' q=' num2str(q,'%.1f')];
title(name)
xlabel('t [s]')
ylabel('v,[m/s]')
legend('ball 1','ball 2', 'ball 3')
Title = ['plot_NewtonsCradle_k' num2str(k,'%g') '_q' num2str(q,'%g') '.png'];
print(f,'-dpng', Title) % save figure as png with name 'name'

end