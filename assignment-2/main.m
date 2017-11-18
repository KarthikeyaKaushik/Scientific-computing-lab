%% explicit definition of function. ie., population curve
clc;
close all;
clear all;
dt = [1,1/2,1/4,1/8];
tend = 5;
y0 = 1;    
t = 0:0.001:tend;                  % small time step for exact function
p = 10./(1+9*exp(-t));             % exact value

%% %% Plotting and computing Euler method results
Euler_E = zeros(1,size(dt,2));
app_euler = zeros(1,size(dt,2));
figure(1)
plot(t,p,'DisplayName','Exact function')
title('Euler method results')
axis([0 tend 0 15]);
hold on;
best_euler = euler(y0,dt(4),tend);% To compute approximate error
for i = 1:size(dt,2)
    yeuler = euler(y0,dt(i),tend);
    figure(1)
    plot(0:dt(i):tend,yeuler,'DisplayName',strcat('dt = ',string(dt(i))));
    hold on;
    pexact = 10./(1+9*exp(-(0:dt(i):tend)));
    app_euler(i) = sqrt((dt(i)/tend)*sum((yeuler - downsample(best_euler,2^(4-i))).^2));%computing approximate error
    Euler_E(i) = sqrt((dt(i)/tend)*sum((yeuler - pexact).^2));
end
legend('show')


%% %% Plotting and computing Heun method results
Heun_E = zeros(1,size(dt,2));
app_heun = zeros(1,size(dt,2));
figure(2)
plot(t,p,'DisplayName','Exact function')
title('Heun method results')
axis([0 tend 0 15]);
hold on;
best_heun = heun(y0,dt(4),tend);% To compute approximate error
for i = 1:size(dt,2)
    yheun = heun(y0,dt(i),tend);
    figure(2)
    plot(0:dt(i):tend,yheun,'DisplayName',strcat('dt = ',string(dt(i))));
    hold on;
    pexact = 10./(1+9*exp(-(0:dt(i):tend)));
    app_heun(i) = sqrt((dt(i)/tend)*sum((yheun - downsample(best_heun,2^(4-i))).^2));%computing approximate error
    Heun_E(i) = sqrt((dt(i)/tend)*sum((yheun - pexact).^2));
end
legend('show')

%% %% Plotting and computing Runge-Kutta method results
rk_E = zeros(1,size(dt,2));
app_rk = zeros(1,size(dt,2));
figure(3)
plot(t,p,'DisplayName','Exact function')
title('Runge Kutta method results')
axis([0 tend 0 15]);
hold on;
best_rk = rk(y0,dt(4),tend);% To compute approximate error
for i = 1:size(dt,2)
    yrk = rk(y0,dt(i),tend);
    figure(3)
    plot(0:dt(i):tend,yrk,'DisplayName',strcat('dt = ',string(dt(i))));
    hold on;
    pexact = 10./(1+9*exp(-(0:dt(i):tend)));
    app_rk(i) = sqrt((dt(i)/tend)*sum((yrk - downsample(best_rk,2^(4-i))).^2));%computing approximate error
    rk_E(i) = sqrt((dt(i)/tend)*sum((yrk - pexact).^2));
end
legend('show')

%% %%Function definition for computing p'

function differential = diff(x)
differential = (1-x/10)*x;
end


%%  %%Function definition 1: Euler method

function y = euler(y0,dt,tend)
y = zeros(1,tend/dt+1);
y(1) = y0;
for i = 2:(tend/dt+1)
    y(i) = y(i-1) + dt*(diff(y(i-1)));
end
end

%% %% Function definition 2 : Heun method
function y = heun(y0,dt,tend)
y = zeros(1,tend/dt+1);
y(1) = y0;
for i = 2:(tend/dt+1)
    y_temp = y(i-1) + dt*(diff(y(i-1)));
    y(i) = y(i-1) + (dt/2)*(diff(y(i-1)) + diff(y_temp));
end
end

%% %% Function definition 2 : Runge-Kutta method
function y = rk(y0,dt,tend)
y = zeros(1,tend/dt+1);
y(1) = y0;
for i = 2:(tend/dt+1)
    y1 = diff(y(i-1));
    y2 = diff(y(i-1)+(dt/2)*y1);
    y3 = diff(y(i-1)+(dt/2)*y2);
    y4 = diff(y(i-1)+dt*y3);
    y(i) = y(i-1) + (dt/6)*(y1 + 2*y2 + 2*y3 + y4);
end
end