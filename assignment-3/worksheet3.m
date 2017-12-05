%% explicit definition of function. ie., population curve
clc;
close all;
clear all;
dt = [1/2,1/4,1/8,1/16,1/32];
tend = 5;
y0 = 20;    
t = 0:0.001:tend;                  % small time step for exact function
p = 200./(20-10*exp(-7*t));             % exact value

%% %% Plotting and computing Euler method results
Euler_E = zeros(1,size(dt,2));
figure(1)
plot(t,p,'DisplayName','Exact function')
title('Euler method results')
axis([0 tend 0 15]);
hold on;
best_euler = euler(y0,dt(end),tend);% To compute approximate error
for i = 1:size(dt,2)
    yeuler = euler(y0,dt(i),tend);
    figure(1)
    plot(0:dt(i):tend,yeuler,'DisplayName',strcat('dt = ',string(dt(i))));
    hold on;
    pexact = 200./(20-10*exp(-7*(0:dt(i):tend)));
    Euler_E(i) = sqrt((dt(i)/tend)*sum((yeuler - pexact).^2));
end
legend('show')
error_reduction = zeros(1,size(dt,2));
for i = 2:size(dt,2)
    error_reduction(i) = Euler_E(i-1)/Euler_E(i);
end
Error_mat_euler = [dt;Euler_E;error_reduction];
error_table_euler = array2table(Error_mat_euler,'RowNames',{'dt','error','error red'})

%% %% Plotting and computing Heun method results
Heun_E = zeros(1,size(dt,2));
figure(2)
plot(t,p,'DisplayName','Exact function')
title('Heun method results')
axis([0 tend 0 15]);
hold on;
best_heun = heun(y0,dt(end),tend);% To compute approximate error
for i = 1:size(dt,2)
    yheun = heun(y0,dt(i),tend);
    figure(2)
    plot(0:dt(i):tend,yheun,'DisplayName',strcat('dt = ',string(dt(i))));
    hold on;
    pexact = 200./(20-10*exp(-7*(0:dt(i):tend)));
    Heun_E(i) = sqrt((dt(i)/tend)*sum((yheun - pexact).^2));
end
legend('show')
error_reduction = zeros(1,size(dt,2));
for i = 2:size(dt,2)
    error_reduction(i) = Heun_E(i-1)/Heun_E(i);
end
Error_mat_heun = [dt;Heun_E;error_reduction];
error_table_heun = array2table(Error_mat_heun,'RowNames',{'dt','error','error red'})

%% %% Computing Implicit Euler method results
imp_euler_E = zeros(1,size(dt,2));
figure(3)
plot(t,p,'DisplayName','Exact function')
title('Implicit Euler method results')
axis([0 tend 0 15]);
hold on;
for i = 1:size(dt,2)
    yimpeuler = imp_euler(y0,dt(i),tend);
    figure(3)
    plot(0:dt(i):tend,yimpeuler,'DisplayName',strcat('dt = ',string(dt(i))));
    hold on;
    pexact = 200./(20-10*exp(-7*(0:dt(i):tend)));
    imp_euler_E(i) = sqrt((dt(i)/tend)*sum((yimpeuler - pexact).^2));
end
legend('show')
error_reduction = zeros(1,size(dt,2));
for i = 2:size(dt,2)
    error_reduction(i) = imp_euler_E(i-1)/imp_euler_E(i);
end
Error_mat_imp_euler = [dt;imp_euler_E;error_reduction];
error_table_imp_euler = array2table(Error_mat_imp_euler,'RowNames',{'dt','error','error red'})


%% %% Computing Adams Moulton's method results
ad_mou_E = zeros(1,size(dt,2));
figure(4)
plot(t,p,'DisplayName','Exact function')
title('Adams Moulton method results')
axis([0 tend 0 15]);
hold on;
for i = 1:size(dt,2)
    yam = ad_mou(y0,dt(i),tend);
    figure(4)
    plot(0:dt(i):tend,yam,'DisplayName',strcat('dt = ',string(dt(i))));
    hold on;
    pexact = 200./(20-10*exp(-7*(0:dt(i):tend)));
    ad_mou_E(i) = sqrt((dt(i)/tend)*sum((yam - pexact).^2));
end
legend('show')
error_reduction = zeros(1,size(dt,2));
for i = 2:size(dt,2)
    error_reduction(i) = ad_mou_E(i-1)/ad_mou_E(i);
end
Error_mat_ad_mou = [dt;ad_mou_E;error_reduction];
error_table_ad_mou = array2table(Error_mat_ad_mou,'RowNames',{'dt','error','error red'})


%% %% Computing Adams Moulton's method results (Linearised 1)
ad_mou_L1E = zeros(1,size(dt,2));
figure(5)
plot(t,p,'DisplayName','Exact function')
title('Adams Moulton (Linearised 1) method results')
axis([0 tend 0 15]);
hold on;
for i = 1:size(dt,2)
    yam = ad_mou_L1(y0,dt(i),tend);
    figure(5)
    plot(0:dt(i):tend,yam,'DisplayName',strcat('dt = ',string(dt(i))));
    hold on;
    pexact = 200./(20-10*exp(-7*(0:dt(i):tend)));
    ad_mou_L1E(i) = sqrt((dt(i)/tend)*sum((yam - pexact).^2));
end
legend('show')
error_reduction = zeros(1,size(dt,2));
for i = 2:size(dt,2)
    error_reduction(i) = ad_mou_L1E(i-1)/ad_mou_L1E(i);
end
Error_mat_ad_mouL1 = [dt;ad_mou_L1E;error_reduction];
error_table_ad_mouL1 = array2table(Error_mat_ad_mouL1,'RowNames',{'dt','error','error red'})


%% %% Computing Adams Moulton's method results (Linearised 2)
ad_mou_L2E = zeros(1,size(dt,2));
figure(6)
plot(t,p,'DisplayName','Exact function')
title('Adams Moulton (Linearised 2) method results')
axis([0 tend -50 50]);
hold on;
for i = 1:size(dt,2)
    yam = ad_mou_L2(y0,dt(i),tend);
    figure(6)
    plot(0:dt(i):tend,yam,'DisplayName',strcat('dt = ',string(dt(i))));
    hold on;
    pexact = 200./(20-10*exp(-7*(0:dt(i):tend)));
    ad_mou_L2E(i) = sqrt((dt(i)/tend)*sum((yam - pexact).^2));
end
legend('show')
error_reduction = zeros(1,size(dt,2));
for i = 2:size(dt,2)
    error_reduction(i) = ad_mou_L2E(i-1)/ad_mou_L2E(i);
end
Error_mat_ad_mouL2 = [dt;ad_mou_L2E;error_reduction];
error_table_ad_mouL2 = array2table(Error_mat_ad_mouL2,'RowNames',{'dt','error','error red'})

%% %% Function to compute using implicit Euler method
function y = imp_euler(y0,dt,tend)
y = zeros(1,tend/dt+1);
y(1) = y0;
for i = 2:(tend/dt+1)
    y(i) = newton(10^-4,y(i-1),dt,@G,@Gp);
end
end

%% %% Function to compute using Adams Moulton method
function y = ad_mou(y0,dt,tend)
y = zeros(1,tend/dt+1);
y(1) = y0;
for i = 2:(tend/dt+1)
    y(i) = newton(10^-4,y(i-1),dt,@F,@Fp);
end
end

%% %% Function definition for computing the roots using Newton method
function root = newton(accuracy,yn,dt,f,fp)
x = yn + 25*rand;                 % Just a guess. Shouldnt matter much where yn is.
error = 10;                                % Just to go inside the loop
iterations = 0;
if fp(x,dt) < 10*eps
    x = x + 5;               
end
while (error>accuracy) && iterations < 15
    x = x - f(x,yn,dt)/fp(x,dt);
    error = abs(f(x,yn,dt));
    iterations = iterations + 1;
end
if error > accuracy
    root = NaN;
else
    root = x;
end
end

%% %% Function definition for F prime
function ddifferential = Fp(x,dt)
p = [0.7*dt 1-3.5*dt];
ddifferential = polyval(p,x);
end

%% %% Function definition for F
function value = F(x,k,dt)
value = x - k - (dt / 2) * (diff(k) + diff(x));
end

%% %% Function definition for computing p'' (or G prime) (needed for Newton method (Implicit Euler))
function ddifferential = Gp(x,dt)  
p = [1.4*dt 1-7*dt];
ddifferential = polyval(p,x);
end

%% %% Function to compute G(x) (Needed for Newton method (Implicit Euler))
function value = G(x,k,dt)
value = x - dt*diff(x) - k;
end

%% %% Function to compute linearised Adams Moulton method (1)
function y = ad_mou_L1(y0,dt,tend)
y = zeros(1,tend/dt+1);
y(1) = y0;
for i = 2:(tend/dt+1)
    y(i) = (y(i-1) + (dt/2)*(diff(y(i-1)) + 7*y(i-1)))/(1 + 0.35*dt*y(i-1));
end
end

%% %% Function to compute linearised Adams Moulton method (2)
function y = ad_mou_L2(y0,dt,tend)
y = zeros(1,tend/dt+1);
y(1) = y0;
for i = 2:(tend/dt+1)
    y(i) = (y(i-1) + (dt/2)*(diff(y(i-1))))/(1 - dt*3.5*(1 - 0.1*y(i-1)));
end
end


%% %%Function definition for computing p'

function differential = diff(x)
differential = 7*(1-x/10)*x;
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
