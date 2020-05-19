clc; clear all; close all;
%% Constants
xi = 0.1;
sigma = 1/5.1;
gamma = 0.154;
alpha = 0.034;
rho = 1/17.8;
beta = 2.2/6.5;
N = 8*1e9;
lambda = 0.0000;
% Controller Parameters
alpha_1 = 1;
alpha_2 = 1;
% Scaling parameter
k = 1e-6;
k2 = 1e-5;
k1 = 1e-5;
k3 = 0;


%% Intial conditions
I0 = 0.0015*N; % 10 infected
E0 = 0.005*N;
S0 = N - I0 - E0;
R0 = 0;
D0 = 0;
V0 = 0;
M0 = 1e6;
beta0 = 2.2/6.5;

x0 = [S0, E0, I0, R0, D0, M0, beta0]';

%% Time length
TOTAL_TIME = 100; % days
dt = 0.1;
TOTAL_STEPS = length(0:dt:TOTAL_TIME);

%% Simulate
x = zeros(7,TOTAL_STEPS);
x(:,1) = x0;
u = zeros(1,TOTAL_STEPS);
%maximum_beta = 2.2/6.5;
umax = 0.75;
for i=2:TOTAL_STEPS
    i
    if i > TOTAL_STEPS/2 && i < TOTAL_STEPS*3/4
        k3 = 100000;
    end
    %%% feedback linearization controller u
    u(i) = feedback(x(:,i-1), sigma, gamma, xi, alpha_1, alpha_2, N,k1,k3);
    %%%
    
    u(i) = min(u(i), umax);
    
    x(:,i) = x(:,i-1) + dynamics(x(:,i-1), xi, sigma, gamma,rho, alpha, N,k1,k3,u(i)) * dt;

    %x(7,i) = max(x(7,i),0);
    %x(6,i) = max(x(6,i),0);
    k3 = 0;
end

subplot(2,4, 1)
plot(x(1,:), 'displayname', 'S')
title('S')

subplot(2, 4, 2)
plot(x(2,:), 'displayname', 'E')
title('E')

subplot(2, 4, 3)
plot(x(3,:), 'displayname', 'I')
title('I')

subplot(2, 4, 4)
plot(x(4,:), 'displayname', 'R')
title('R')

subplot(2, 4, 5)
plot(x(5,:), 'displayname', 'D')
title('D')

subplot(2, 4, 6)
plot(x(6,:), 'displayname', 'M')
title('M')

subplot(2, 4, 7)
plot(x(7,:), 'displayname', 'beta')
title('\beta')

subplot(2, 4, 8)
plot(u, 'displayname', 'u')
title('u')