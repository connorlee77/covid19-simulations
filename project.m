%% Constants
zeta = 0.2;
sigma = 1/5.1;
gamma = 0.154;
alpha = 0.034;
rho = 1/17.8;
beta = 2.2/6.5;
N = 1e5;

%% Intial conditions
I0 = 10; % 10 infected
S0 = N - I0;
E0 = 0;
R0 = 0;
D0 = 0;

x0 = [S0, E0, I0, R0, D0]';

%% Time length
TOTAL_TIME = 900; % days
dt = 0.1;
TOTAL_STEPS = length(0:dt:TOTAL_TIME);

%% Simulate
x = zeros(5,TOTAL_STEPS);
x(:,1) = x0;
for i=2:TOTAL_STEPS
    
    %%% feedback linear beta
    % beta = feedback(...)
    %%%
    
    x(:,i) = x(:,i-1) + dynamics(x(:,i-1), zeta, sigma, gamma, alpha, rho, N, beta)*dt;
end

subplot(3, 2, 1)
plot(x(1,:), 'displayname', 'S')
title('S')

subplot(3, 2, 2)
plot(x(2,:), 'displayname', 'E')
title('E')

subplot(3, 2, 3)
plot(x(3,:), 'displayname', 'I')
title('I')

subplot(3, 2, 4)
plot(x(4,:), 'displayname', 'R')
title('R')

subplot(3, 2, 5)
plot(x(5,:), 'displayname', 'D')
title('D')