%% Constants
xi = 0.0000;
sigma = 1/5.1;
gamma = 0.154;
alpha = 0.034;
rho = 1/17.8;
beta = 2.2/6.5;
N = 8*1e9;

%% CLF parameters
scale = 1;
k1 = 1;
k2 = 1;

p1 = (k1+k1^2+k2^2)/(2*k1*k2);
p2 = 1/(2*k1);
p3 = (k1+1)/(2*k1*k2);
P = [p1, p2; 
    p2, p3];

Q = eye(2);
clf_alpha = 1/max(eig(P));


%% Intial conditions
I0 = 10; % 10 infected
S0 = N - I0;
E0 = 0;
R0 = 0;
D0 = 0;
V0 = 0;

x0 = [S0, E0, I0, R0, D0, V0]';

%% Functions
F = @(S, I, E, R, V, xi, sigma, gamma, alpha, rho, lambda) [
    xi*(R+V) - S*lambda; 
    -sigma*E; 
    sigma*E - gamma*I;
    (1-alpha)*gamma*I - xi*R;
    alpha*rho*I;
    -xi*V + S*lambda];

G = @(S, I, N) [
    -S*I/N; 
    S*I/N; 
    0; 
    0; 
    0;
    0];

eta = @(I, E, sigma, gamma) [I; sigma*E-gamma*I];
V = @(I, E, sigma, gamma) eta(I, E, sigma, gamma)'*P*eta(I, E, sigma, gamma);

dvdeta = @(I, E, sigma, gamma) 2*eta(I, E, sigma, gamma)'*P;
detadx = @(sigma, gamma) [0,0,1,0,0,0;0,sigma,(-1).*gamma,0,0,0];
dvdx = @(I, E, sigma, gamma) dvdeta(I, E, sigma, gamma)*detadx(sigma, gamma);
LfV = @(S, E, I, R, D, V, xi, sigma, gamma, alpha, rho, N, lambda) dvdx(I, E, sigma, gamma)*F(S, I, E, R, V, xi, sigma, gamma, alpha, rho, lambda);
LgV = @(S, E, I, R, D, V, xi, sigma, gamma, alpha, rho, N, lambda) dvdx(I, E, sigma, gamma)*G(S, I, N);

Vdot = @(S, E, I, R, D, V, xi, sigma, gamma, alpha, rho, N, lambda, u) LfV(S, E, I, R, D, V, xi, sigma, gamma, alpha, rho, N, lambda) + LgV(S, E, I, R, D, V, xi, sigma, gamma, alpha, rho, N, lambda)*u;

%% Time length
TOTAL_TIME = 500; % days
dt = 1;
TOTAL_STEPS = length(0:dt:TOTAL_TIME);

%% Simulate
x = zeros(6,TOTAL_STEPS);
betas = zeros(1, TOTAL_STEPS);
x(:,1) = x0;
maximum_beta = 2.2/6.5;
beta = 2.2/6.5;
D = 0;
c = 0;
for i=2:TOTAL_STEPS
    
    if mod(i, 100) == 0
        lambda = 0.5;
    else
        lambda = 0.0;
    end
    
    cvx_begin quiet

        variables c;
        
        cost = c^2;
        
        minimize(cost);
        subject to
            
            % Get current states
            S = x(1,i-1);
            E = x(2,i-1); 
            I = x(3,i-1);
            R = x(4,i-1);
            D = x(5,i-1);
            Vac = x(6,i-1);
            % Constraints

            Vdot(S, E, I, R, D, Vac, xi, sigma, gamma, alpha, rho, N, lambda, beta) <= -clf_alpha*V(I, E, sigma, gamma);
    cvx_end

    beta = max(beta, 2.2/6.5);
    x(:,i) = x(:,i-1) + dynamics(x(:,i-1), xi, sigma, gamma, alpha, rho, N, lambda, beta)*dt;
    betas(i) = beta;
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

subplot(3, 2, 6)
plot(x(6,:), 'displayname', 'D')
title('V')

figure;
plot(betas)
title('beta');