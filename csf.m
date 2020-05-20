clc; clear all; close all;

%% Constants
xi = 0.1;
sigma = 1/5.1;
gamma = 0.154;
alpha = 0.034;
rho = 1/17.8;
N = 8*1e9;
k1 = 1e-5;
k3 = 0;
%% CLF parameters
scale = 1;
a1 = 1;
a2 = 1;

p1 = (a1+a1^2+a2^2)/(2*a1*a2);
p2 = 1/(2*a1);
p3 = (a1+1)/(2*a1*a2);
P = [p1, p2; 
    p2, p3];

Q = eye(2);
clf_alpha = 1/max(eig(P));
slack_alpha = 0.1;

%% CSF functions
E_max = 1.5*1e9;
% Defining Safety function
h = @(S, Ex, Ix, R, D, M, beta) E_max - Ex;
h_dot = @(S, Ex, Ix, R, D, M, beta) -(S*Ix/N*beta - sigma*Ex);
h_dot_dot = @(S, Ex, Ix, R, D, M, beta) -(1/N*(0.05*S*Ix+xi*beta*R*Ix-beta^2*S*Ix^2/N+beta*sigma*S*Ex-gamma*Ix*S*beta+sigma*beta*S*Ix+sigma^2*Ex*N)-k1*beta*S*Ix*M/N*u);
L_fh = @(S, Ex, Ix, R, D, M, beta) h_dot(S, Ex, Ix, R, D, M, beta);
L_f2h = @(S, Ex, Ix, R, D, M, beta) -(1/N*(0.05*S*Ix+xi*beta*R*Ix-beta^2*S*Ix^2/N+beta*sigma*S*Ex-gamma*Ix*S*beta+sigma*beta*S*Ix+sigma^2*Ex*N));
L_gL_fh = @(S, Ex, Ix, R, D, M, beta) k1*beta*S*Ix*M/N;
H = @(lambda) atan(lambda)+pi/2;%(1+exp(-lambda))^(-1);
dHdlam = @(lambda) 1/(lambda^2+1);%exp(-lambda)/(1+exp(-lambda))^2;
L_fhr = @(S, Ex, Ix, R, D, M, beta) L_fh(S, Ex, Ix, R, D, M, beta)*(1+H(L_fh(S, Ex, Ix, R, D, M, beta)))+h(S, Ex, Ix, R, D, M, beta)*dHdlam(L_fh(S, Ex, Ix, R, D, M, beta))*L_f2h(S, Ex, Ix, R, D, M, beta);
L_ghr = @(S, Ex, Ix, R, D, M, beta) h(S, Ex, Ix, R, D, M, beta)*dHdlam(L_fh(S, Ex, Ix, R, D, M, beta))*L_gL_fh(S, Ex, Ix, R, D, M, beta);
h_r   = @(S, Ex, Ix, R, D, M, beta) h(S, Ex, Ix, R, D, M, beta)*(1+H(L_fh(S, Ex, Ix, R, D, M, beta)));
h_r_dot = @(S, Ex, Ix, R, D, M, beta,u) L_fhr(S, Ex, Ix, R, D, M, beta)+L_ghr(S, Ex, Ix, R, D, M, beta)*u;

csf_beta = 0.05;
%% Intial conditions
I0 = 0.0001*N; % 10 infected
E0 = 0.005*N;
S0 = N - I0 - E0;
R0 = 0;
D0 = 0;
M0 = 1e6;
beta0 = 2.2/6.5;
x0 = [S0, E0, I0, R0, D0, M0, beta0]';

%% Functions
F = @(S, Ix, Ex, R, M, beta) [(-1).*beta.*Ix.*N.^(-1).*S+R.*xi,beta.*Ix.*N.^(-1).*S+(-1).*Ex.* ...
  sigma,(-1).*gamma.*Ix+Ex.*sigma,(1+(-1).*alpha).*gamma.*Ix+(-1).* ...
  R.*xi,alpha.*Ix.*rho,k3,0.5E-1]';

G = @(M, beta) [0,0,0,0,0,(-1).*M,(-1).*beta.*k1.*M]';

eta = @(S, Ex, Ix, beta) [Ex,beta.*Ix.*N.^(-1).*S+(-1).*Ex.*sigma]';
V = @(S, Ex, Ix, beta) eta(S, Ex, Ix, beta)'*P*eta(S, Ex, Ix, beta);

dvdeta = @(S, Ix, Ex, beta) 2*eta(S, Ex, Ix, beta)'*P;
detadx = @(S,Ix,beta) [0,1,0,0,0,0,0;beta.*Ix.*N.^(-1),(-1).*sigma,beta.*N.^(-1).*S,0,0, ...
  0,Ix.*N.^(-1).*S];
dvdx = @(S, Ix, Ex, beta) dvdeta(S, Ix, Ex, beta)*detadx(S,Ix,beta);
LfV = @(S, Ex, Ix, R, D, M, beta) dvdx(S, Ix, Ex, beta)*F(S, Ix, Ex, R, M, beta);
LgV = @(S, Ix, Ex, M, beta) dvdx(S, Ix, Ex, beta)*G(M, beta);

Vdot = @(S, Ex, Ix, R, D, M, beta, u) LfV(S, Ex, Ix, R, D, M, beta) + LgV(S, Ix, Ex, M, beta)*u;

%% Time length
TOTAL_TIME = 100; % days
dt = 0.1;
TOTAL_STEPS = length(0:dt:TOTAL_TIME);

%% Simulate
x = zeros(7,TOTAL_STEPS);
us = zeros(1, TOTAL_STEPS);
x(:,1) = x0;

for i=2:TOTAL_STEPS
    i
    
    if i > TOTAL_STEPS/2 && i < TOTAL_STEPS*3/4
        k3 = 100000;
    end
    cvx_begin quiet

        variables u dh;
        
        cost = u^2 + dh^2;
        
        minimize(cost);
        subject to
            
            % Get current states
            S = x(1,i-1);
            Ex = x(2,i-1); 
            Ix = x(3,i-1);
            R = x(4,i-1);
            D = x(5,i-1);
            M = x(6,i-1);
            beta = x(7,i-1);
            % Constraints

            Vdot(S, Ex, Ix, R, D, M, beta, u) <= -clf_alpha*V(S, Ex, Ix, beta) + dh;
            h_r_dot(S, Ex, Ix, R, D, M, beta,u) >= -csf_beta*h_r(S, Ex, Ix, R, D, M, beta);
            M > 0;
    cvx_end
%     if i < TOTAL_STEPS*0.25
%         u = 0;
%     end
    M
    u = min(0.75, u);
    x(:,i) = x(:,i-1) + dynamics(x(:,i-1), xi, sigma, gamma, alpha, rho, N, k1, k3, u)*dt;
    us(i) = u;
end

subplot(2, 4, 1)
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
title('beta')

subplot(2, 4, 8)
plot(us, 'displayname', 'u')
title('u')