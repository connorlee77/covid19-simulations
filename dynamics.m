function dxdt = dynamics(x, xi, sigma, gamma, alpha, rho, N, k, u)
    dxdt = zeros(5,1);
    
    S = x(1);
    E = x(2);
    I = x(3);
    R = x(4);
    D = x(5);
    M = x(6);
    beta = x(7);
    
    sdot = xi*R - S*I/N*beta;
    edot = -sigma*E + S*I/N*beta;
    idot = sigma*E - gamma*I;
    rdot = (1 - alpha)*gamma*I - xi*R;
    ddot = alpha*rho*I;
    mdot = (1-beta)-u*M;
    betadot = -(k*u*M - 0.5);
    
    dxdt(1) = sdot;
    dxdt(2) = edot;
    dxdt(3) = idot;
    dxdt(4) = rdot;
    dxdt(5) = ddot;
    dxdt(6) = mdot;
    dxdt(7) = betadot;
end