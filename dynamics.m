function dxdt = dynamics(x, zeta, sigma, gamma, alpha, rho, N, beta)
    dxdt = zeros(5,1);
    
    S = x(1);
    E = x(2);
    I = x(3);
    R = x(4);
    D = x(5);
    
    sdot = zeta*R - S*I/N*beta;
    edot = -sigma*E + S*I/N*beta;
    idot = sigma*E - gamma*I;
    rdot = (1 - alpha)*gamma*I - zeta*R;
    ddot = alpha*rho*I;
    
    dxdt(1) = sdot;
    dxdt(2) = edot;
    dxdt(3) = idot;
    dxdt(4) = rdot;
    dxdt(5) = ddot;
end