function u = feedback(x, sigma, gamma, xi, alpha_1,alpha_2, N,k1, k3)
    %beta = N/x(1)/x(3)/sigma*(sigma^2*x(2)+gamma*sigma*x(2)-gamma^2*x(3)-alpha_1*x(3)-alpha_2*(sigma*x(2)-gamma*x(3)));  
    S = x(1);
    E = x(2);
    I = x(3);
    R = x(4);
    D = x(5);
    M = x(6);
    beta = x(7);
    
    eta_1 = E;
    eta_2 = beta*S*I/N-sigma*E;
    
    LfLg = -k1*beta*S*I*M/N;
    Lf2 = 1/N*(0.05*S*I+xi*beta*R*I-beta^2*S*I^2/N+beta*sigma*S*E-gamma*I*S*beta+sigma*beta*S*I+sigma^2*E*N);
    u = LfLg^(-1)*(Lf2-alpha_1*eta_1-alpha_2*eta_2);
end