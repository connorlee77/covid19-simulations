function beta = feedback(x, sigma, gamma, alpha_1,alpha_2, N)

beta = N/x(1)/x(3)/sigma*(sigma^2*x(2)+gamma*sigma*x(2)-gamma^2*x(3)-alpha_1*x(3)-alpha_2*(sigma*x(2)-gamma*x(3)));
    
end