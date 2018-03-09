function price = merton_eurp(cp,S0,K,t,r,q,sigma,mu_y,sigma_y,lambda)
%MERTON_EURP Calculates the European-style option price for the Merton 
%   model using the closed-form solution.

% References:
%   Merton, R. C. (1976). Option pricing when underlying stock returns are 
%       discontinuous. Journal of Financial Economics, 3(1-2), 125-144.
%

J=170;
c = zeros(J,1);

for j=0:J-1    
    lambda_j = lambda*exp(mu_y + (0.5*sigma_y^2));
    r_j = r-lambda*(exp(mu_y+0.5*sigma_y^2)-1)+((j*(mu_y+0.5*sigma_y^2))/t);
    sigma_j = sqrt(sigma^2 + (j*sigma_y^2)/t);    
    bs = gbm_eurp(1,S0,K,t,r_j,q,sigma_j);
    c(j+1) = (1/factorial(j))*((lambda_j*t)^j)*exp(-lambda_j*t)*bs;        
end

price = sum(c);

% Return put value using put-call parity.
if cp~=1
    price = pcparity(cp,price,S0,K,r,q,t);
end

