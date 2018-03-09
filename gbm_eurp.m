function price = gbm_eurp(cp,S0,K,t,r,q,sigma)
% GBM_EURP Calculates the European-style option price for the Geometric 
%   Brownian Motion using the closed-form solution.

% Calculate d1.
d1 = (log(S0 / K) + (r - q + 0.5*(sigma^2))*t) / (sigma * sqrt(t));

% Calculate d2.
d2 = d1 - sigma * sqrt(t);

% Calculate call value.
c = (S0 * exp(-q * t) * normcdf(d1)) - (K * exp(-r * t) * normcdf(d2));

if cp==1
    price = c;
else
    price = pcparity(cp,c,S0,K,r,q,t);
end