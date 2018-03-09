function price = bm_eurp(model,method,cp,S0,K,T,r,q,varargin)
%BM_EUROP Calculates the European-style option price using the Bakshi-Madan 
%   method.

logS = log(S0);
logK = log(K);

phi = @(u) cf_log_s(u,model,logS,T,r-q,varargin{:});

integrandQs = @(u) real(exp(-1i*u*logK).*phi(u-1i) ./ (1i*u*phi(-1i)));
integrandQ = @(u) real(exp(-1i*u*logK).*phi(u) ./ (1i*u));

Qs = (1/2) + (1/pi) * intgr(integrandQs,method);
Q = (1/2) + (1/pi) * intgr(integrandQ,method);

price = S0*exp(-q*T)*Qs - exp(-r*T)*K*Q;

% Return put value using put-call parity.
if cp~=1
    price = pcparity(cp,price,S0,K,r,q,T);
end

end

% Function that performs the integration using the specified method.
function y = intgr(f,method)
    switch method
        case 'quad'            
            % Simpson's rule
            y = quad(f,0,Inf);
        case 'quadl'
            % Adaptive Lobatto quadrature
            y = quadl(f,0,Inf,1e-8);     
        case 'quadgk'
            % Adaptive Gauss-Kronrod quadrature
            y = quadgk(f,0,Inf,'MaxIntervalCount',1e5);      
    end
            
end