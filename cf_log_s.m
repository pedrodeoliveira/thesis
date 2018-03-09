function y = cf_log_s(u,model,logS,T,mu,varargin)
%CF_LOG_S Calculates the characteristic function value for the log-spot
%   price.

switch model
    case {'GBM','MJD','KJD'}
        varargin = [{0} varargin];
end

% CF of process Xt
cf_x = cf(u,model,T,varargin{:});

% Convexity correction
cf_minus_i = cf(-1i,model,T,varargin{:});

% Final value
y = exp(1i*u*(logS + mu*T - log(cf_minus_i))) .* cf_x;
% y = exp(1i*u*(mu*T - log(cf_minus_i))) .* cf_x;