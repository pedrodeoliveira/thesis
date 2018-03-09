function cf = cf(u,model,T,varargin)
%---------------------------------------------------------
% Characteristic Function Library for the Levy process Xt
% for the following models:
%---------------------------------------------------------
% GBM - Geometric Brownian Motion or Black-Scholes Model
% MJD - Merton Jump-diffusion
% KJD - Kou Jump-diffusion
% VG - Variance Gamma
% CGMY
%---------------------------------------------------------

% Exceptions definitions
MEargs = MException('VerifyInput:InvalidNrOfArguments',...
    'Invalid number of Input arguments');

% Models: names, functions and number of expected args.
models.names = {'GBM', 'MJD', 'KJD', 'VG', 'CGMY'};
models.funcs = {@ce_gbm, @ce_merton, @ce_kou, @ce_vg, @ce_cgmy};
models.nbrargs = [2, 5, 6, 3, 4];

% Iterate over list of models.
for i=1:length(models.names)
    if strcmp(model, models.names{i})
        if length(varargin(:)) == models.nbrargs(i)
%         if nargin == 5
            funobj = models.funcs{i};
            ce = feval(funobj,u,varargin{:}); 
            cf = exp(T*ce);
            return;
        else 
            throw(MEargs)
        end
    end
end

MEmodel = MException('VerifyInput:InvalidModel',...
    sprintf('Undefined Model: %s', model));
throw(MEmodel)

end

%% Explicit implementation of the characteristic exponent
%-----------------------------------------------------------------------

% Geometric Brownian Motion (GBM) - Black-Scholes
function y = ce_gbm(u,mu,sigma)
    y = 1i*mu*u - 0.5*sigma*sigma*u.*u;
end

% Merton Jump-diffusion
function y = ce_merton(u,mu,sigma,mu_y,sigma_y,lambda)
    y = ce_gbm(u,mu, sigma) + lambda*ce_lognormal(u,mu_y,sigma_y);
end

% Kou Jump-diffusion
function y = ce_kou(u,mu,sigma,p,eta1,eta2,lambda)
    y = ce_gbm(u,mu, sigma) + lambda*ce_double_exp(u,p,eta1,eta2);
end

% Variance Gamma (VG)
function y = ce_vg(u,sigma,nu,theta)
    y = (-1/nu)*(log(1 - 1i*theta*nu*u + 0.5*sigma*sigma.*u.*u*nu));
end

% CGMY
function y = ce_cgmy(u,C,G,M,Y)
    y = C*gamma(-Y)*((M-1i*u).^Y - M^Y + (G+1i*u).^Y - G^Y);
end

% Lognormal
function y = ce_lognormal(u,mu_y,sigma_y)
    y = exp(1i*mu_y.*u - 0.5*sigma_y*sigma_y.*u.*u) - 1;
end

% Double-exponential
function y = ce_double_exp(u,p, eta1, eta2)
    y = 1i.*u.*( (p./(eta1 - 1i.*u)) - ((1-p)./(eta2 + 1i.*u)));
end
