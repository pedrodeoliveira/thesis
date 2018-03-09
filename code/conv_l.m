function L = conv_l(model,delta,t,varargin)
% CONV_L Calculates the value of integration length to use in the
%   convolution (CONV) algorithm for the specified model. The value is
%   calculated using a rule of thumb, i.e., it is a multiple of the 
%   standard deviation of the asset price process.

switch model
    case 'GBM'
        std = gbm_std(t,varargin{:});
    case 'VG'
        std = vg_std(t,varargin{:});
    case 'MJD'
        std = mjd_std(t,varargin{:});
    case 'KJD'
        std = kjd_std(t,varargin{:});
    case 'CGMY'
        std = cgmy_std(t,varargin{:});
    otherwise
        disp('Unknown model.')
end

L = delta*std;

end


%% Explicit implementation of the standard deviation functions
%-----------------------------------------------------------------------

% Geometric Brownian Motion (GBM) - Black-Scholes
function std = gbm_std(t,sigma)
    std = sigma * sqrt(t);
end

% Variance Gamma (VG)
function std = vg_std(t,sigma,theta,nu)
    std = sqrt(t *(sigma^2 + nu*theta^2)); 
end

% Merton Jump-diffusion
function std = mjd_std(t,sigma,mu_y,sigma_y,lambda)
    std = sqrt(t *(sigma^2 + lambda*sigma_y^2 + lambda*mu_y^2)); 
end

% Kou Jump-diffusion
function std = kjd_std(t,sigma,p,eta1,eta2,lambda)
    std = sqrt(t *(sigma^2 + lambda*(p/eta1^2) + lambda*((1-p)/eta2^2))); 
end

% CGMY
function std = cgmy_std(t,C,G,M,Y)
    std = sqrt(t *(C*(M^(Y-2) + G^(Y-2))*gamma(2-Y)));
end

