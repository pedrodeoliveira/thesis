function y = closedf_eurp(model,cp,S0,K,t,r,q,varargin)
% CLOSEDF_EURP Calculates the European-style option price using the the 
%   closed-form solution, when available, for the specified model.

switch model
    case 'GBM'
        y = gbm_eurp(cp,S0,K,t,r,q,varargin{:});
    case 'VG'
        y = vg_eurp(cp,S0,K,t,r,q,varargin{:});
    case 'MJD'
        y = merton_eurp(cp,S0,K,t,r,q,varargin{:});
    case 'KJD'
        y = kou_eurp(cp,S0,K,t,r,q,varargin{:});
    case 'CGMY'
        % No closed-form available, using CarrMadan method.
        y = cm_fft(20,model,cp,S0,K,t,r,q,varargin{:});
    otherwise
        disp('Unknown model.')
end

