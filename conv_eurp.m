function price = conv_eurp(n,L,alpha,dtype,model,cp,S0,K,t,r,q,varargin)
%CONV_EURP Calculates the European-style option price using the convolution 
%   (CONV) method.

% Lord, R., Fang, F., Bervoets, F., & Oosterlee, C. W. (2008, April). A 
%   fast and accurate fft-based method for pricing early-exercise options 
%   under Levy processes. SIAM J. Sci. Comput., 30(4), 1678-1705.


% CREATE THE GRID
N = 2^n;                                        % N gridpoints    
rdt = r * t;                                    % riskfree times t

% Grid steps and index.
Delta_y = L / N;
Delta_x = Delta_y;
Delta_u = (2 * pi) / L;
Grid_i = (0:N-1)';                              % Gridindex
Grid_m = (-1).^Grid_i;                          % (-1)^Grid_i
u = (Grid_i - N/2) * Delta_u;

% DEPENDING ON THE DISCRETIZATION TYPE ADJUST THE GRID
if dtype == 1
    eps_y = 0;
elseif dtype == 2
    eps_y = log(K/S0);                          % Set K on the grid.
end
eps_x = 0;                                      % Set S0 on the grid.
x = eps_x + (Grid_i - N/2) * Delta_x;
y = eps_y + (Grid_i - N/2) * Delta_y;

% Calculate the payoff and dampened the option value.
V = max(cp .* (S0*exp(y) - K), 0);         
v = V .* exp(alpha .* y);                       

% Create the coefficients.
w = ones(N,1); 
w(1) = 0.5; 
w(N) = 0.5;

% Calculate the inner transform and the vector to be transformed.
FT_Vec = ifft( (Grid_m) .* w .* v );            
FT_Vec_tbt = exp( 1i .* Grid_i .* (y(1) - x(1)) .* Delta_u ) ...
    .* feval(@cf_log_s,-(u -(1i*alpha)),model,0,t,r-q,varargin{:}) ...   
    .* FT_Vec;                                  

% Calculate the final value (check value at the given pos.)
% c = abs(exp(-rdt - (alpha .* x) + (1i .* u .* (y(1) - x(1))) ) ...
%     .* (Grid_m) .* fft(FT_Vec_tbt));           

c = abs(exp(-rdt + (1i .* u .* (y(1) - x(1))) ) ...
    .* (Grid_m) .* fft(FT_Vec_tbt));           

% Undampen value.
% C = c;
C = exp(-alpha .* x) .* c;


% Return price
price = double(C(N/2+1, 1));              

