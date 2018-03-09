 function price = conv_berm(n,L,alpha,dtype,M,model,cp,S0,K,t,r,q,varargin)
%CONV_BERM Calculates the Bermudan-style option price using the convolution 
%   (CONV) method.
%   
%   Input arguments:
%
% n           -> points used for the grid construction, N = 2^n
% L           -> width of the interval approx the density
% alpha       -> dampening factor
% dtype       -> discretization type I(=1) or II(=2)
% cp          -> call cp=1; put cp=-1;
% model       -> string for finding model from CF lib
% t           -> maturity
% r           -> riskfree rate
% q           -> dividend yield
% K           -> K of the option
% parameters  -> model parameters
% M           -> number of exercise dates

% Lord, R., Fang, F., Bervoets, F., & Oosterlee, C. W. (2008, April). A 
%   fast and accurate fft-based method for pricing early-exercise options 
%   under Levy processes. SIAM J. Sci. Comput., 30(4), 1678-1705.

% CREATE THE GRID
N = 2^n;                                    % N grid points
dt = t / M;                                 % time step
rdt = r * dt;                               % riskfree times dt

% Grid steps and index.
Delta_y = L / N;
Delta_x = Delta_y;
Delta_u = (2 * pi) / L;
Grid_i = (0:N-1)';                          % Gridindex

% Adjust grid using the specified discretization type

if dtype == 1
    eps_y = 0;
elseif dtype == 2
    % ensure that K is on the grid.
    dm = log(K / S0);                               
    eps_y = dm - ceil(dm / Delta_y) * Delta_y;
end                      
eps_x = 0;                               % ensure that S0 is on the grid.
x = eps_x + (Grid_i - N/2) * Delta_x;
y = eps_y + (Grid_i - N/2) * Delta_y;
u = (Grid_i - N/2) * Delta_u;

% Calculate the option value at t(N) and dampened the option value.
V = max(cp .* (S0*exp(y) - K), 0);         
v = V .* exp(alpha .* y);                       

% Create the coefficients.
w = ones(N,1); 
w(1) = 0.5; 
w(N) = 0.5;

% Create the continuation value handle.
% cval = @(v,x,y) cont_val(Grid_i, rdt, w, v, Delta_u, model, u, x, y, ...
%     alpha, dt, r, q, varargin{:});

% BACKWARD INDUCTION TO CALCUALTE OPTION PRICES
for m = M-1:-1:1         
    
    % Calculate the continuation value using (27)
    % Inner transform
    FT_Vec = ifft( ((-1) .^ Grid_i) .* w .* v );
    % Outer transform
    FT_Vec_tbt = exp( 1i .* Grid_i .* (y(1) - x(1)) .* Delta_u ) ...
        .* feval(@cf_log_s,-(u-(1i*alpha)),model,0,dt,r-q,varargin{:})...
        .* FT_Vec;
    C = abs(exp(-rdt-(alpha .* x) + (1i .* u .* (y(1) - x(1))) ) ...
        .* ((-1).^Grid_i) .* fft(FT_Vec_tbt));
    
    % Calculate the exercise (payoff) at t_m
    E = max(cp .* (S0*exp(x) - K), 0);
    % The option value is the maxium of the previous values.
    V = max(C, E);        

    % Locate the exercise boundary j
    j = find(V==E,1,'last');               
    if isempty(j)
        if cp == 1
            j = size(E,1) - 1;
        else
            j = 1;
        end
    end
    if cp == 1
        j = j - 1;
    end

    % Approximate dm using (45)
    dm = ( x(j+1)*(C(j) - E(j)) - x(j)*(C(j+1) - E(j+1)) ) ...
        / ( (C(j) - E(j)) - (C(j+1) - E(j+1)) );

    % Recalculate the x-grid
    % QUESTION: eps_dm = dm or replace the new dm in eps_x equation??
    eps_x = dm - (ceil(dm / Delta_x) * Delta_x);
    x = eps_x + (Grid_i - N/2) * Delta_x;
    
    % Recompute the continuation value, the inner transform does not need
    % to be recalculated.
    FT_Vec_tbt = exp( 1i .* Grid_i .* (y(1) - x(1)) .* Delta_u ) ...
        .* feval(@cf_log_s,-(u-(1i*alpha)),model,0,dt,r-q,varargin{:})...
        .* FT_Vec;
    C = abs(exp(-rdt-(alpha .* x) + (1i .* u .* (y(1) - x(1))) ) ...
        .* ((-1).^Grid_i) .* fft(FT_Vec_tbt));
    E = max(cp .* (S0*exp(x) - K), 0);
    V = max(C, E);

    % Set y-grid = x-grid
    y = x;                                  
    v = V .* exp(alpha .* y);               % dampened option value
end

% LAST STEP: 
% Set eps_x such thath S0 is on the grid and calculate the final value.
eps_x = 0;
x = eps_x + (Grid_i - N/2) * Delta_x;
FT_Vec = ifft( ((-1) .^ Grid_i) .* w .* v );
% Outer transform
FT_Vec_tbt = exp( 1i .* Grid_i .* (y(1) - x(1)) .* Delta_u ) ...
        .* feval(@cf_log_s,-(u-(1i*alpha)),model,0,dt,r-q,varargin{:})...
        .* FT_Vec;
C = abs(exp(-rdt-(alpha .* x) + (1i .* u .* (y(1) - x(1))) ) ...
    .* ((-1).^Grid_i) .* fft(FT_Vec_tbt));                       
price = C(N/2 + 1, 1);                  

