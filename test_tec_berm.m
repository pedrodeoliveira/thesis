% This script tests the CPU time, error and convergence rate of the CONV
% method for Bermudan-style options using Discretization II for a given 
% range of grid sizes. The test
% parameters to consider are T1-GBM and T2-VG. The CPU times are determined 
% after averaging the times of 1.000 experiments.

clear all, close all

n = 7:12;
prices = zeros(length(n),2);
error = prices;
time = prices;
convergence = zeros(length(n),2);
p = convergence;
Navg = 1000;

% T1 - GBM
% model = 'GBM'; S0 = 100; r = 0.1; q = 0; sigma = 0.25;
% delta = 20;
% modparams = {sigma};

% T2 - VG
model = 'VG'; S0 = 100; r = 0.1; q = 0; sigma = 0.12; 
theta = -0.14; nu = 0.2;
delta = 40;
modparams = {sigma, nu, theta};

% Other option parameters
t = 1;                            % time to maturity
cp = -1;                          % call (1), put (-1)
K = 110;                          % strike price
M = 10;

% CONV parameters                 
L = conv_l(model,delta,t,modparams{:});
dt = 2;                                 % discr. type (1) or (2)
alpha = 0;

% Create pricing function for the conv method
pricefcn = @(x) conv_berm(x,L,alpha,dt,M,model,cp,S0,K,t,r,q,modparams{:});

vref = pricefcn(20);
fprintf('Model = %s \t Vref = %.8f\n',model,vref)

for j=1:length(n)
    v = 0;
    start_t = cputime;    
    for k=1:Navg
        v = pricefcn(n(j));
    end
    end_t = cputime;
    time(j,1) = ((end_t - start_t) / Navg) * 1000;
    error(j,1) = v - vref;       
    if j > 1
        convergence(j,1) = abs(error(j-1))/abs(error(j));
        p(j,1) = log(convergence(j,1))/log(2^n(j)/(2^n(j-1)));
    end

%     fprintf('n = %2d, time(ms) = %1.2f, error = %1.2e, conv. = %3.1f, p = %2.1f\n',...
%         n(j),time(j,1),error(j,1),convergence(j,1),p(j,1))
%     fprintf('$%2d$ & $%1.2f$ & $%1.2e$ & $%2.1f$ \n',...
%         n(j),time(j,1),error(j,1),p(j,1))    
    fprintf('& $%1.2f$ & $%1.2e$ & $%2.1f$ \\\\ \n',time(j,1),...
        error(j,1),p(j,1))      
end

