% This script tests the effect of the dampening factor "alpha" parameter on
% the error of the CONV method for Bermudan-style options using as 
% reference the CONV method with 2^20 points. Plots for different values 
% of n are given.

clc, clear all, close all

alpha = -10:0.25:10;                % choose dampening parameter
num = length(alpha);                % get the size of alpha array
prices = zeros(size(alpha));        % init the price vector
error = prices;

% Discretization parameters
n = 5:4:13;                     % 2^N grid points
n_str = {'b--','ro-','k-'};
n_legend = cell(size(n));

% T2 - VG
model = 'VG'; S0 = 100; r = 0.1; q = 0; sigma = 0.12; 
theta = -0.14; nu = 0.2;
modparams = {sigma, nu, theta};

% Other opttion parameters
t = 1;                          % time to maturity
cp = -1;                          % call (1), put (-1)
K = 110;                          % strike price
M = 10;                           % Number of exercise opportunities.

% CONV parameters
delta = 40;                     
L = conv_l(model,delta,t,modparams{:});
dt = 2;                                 % discr. type (1) or (2)

% Create pricing function for the conv method
pricefunc = @(x,y) conv_berm(x,L,y,dt,M,model,cp,S0,K,t,r,q,modparams{:});

% Calculate reference value
vref = pricefunc(20,0);

% Plot properties
color = [1 1 1];
fontsize = 14;
linewidth = 1;
h = figure;
set(h,'Color',color)
axesh = axes('Parent',h);
set(axesh,'FontSize',fontsize);

hold on
% Cycle through the gridsizes array
for k=1:length(n)
    % Cycle through alpha factors    
    n_legend{k} = sprintf('2^{%d}',n(k));
    for j = 1:num
        prices(j) = pricefunc(n(k),alpha(j));
        error(j) = prices(j) - vref;
    end
    if k == 3
        linewidth = 2;
    end
    plot(alpha,log10(abs(error)),n_str{k},'LineWidth',linewidth);   
end

legh = legend(axesh,n_legend);
set(legh,'FontSize',fontsize)
xlabel('\alpha','FontSize',fontsize)
ylabel('log_{10}|error|','FontSize',fontsize)
axis([-10 10 -7 7])

