% This script tests the effect of the choice of "L" on the error of the 
% CONV method for European-style options using as reference the CarrMadan 
% method. Plots for different values of n are given. Model used  T6 - KJD.

clc, clear all, close all

L = 0.1:0.1:10;                     % choose L
num = length(L);                    % get the size of L array
prices = zeros(size(L));            % init the price vector
error = prices;

% Discretization parameters
n = 5:4:13;                     % 2^N grid points
n_str = {'b--','ro-','k-'};
n_legend = cell(size(n));

% T6 - KJD
model='KJD'; S0 = 100; K = 98; r = 0.05; q = 0; t = 0.5; cp = 1;
eta1 = 10; eta2 = 5; lambda = 1; p = 0.4; sigma = 0.16;
modparms = {sigma, p, eta1, eta2, lambda};

% Other opttion parameters
% t = 0.1;                          % time to maturity
% cp = -1;                          % call (1), put (-1)
% K = 110;                          % strike price

% CONV parameters
% delta = 40;                     
% L = conv_l(model,delta,t,modparms{:});
dt = 2;                                 % discr. type (1) or (2)
alpha=0;
delta = 40;
L_RT = conv_l(model,delta,t,modparms{:});

% Create pricing function for the conv method
pricefunc = @(x,y) conv_eurp(x,y,alpha,dt,model,cp,S0,K,t,r,q,modparms{:});

% Calculate reference value
% vref = closedf_eurp(model,cp,S0,K,t,r,q,modparams{:});
vref = cm_fft(20,model,cp,S0,K,t,r,q,modparms{:});
% vref = pricefunc(20,0);
% vref = bm_eurp(model,'quadgk',cp,S0,K,t,r,q,modparams{:});

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
        prices(j) = pricefunc(n(k),L(j));
        error(j) = prices(j) - vref;
    end
    if k == 3
        linewidth = 1.5;
    end
    plot(L,log10(abs(error)),n_str{k},'LineWidth',linewidth);       
end

for k=-10:0.1:10
    plot(L_RT,k,'-');    
end
text(L_RT+0.1,-9,sprintf('L_{RT} = %1.2f',L_RT),'FontSize',fontsize)

legh = legend(axesh,n_legend);
set(legh,'FontSize',fontsize)
xlabel('L','FontSize',fontsize)
ylabel('log_{10}|error|','FontSize',fontsize)
axis([0 10 -10 10])

