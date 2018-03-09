% This script tests the effect of the choice of "L" on the error of the 
% CONV method for Bermudan-style options using as reference the CONV method 
% with n=20. Plots for different values of n are given. Model used  
% T1 - GBM.

clc, clear all, close all

L = 0.1:0.1:10;                     % choose L
num = length(L);                    % get the size of L array
prices = zeros(size(L));            % init the price vector
error = prices;

% Discretization parameters
n = 5:4:13;                     % 2^N grid points
n_str = {'b--','ro-','k-'};
n_legend = cell(size(n));

% T1 - GBM
model='GBM'; cp=1; S0=100; K=90; r=0.1; q=0; sigma=0.25;
modparms = {sigma};

% Other opttion parameters
t = 1;                          % time to maturity
M = 10;

% CONV parameters
% delta = 40;                     
% L = conv_l(model,delta,t,modparms{:});
dt = 2;                                 % discr. type (1) or (2)
alpha=0;
delta = 20;
L_RT = conv_l(model,delta,t,modparms{:});

% Create pricing function for the conv method
pricefcn = @(x,y) conv_berm(x,y,alpha,dt,M,model,cp,S0,K,t,r,q,...
    modparms{:});

% Calculate reference value
vref = pricefcn(20,L_RT);

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
        prices(j) = pricefcn(n(k),L(j));
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
% saveas(h,'figures/eurp_VG_alpha_errors.eps')

