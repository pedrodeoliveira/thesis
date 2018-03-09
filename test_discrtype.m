% This script tests the convergence of the two discretization methods for
% pricing European call options at various strikes.

clear all, close all

K = 80:10:120;
n = 5:14;
prices = zeros(length(n),length(K));
error = prices;
k_legend = cell(size(K));
k_str = {'gs-','bx-','ro-','kd-','m<-'};

% T2 - VG
model = 'VG'; S0 = 100; r = 0.1; q = 0; sigma = 0.12; 
theta = -0.14; nu = 0.2;
modparams = {sigma, nu, theta};

% Other opttion parameters
t = 1;                              % time to maturity
cp = 1;                             % call (1), put (-1)

% CONV parameters
delta = 40;                     
L = conv_l(model,delta,t,modparams{:});
dt = 2;                                 % discr. type (1) or (2)
alpha = 0;

% Calculate reference value
% vreffcn = @(x) closedf_eurp(model,cp,S0,x,t,r,q,modparams{:});
vreffcn = @(x) cm_fft(20,model,cp,S0,x,t,r,q,modparams{:});
% vreffcn = @(x) bm_eurp(model,'quadgk',cp,S0,x,t,r,q,modparams{:});

% Create pricing function for the conv method
pricefcn = @(x,y) conv_eurp(x,L,alpha,dt,model,cp,S0,y,t,r,q,...
    modparams{:});

% Plot properties
color = [1 1 1];
fontsize = 14;
% linewidth = 1;
h = figure;
set(h,'Color',color)
axesh = axes('Parent',h);
set(axesh,'FontSize',fontsize);

hold on
for j=1:length(K)
    k_legend{j} = sprintf('K = %d',K(j));
    vref = vreffcn(K(j));        
    for i=1:length(n)
        prices(i,j) = pricefcn(n(i),K(j));
        error(i,j) = prices(i,j) - vref;
    end
    if j == 3
        linewidth = 1.5;
    else 
        linewidth = 1;
    end
    plot(n,log10(abs(error(:,j))),k_str{j},'LineWidth',linewidth);   
end
legh = legend(axesh,k_legend);
set(legh,'FontSize',fontsize)
xlabel('n, N = 2^n','FontSize',fontsize)
ylabel('log_{10}|error|','FontSize',fontsize)
