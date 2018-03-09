% This script tests the effect of the dampening factor "alpha" parameter on
% the error of the CONV method for European-style options using as 
% reference the CarrMadan method for a given n=9 and for a set of models
% being studied GBM,VG,CGMY,MJD and KJD.

clc, clear all, close all

alpha = -10:0.25:10;            % choose dampening parameter
num = length(alpha);            % get the size of alpha array
prices = zeros(size(alpha));    % init the price vector
error = prices;

% Discretization parameters
n = 9;
models = cell(5,1);
models_str = {'m-','g--','k-.','bo-','rx-'};
n_legend = cell(size(models));

% Other option parameters
t = 0.1;                          % time to maturity
cp = -1;                          % call (1), put (-1)
K = 110;                          % strike price

S0=100;                           % spot price
r=0.1;
q=0;

% CONV parameters             
dt = 2;                                 % discr. type (1) or (2)

% T1 - GBM
model = 'GBM'; sigma=0.25; delta = 20;
modparms = {sigma};
models{1}.name=model;
models{1}.parms = modparms;
models{1}.L = conv_l(model,delta,t,modparms{:});

delta = 40;

% T2 - VG
model = 'VG'; sigma = 0.12; theta = -0.14; nu = 0.2;
modparms = {sigma, nu, theta};
models{2}.name=model;
models{2}.parms = modparms;
models{2}.L = conv_l(model,delta,t,modparms{:});

% T3 - CGMY
model='CGMY'; C=1; G=5; M=5; Y=0.5;
modparms = {C, G, M, Y};
models{3}.name=model;
models{3}.parms = modparms;
models{3}.L = conv_l(model,delta,t,modparms{:});

% T5 - KJD
model='KJD'; eta1 = 10; eta2 = 5; lambda = 1; p = 0.4; sigma = 0.16;
modparms = {sigma, p, eta1, eta2, lambda};
models{4}.name=model;
models{4}.parms = modparms;
models{4}.L = conv_l(model,delta,t,modparms{:});

% T6 - MJD
model='MJD'; sigma = 0.2; lambda = 0.01; mu_y = -0.2;  sigma_y = 0.6;
modparms = {sigma, mu_y, sigma_y, lambda};
models{5}.name=model;
models{5}.parms = modparms;
models{5}.L = conv_l(model,delta,t,modparms{:});

% Plot properties
color = [1 1 1];
fontsize = 14;
linewidth = 2;
h = figure;
set(h,'Color',color)
axesh = axes('Parent',h);
set(axesh,'FontSize',fontsize);

hold on
% Cycle through the gridsizes array
for k=1:length(models)    
    model = models{k}.name;
    n_legend{k} = model;
    modparms = models{k}.parms;
    vref = cm_fft(20,model,cp,S0,K,t,r,q,modparms{:});
    L = models{k}.L;
    pricefcn = @(y) conv_eurp(n,L,y,dt,model,cp,S0,K,t,r,q,modparms{:});
    % Cycle through alpha factors 
    for j = 1:num        
        prices(j) = pricefcn(alpha(j));
        error(j) = prices(j) - vref;
    end
    if k > 3
        linewidth = 1;
    end
    plot(alpha,log10(abs(error)),models_str{k},'LineWidth',linewidth);   
end

legh = legend(axesh,n_legend);
set(legh,'FontSize',fontsize)
xlabel('\alpha','FontSize',fontsize)
ylabel('log_{10}|error|','FontSize',fontsize)
axis([-10 10 -7 7])

