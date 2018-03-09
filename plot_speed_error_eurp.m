% This script tests the 

clc, clear all, close all

% Discretization parameters
n = 7:18;                     % 2^N grid points
method_str = {'bo','rx'};
method_legend = {'CM','CONV'};
Navg = 10;

time = zeros(length(n),length(method_legend));
error = time;

% T5 - KJD
% model='KJD'; S0 = 100; K = 98; r = 0.05; q = 0; t = 0.5; cp = 1;
% eta1 = 10; eta2 = 5; lambda = 1; p = 0.4; sigma = 0.16;
% modparms = {sigma, p, eta1, eta2, lambda};

% T6 - MJD
model='MJD'; S0 = 100; sigma = 0.2; r = 0.0075; q = 0; lambda = 0.01; 
mu_y = -0.2; sigma_y = 0.6; t = 0.5; cp = 1; K = 80;
modparms = {sigma, mu_y, sigma_y, lambda};

% Other opttion parameters
% t = 0.1;                          % time to maturity
% cp = -1;                          % call (1), put (-1)
% K = 110;                          % strike price

% CONV parameters
delta = 40;    
alpha = 0;
L = conv_l(model,delta,t,modparms{:});
dt = 2;                                 % discr. type (1) or (2)

% Create pricing function for the conv method
pricefcn = @(x) conv_eurp(x,L,alpha,dt,model,cp,S0,K,t,r,q,modparms{:});

% Calculate reference value
vref = closedf_eurp(model,cp,S0,K,t,r,q,modparms{:});

% Cycle through the gridsizes array
for k=1:length(n)           
    
    % CarrMadan (CM)
    start_t = cputime;    
    for j=1:Navg
        v = cm_fft(n(k),model,cp,S0,K,t,r,q,modparms{:});
    end
    end_t = cputime;
    time(k,1) = ((end_t - start_t) / Navg) * 1000;
    error(k,1) = v - vref; 
    
    % CONV
    start_t = cputime;    
    for j=1:Navg
        v = pricefcn(n(k));
    end
    end_t = cputime;
    time(k,2) = ((end_t - start_t) / Navg) * 1000;
    error(k,2) = v - vref; 
end

% Plot properties
color = [1 1 1];
fontsize = 14;
linewidth = 1;
h = figure;
set(h,'Color',color)
axesh = axes('Parent',h);
set(axesh,'FontSize',fontsize);

hold on

for j=1:length(method_legend)
    plot(time(:,j),log10(abs(error(:,j))),method_str{j},'LineWidth',...
        linewidth);   
end

legh = legend(axesh,method_legend);
set(legh,'FontSize',fontsize)
xlabel('Time (ms)','FontSize',fontsize)
ylabel('log_{10}|error|','FontSize',fontsize)
axis([0 200 -8 8])
set(axesh,'Box','on')


h2 = figure;
set(h2,'Color',color)
axesh2 = axes('Parent',h2);
set(axesh2,'FontSize',fontsize);

N=2.^n;
plot(N,time(:,2),'o-')
xlabel('Problem Size N','FontSize',fontsize)
ylabel('Time (ms)','FontSize',fontsize)


