% This script tests the 

clc, clear all

prices = zeros(5,5);    % init the price vector
time = prices;
error = zeros(5,4);

% Discretization parameters
n = 9;
models = cell(5,1);
Navg = 10;

% Other option parameters
t = 0.1;                          % time to maturity
cp = -1;                          % call (1), put (-1)
K = 110;                          % strike price

S0=100;                           % spot price
r=0.1;
q=0;

% CONV parameters             
dt = 2;                                 % discr. type (1) or (2)
alpha = 0;

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

% Cycle through the models
for k=1:length(models)    
    model = models{k}.name;    
    modparms = models{k}.parms;
    
    % Closed-form
    start_t = cputime;    
    for j=1:Navg
        vref = closedf_eurp(model,cp,S0,K,t,r,q,modparms{:});
    end
    end_t = cputime;
    time(k,1) = ((end_t - start_t) / Navg) * 1000;    
    prices(k,1) = vref;
    
    % Bakshi-Madan (Inversion of CF)
    start_t = cputime;    
    for j=1:Navg
        v = bm_eurp(model,'quadgk',cp,S0,K,t,r,q,modparms{:});
    end
    end_t = cputime;
    time(k,2) = ((end_t - start_t) / Navg) * 1000;      
    prices(k,2) = v;
    error(k,1) = v - vref;
    
    % CarrMadan (FFT)
    start_t = cputime;    
    for j=1:Navg
        v = cm_fft(n,model,cp,S0,K,t,r,q,modparms{:});
    end
    end_t = cputime;
    time(k,3) = ((end_t - start_t) / Navg) * 1000;        
    prices(k,3) = v;
    error(k,2) = v - vref;
    
    % CONV
    L = models{k}.L;
    start_t = cputime;    
    for j=1:Navg
        v = conv_eurp(n,L,alpha,dt,model,cp,S0,K,t,r,q,modparms{:});
    end
    end_t = cputime;
    time(k,4) = ((end_t - start_t) / Navg) * 1000;        
    prices(k,4) = v;
    error(k,3) = v - vref;  
    
    % CONV n=20
    L = models{k}.L;
    start_t = cputime;    
    for j=1:Navg
        v = conv_eurp(20,L,alpha,dt,model,cp,S0,K,t,r,q,modparms{:});
    end
    end_t = cputime;
    time(k,5) = ((end_t - start_t) / Navg) * 1000;        
    prices(k,5) = v;
    error(k,4) = v - vref;    
end

methods = {'CF','BM','CM','CONV','CONV(n=20)'};
% Print prices, times and errors for each method/model
for k=1:length(models)   
    model = models{k}.name;
    fprintf('%s & ',model)
    for j=1:length(methods)-1
        method = methods{j};
        if j==1
            fprintf('$%3.2f$ & $%10.7f$ &',time(k,j),...
                prices(k,j))
        elseif j==length(methods)-1
            fprintf('$%3.2f$ & $%3.2e$ \\\\ \n',time(k,j),...
                error(k,j-1))
        else
            fprintf('$%3.2f$ & $%3.2e$ &',time(k,j),...
                error(k,j-1))            
        end
    end
end




