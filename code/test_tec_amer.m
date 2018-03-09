% This script tests the CPU time, error and convergence rate of the CONV
% method for American-style options using Discretization II for a given 
% range of grid sizes. The test
% parameters to consider are T1-GBM and T2-VG. The CPU times are determined 
% after averaging the times of 1.000 experiments.

clear all, close all

n = 7:12;
prices = zeros(length(n),2);
error = prices;
time = prices;
convergence = prices;
p = prices;
Navg = 10;

% T1 - GBM
% model = 'GBM'; S0 = 100; r = 0.1; q = 0; sigma = 0.25;
% delta = 20;
% modparams = {sigma};

% T2 - VG
model = 'VG'; S0 = 100; r = 0.1; q = 0; sigma = 0.12; 
theta = -0.14; nu = 0.2;
delta = 40;
modparams = {sigma, nu, theta};

% Other option parameters T1 and T2
t = 1;                            % time to maturity
cp = -1;                          % call (1), put (-1)
K = 110;                          % strike price

% T3 - CGMY
% model='CGMY'; K = 1; S0=1; r=0.1; q=0; C=1; G=5; M=5; Y=0.5;
% modparams = {C, G, M, Y};

% T4 - CGMY
model='CGMY'; S0=90; r=0.06; q=0; C=0.42; G=4.37; M=191.2; Y=1.0102;
t=0.25; K = 98;
modparams = {C, G, M, Y};


% CONV parameters                 
L = conv_l(model,delta,t,modparams{:});
dt = 2;                                 % discr. type (1) or (2)
alpha = 0;

% Extrapolation
m = 3;                       % number of repetitions
d = 5;

% Create pricing function for the conv method
bermfcn = @(x,y) conv_berm(x,L,alpha,dt,y,model,cp,S0,K,t,r,q,modparams{:});
amerfcn = @(x,y) conv_amer(x,L,alpha,m,y,model,cp,S0,K,t,r,q,...
    modparams{:});

vref = amerfcn(14,7);
fprintf('Model = %s \t Vref = %.8f\n',model,vref)

for j=1:length(n)
    v = 0;
%     start_t = cputime;    
%     for k=1:Navg
%         v = bermfcn(n(j),2^(n(j)-1));
%     end
%     end_t = cputime;
%     time(j,1) = ((end_t - start_t) / Navg) * 1000;
%     error(j,1) = v - vref;       
%     if j > 1
%         convergence(j,1) = abs(error(j-1))/abs(error(j));
%         p(j,1) = log(convergence(j,1))/log(2^n(j)/(2^n(j-1)));
%     end
    
%     fprintf('n = %2d, time(ms) = %1.2f, error = %1.2e, conv. = %3.1f, p = %2.1f\t',...
%         n(j),time(j,1),error(j,1),convergence(j,1),p(j,1))   
%     fprintf('$%2d$ & $%1.2f$ & $%1.2e$ & $%2.1f$ ',...
%         n(j),time(j,1),error(j,1),p(j,1))  
    
    start_t = cputime;    
    for k=1:Navg
        v = amerfcn(n(j),d);
    end
    end_t = cputime;
    time(j,2) = ((end_t - start_t) / Navg) * 1000;
    error(j,2) = v - vref;       
    if j > 1
        convergence(j,2) = abs(error(j-1,2))/abs(error(j,2));
        p(j,2) = log(convergence(j,2))/log(2^n(j)/(2^n(j-1)));
    end    
%     fprintf('time(ms) = %1.2f, error = %1.2e, conv. = %3.1f, p = %2.1f\n',...
%         time(j,2),error(j,2),convergence(j,2),p(j,2))   
%     fprintf('& $%1.2f$ & $%1.2e$ & $%2.1f$ \\\\ \n',time(j,2),...
%         error(j,2),p(j,2))  
    
    fprintf('& $%1.2f$ & $%1.2e$ \n',time(j,2),error(j,2))      
end

