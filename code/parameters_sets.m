% This file contains the parameters sets in the numerical experiments.

%% Lord2008 test parameters:

% T1 - GBM
% model='GBM'; cp=1; S0=100; K=90; r=0.1; q=0; sigma=0.25;
% modparms = {sigma};

% T2 - VG
% model='VG'; cp=1; S0=100; K=90; r=0.1; q=0; sigma=0.12; theta=-0.14; 
% nu=0.2; t=0.1;
% modparms = {sigma, theta, nu};

% T3 - CGMY
% model='CGMY'; cp=1; S0=1; r=0.1; q=0; C=1; G=5; M=5; Y=0.5;
% t=0.25;

% T4 - CGMY
% model='CGMY'; cp=1; S0=90; r=0.06; q=0; C=0.42; G=4.37; M=191.2; Y=1.0102;
% t=0.25;

%% Other sets of parameters:

% Kou test
% model='KJD'; S0 = 100; K = 98; r = 0.05; q = 0; t = 0.5; cp = 1;
% eta1 = 10; eta2 = 5; lambda = 1; p = 0.4; sigma = 0.16;

% Merton Test
%   S0 = 100; K = 80; T = 0.5;
%   sigma = 0.2; cp = 1; r = 0.0075; q = 0; lambda = 0.01; a = -0.2; 
%   b = 0.6; n = 50;
% model='MJD'; S0 = 100; sigma = 0.2; r = 0.0075; q = 0; lambda = 0.01; mu_y = -0.2; 
%   sigma_y = 0.6;
% modparams = {sigma, mu_y, sigma_y, lambda};
