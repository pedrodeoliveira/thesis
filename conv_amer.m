 function price = conv_amer(n,L,alpha,m,d,model,cp,S0,K,t,r,q,varargin)
%CONV_AMER Calculates the American-style option price using the convolution 
%   (CONV) method, which uses the Repeated Ricardson Extrapolation
%   technique for extrapolating the prices of Bermudan-style options into
%   one of the American-style.

% Lord, R., Fang, F., Bervoets, F., & Oosterlee, C. W. (2008, April). A 
%   fast and accurate fft-based method for pricing early-exercise options 
%   under Levy processes. SIAM J. Sci. Comput., 30(4), 1678-1705.


V = @(dt) conv_berm(n,L,alpha,2,t/dt,model,cp,S0,K,t,r,q,varargin{:});
% h = @(i) t/2^(i-1);
% d = 5;

A = rre(V,'GP',t,m,m,d);

% 4 times RRE (i=1)
% A2 = (1/21)*(64*V(h(4))-56*V(h(3))+14*V(h(2))-V(h(1)));
% price = A2;

price = A(1,m);
