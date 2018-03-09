function y = pcparity(cp,v,S,K,r,q,t)
% PCPARITY Calculates the option value using the Put-Call Parity, if cp=1
%   the function will return the Put value, otherwise it will return the 
%   Call value.

if cp==1
    y = v + S*exp(-q*t) - K*exp(-r*t);
else
    y = v - S*exp(-q*t) + K*exp(-r*t);
end