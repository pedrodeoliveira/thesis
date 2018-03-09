function price = kou_eurp(cp,S0,K,t,r,q,sigma,p,eta1,eta2,lambda)
%KOU_EURP Calculates the European-style option price for the Kou model 
%   using the closed-form solution.

% References:
%   Kou, S. G. (2002, August). A jump-diffusion model for option pricing. 
%       Manage. Sci., 48(8), 1086-1101.
%

qsi = ((p*eta1)/(eta1-1)) + (((1-p)*eta2)/(eta2+1)) - 1;
p_tilde = (p/(1+qsi))*(eta1/(eta1-1));
eta1_tilde = eta1-1;
eta2_tilde = eta2+1;
lambda_tilde = lambda*(qsi+1);

mu = r-q + (0.5*sigma^2) - lambda*qsi;
a = log(K/S0);
Qs = upsilon(mu,sigma,lambda_tilde,p_tilde,eta1_tilde,eta2_tilde,a,t);

mu = r-q - (0.5*sigma^2) - lambda*qsi;
Q = upsilon(mu,sigma,lambda,p,eta1,eta2,a,t);

c = S0*exp(-q*t)*Qs - K*exp(-r*t)*Q;

if cp==1
    price = c;
else
    price = pcparity(cp,c,S0,K,r,q,t);
end

end

% Calculates the upsilon or P(Z(t) >= a).
function y = upsilon(mu,sigma,lambda,p,eta1,eta2,a,t)

N=10;
y2 = zeros(N,1);
y3 = y2;
y3aux = y2;

% 1st summation term.
y1 = exp(0.5*t*(sigma*eta1)^2)/(sigma*sqrt(2*pi*t));
y2 = poisspdf(1:N,lambda*t);
for n=1:N
    for k=1:n
        y3aux(k) = P(n,k,eta1,eta2,p)*((sigma*sqrt(t)*eta1).^k)*...
            I(k-1,a-mu*t,-eta1,-(1/(sigma*sqrt(t))),-sigma*eta1*sqrt(t));
    end
    y3(n) = y2(n)*sum(y3aux);
end
y4 = sum(y3);

% 2nd summation term.
y5 = exp(0.5*t*(sigma*eta2)^2)/(sigma*sqrt(2*pi*t));
y6 = zeros(N,1);
y6aux = y6;
for n=1:N
    for k=1:n
        y6aux(k) = Q(n,k,eta1,eta2,p)*((sigma*sqrt(t)*eta2).^k)*...
            I(k-1,a-mu*t,eta2,(1/(sigma*sqrt(t))),-sigma*eta2*sqrt(t));
    end
    y6(n) = y2(n)*sum(y6aux);
end
y7 = sum(y6);

% 3th summation term.
y8 = poisspdf(0,lambda*t)*normcdf(-(a-mu*t)/(sigma*sqrt(t)));

% Finally, add all terms.
y = y1*y4 + y5*y7 + y8;

end

% Calculates the value for the Hh function.
function y = Hh(n,x)

if n==-1
%     y = sqrt(2*pi)*normcdf(x);
    y = exp(-0.5*x^2);
elseif n==0
    y = sqrt(2*pi)*normcdf(-x);
else
    y = (Hh(n-2,x) - x*Hh(n-1,x))/n;
end

end

% Calculates the value of the integration In.
function y = I(n,c,alpha,beta,delta)

y1 = -(exp(alpha*c)/alpha);
y1aux = zeros(n+1,1);
for j=0:n
    y1aux(j+1) = ((beta/alpha).^(n-j)) * Hh(j,beta*c-delta);
end
y1 = y1 * sum(y1aux);
y2 = ((beta/alpha)^(n+1))*(sqrt(2*pi)/beta)*exp((alpha*delta/beta) +...
    ((alpha^2)/(2*beta^2)));

if beta>0 && alpha~=0
    y = y1 + y2*normcdf(-beta*c + delta + (alpha/beta));
elseif beta<0 && alpha<0
    y = y1 - y2*normcdf(beta*c - delta - (alpha/beta));
else
    disp('AHHHHHHHHHHH')
    y = NaN;
end
    
end

% Calculates the Pnk value.
function y = P(n,k,eta1,eta2,p)

q=1-p;

if k==n
    y = p^n;
else
    aux = zeros(n-1-k,1);
    for j=k:n-1
        aux(j) = nchoosek(n-k-1,j-k)*nchoosek(n,j)*((eta1/(eta1+eta2))...
            .^(j-k))*((eta2/(eta1+eta2)).^(n-j))*(p.^j)*(q.^(n-j));
    end
    y = sum(aux);
end

end

% Calculates the Qnk value.
function y = Q(n,k,eta1,eta2,p)

q=1-p;

if k==n
    y = q^n;
else
    aux = zeros(n-1-k,1);
    for j=k:n-1    
        aux(j) = nchoosek(n-k-1,j-k)*nchoosek(n,j)*((eta1/(eta1+eta2))...
            .^(n-j))*((eta2/(eta1+eta2)).^(j-k))*(p.^(n-j))*(q.^j);
    end
    y = sum(aux);
end

end
