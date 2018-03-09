function price = vg_eurp(cp,S0,K,t,r,q,sigma,nu,theta)
%VG_EURP Calculates the European-style option price for Variance Gamma 
%   model using the closed-form solution.

% References:
%   Madan, D. B., Carr, P., & Chang, E. C. (1998). The variance gamma 
%       process and option pricing. European Finance Review, 2, 79?105.
%

% From (16)-(17)
qsi = - theta/sigma^2;
s = sigma/sqrt(1+0.5*nu*(theta/sigma)^2);
alpha = qsi*s;

% From (26)-(28)
c1 = 0.5*nu*(alpha+s)^2;
c2 = 0.5*nu*alpha^2;
d = (1/s)*(log(S0/K) + r*t + (t/nu)*log((1-c1)/(1-c2)));

a1 = d*sqrt((1-c1)/nu);
b1 = (alpha+s)*sqrt(nu/(1-c1));
g = t/nu;
Qs = psi(a1,b1,g);

a2 = d*sqrt((1-c2)/nu);
% b2 = (alpha*s)*sqrt(nu/(1-c2));
b2 = alpha*sqrt(nu/(1-c2));
Q = psi(a2,b2,g);

c = S0*exp(-q*t)*Qs - K*exp(-r*t)*Q;

if cp==1
    price = c;
else
    price = pcparity(cp,c,S0,K,r,q,t);
end

end


% Equation (A11)
function y = psi(a,b,g)

% Auxiliar variables.
c = abs(a)*sqrt(2+b^2);
u = b/sqrt(2+b^2);

y1 = (c^(g+0.5) * exp(sign(a)*c)*(1+u)^g)/(sqrt(2*pi)*gamma(g)*g);
y2 = besselk(g+0.5,c);
y3 = phi(g,1-g,1+g,0.5*(1+u),-sign(a)*c*(1+u));
y4 = sign(a)*((c^(g+0.5)*exp(sign(a)*c)*(1+u)^(g+1))/(sqrt(2*pi)*gamma(g)*(g+1)));
y5 = besselk(g-0.5,c);
y6 = phi(1+g,1-g,2+g,0.5*(1+u),-sign(a)*c*(1+u));
y7 = sign(a)*((c^(g+0.5)*exp(sign(a)*c)*(1+u)^g)/(sqrt(2*pi)*gamma(g)*g));

y = y1*y2*y3 - y4*y5*y6 + y7*y5*y3;

end


% Calculates the Degenerate Hypergeometric function
function y = phi(a,b,g,x,y)

phi_int=@(u) (u.^(a-1)).*((1-u).^(g-a-1)).*((1-u.*x).^(-b)).*exp(u.*y);
y = (gamma(g)/(gamma(a)*gamma(g-a)))*quadgk(phi_int,0,1,'AbsTol',1e-8,...
    'RelTol',1e-8);

end

