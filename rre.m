function A = rre(V,pt,t,I,k,d)
% RRE Performs the k-1 Repeated Richardson Extrapolation on the the given
% function v using I points from v.

% Chang, C.-C., Chung, S.-L., & Stapleton, R. C. (2007). Richardson 
%   extrapolation techniques for the pricing of American-style options. 
%   Journal of Futures Markets, 27(8), 791?817.

if strcmp(pt,'AP') == 1
    h = @(i) t/(d+i);
elseif strcmp(pt,'GP') == 1
    h = @(i) t/2^(d+i-1);
end

A = zeros(max(I,k),k);

for m=0:k-1
    for i=1:(size(A,1)-m)
        A(i,m+1) = rre_step(A,V,h,i,m);
    end
end

end

% The recursive function which calculates each entry of the matrix.
function a = rre_step(A,V,h,i,m)

if m == 0
    if A(i,m+1) == 0
        a = V(h(i));
    else
        a = A(i,m+1);
    end
else
    A(i+1,m) = rre_step(A,V,h,i+1,m-1); 
    A(i,m) = rre_step(A,V,h,i,m-1); 
    a = A(i+1,m) + (A(i+1,m) - A(i,m))/((h(i)/h(i+m))-1);
end

end

