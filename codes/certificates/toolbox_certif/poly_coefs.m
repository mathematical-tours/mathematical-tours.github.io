function M = poly_coefs(P,K)

% poly_coefs - extract polynomial coefficient in 2 variables.
%
%   M = poly_coefs(P,K);
%
%   M has size (K+1,K+1)
%   M(i+1,j+1) is the coefficient in P of x^i*y^j for 0<=i,j<=K
%
%   K should be larger than the degree of P.
%
%   Copyright (c) 2017 Clarice Poon and Gabriel Peyre

syms x y;

M = zeros(K+1,K+1);
for i=0:K
    for j=0:K
        % extract coefficient of order (a,b) in a polynomial
        M(i+1,j+1) = subs( diff(diff(P,x,i),y,j),{x,y},{0,0}) / (factorial(i)*factorial(j));
    end
end

end