function K = poly_order(N)

% poly_order - optimal order for de Boor expansion
%
%   K = poly_order(N);
%
%   Find smallest K such that 
%       (K+1)*(K+2)/2 >= 3*N
%
%   Copyright (c) 2017 Clarice Poon and Gabriel Peyre

syms KK;
v = solve((KK+1)*(KK+2)/2 - 3*N == 0, KK) ;
K = ceil( double(max(v) ) );

end