function G = poly_deriv(F,P)

% poly_deriv - apply a polynomial derivative to a function
%
%   G = poly_deriv(F,P);
%
%   G = P(partial_x,partial_y)[F]
%
%   Copyright (c) 2017 Clarice Poon and Gabriel Peyre

syms x y;

[c,T] = coeffs(P);
G = 0;
for k=1:length(c)
    a = log2(subs(T(k), {x y}, {2 1}));
    b = log2(subs(T(k), {x y}, {1 2}));
    G = G + c(k)*diff(diff(F,x,a),y,b);
end

% G = simplify(G);


end


