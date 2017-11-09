function EtaF = compute_etaf(C,x0,s0)

% compute_etaf - compute eta_F certificate in 1-D and 2-D
%
%   EtaF = compute_etaf(C,x0,s0);
%
%   C(x,y,x1,y1) should be a symbolic function.
%   EtaF(x,y) is a symbolic function.
%   Each x0(i,:) is a point.
%
%   For 1-D problems, simply make the above functions not depend on (y,y1).
%
%   Copyright (c) 2017 Clarice Poon and Gabriel Peyre

syms x y x1 y1;

N = size(x0,1);



evalij = @(F,i,j)double( subs(F,{x,y,x1,y1},{x0(i,1),x0(i,2),x0(j,1),x0(j,2)}) );

% build the Gamma matrix
Gamma = zeros(N,N);
for i=1:N
    for j=1:N
        Gamma(i,j) = evalij(C,i,j);
    end
end

if not(issymmetric(Gamma))
    error('Gamma is not symmetric');
end
if cond(Gamma)>1e6
    warning('Gamma is singular or close to singular.');
end


% coefficients
U = pinv(Gamma)*s0(:);
% build EtaF
subsij = @(F,i)subs(F,{x y x1 y1},{x0(i,1) x0(i,2) x y});
EtaF = 0; 
for i=1:N
    EtaF = EtaF + U(i)*subsij(C,i);
end
EtaF = symfun(EtaF,[x y]);

end