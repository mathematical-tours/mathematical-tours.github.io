function EtaV = compute_etav(C,x0,s0)

% compute_etaw - compute eta_V certificate in 1-D and 2-D
%
%   EtaW = compute_etav(C,x0);
%
%   C(x,y,x1,y1) should be a symbolic function.
%   EtaV(x,y) is a symbolic function.
%   Each x0(i,:) is a point.
%
%   For 1-D problems, simply make the above functions not depend on (y,y1).
%
%   Copyright (c) 2017 Clarice Poon and Gabriel Peyre

syms x y x1 y1;

N = size(x0,1);


Cx = diff(C,x);
Cy = diff(C,y);
Cxx1 = diff(Cx,x1);
Cyy1 = diff(Cy,y1);
% cross
Cyx1 = diff(Cy,x1);
Cxy1 = diff(Cx,y1);


evalij = @(F,i,j)double( subs(F,{x,y,x1,y1},{x0(i,1),x0(i,2),x0(j,1),x0(j,2)}) );

% build the Gamma matrix
Gamma = zeros(3*N,3*N);
for i=1:N
    for j=1:N
        Gamma(i,j) = evalij(C,i,j);
        Gamma(i+N,j+N) = evalij(Cxx1,i,j);
        Gamma(i+2*N,j+2*N) = evalij(Cyy1,i,j);
        Gamma(i+N,j) = evalij(Cx,i,j); Gamma(j,i+N) = Gamma(i+N,j);
        Gamma(i+2*N,j) = evalij(Cy,i,j); Gamma(j,i+2*N) = Gamma(i+2*N,j);
        Gamma(i+2*N,j+N) = evalij(Cyx1,i,j); Gamma(j+N,i+2*N) = Gamma(i+2*N,j+N);   
    end
end

if not(issymmetric(Gamma))
    error('Gamma is not symmetric');
end
if cond(Gamma)>1e6
    warning('Gamma is singular or close to singular.');
end


% coefficients
U = pinv(Gamma)*[s0(:);zeros(2*N,1)];
% build EtaV
subsij = @(F,i)subs(F,{x y x1 y1},{x0(i,1) x0(i,2) x y});
EtaV = 0; 
for i=1:N
    EtaV = EtaV + U(i)*subsij(C,i);
    EtaV = EtaV + U(i+N)*subsij(Cx,i);
    EtaV = EtaV + U(i+2*N)*subsij(Cy,i);
end
EtaV = symfun(EtaV,[x y]);

end