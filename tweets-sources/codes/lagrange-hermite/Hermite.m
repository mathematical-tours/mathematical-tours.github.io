function P = Hermite(Q,x)

n = length(Q);
N = length(x);
x = x(:);

R = Lagrange(Q,x).^2;

% compute derivatives of R
U = 2 ./ ( repmat(Q,[1 n]) - repmat(Q',[n 1]) );
U(U==Inf) = 0;
dR = sum(U,2);


P = zeros(N,2*n);
for i = 1:n
    P(:,i)   = R(:,i) .* (1 - (x-Q(i))*dR(i) );
    P(:,i+n) = R(:,i) .* (x-Q(i));
end



end