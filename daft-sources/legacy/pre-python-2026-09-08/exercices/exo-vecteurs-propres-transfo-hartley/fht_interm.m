function y = fht_interm(x,lambda)
N = length(x);
u1 = sqrt(N)*x+fht(x); u2 = sqrt(N)*x-fht(x);
y = ( sqrt(N)^lambda )*( u1 + (-1)^lambda*u2 );