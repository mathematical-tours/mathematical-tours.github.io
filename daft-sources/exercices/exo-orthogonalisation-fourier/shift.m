function y = shift(x, s)

n = length(x);

y = [x(s+1:n);x(1:s)];