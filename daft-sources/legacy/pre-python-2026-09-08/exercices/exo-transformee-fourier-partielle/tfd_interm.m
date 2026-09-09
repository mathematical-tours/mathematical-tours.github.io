function y = tfd_interm(x,alpha)
n = length(x);
f = (0:n-1)'*(0:n-1);
Omega = exp(2i*f*pi/n);
[V,D] = eig(Omega,'nobalance');
y = V*D^alpha*ctranspose(V)*x;