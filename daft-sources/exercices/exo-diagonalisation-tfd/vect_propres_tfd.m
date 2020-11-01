function V = vect_propres_tfd(n)
x = (0:n-1)'*(0:n-1);
Omega = 1/sqrt(n)*exp(2i*x*pi/n);
d = 2*( cos(2*pi/n*(0:n-1)) - 2 );
S = diag(d,0) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
S(1,n) = 1; S(n,1) = 1;
[V,D] = eig(S);