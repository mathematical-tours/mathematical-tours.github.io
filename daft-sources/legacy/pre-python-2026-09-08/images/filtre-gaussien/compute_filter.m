function f = compute_filter(n,s)

% calcule un filtre gaussien de taille n et de variance s.

x = -1:2/(n-1):1;
[X,Y] = MESHGRID(x,x);

f = exp( -(X.^2+Y.^2)/(2*s) );

f = f / sum(sum(f));