function f = calcul_filtre(n,s)
x = -1:2/(n-1):1;
[X,Y] = meshgrid(x,x);
f = exp( -(X.^2+Y.^2)/(2*s) );
f = f / sum(sum(f));