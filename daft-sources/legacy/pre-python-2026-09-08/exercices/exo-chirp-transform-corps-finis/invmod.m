function y = invmod(x,p)
[u,y,d] = gcd(x,p); y = mod(y,p);