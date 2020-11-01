function r = karatsuba(p,q)
n = length(p)-1;
if(n==0) r=p*q; return; end;
k = floor((n+1)/2);
p0 = p(1:k); p1 = p((k+1):(n+1));
q0 = q(1:k); q1 = q((k+1):(n+1));
r0 = karatsuba(p0,q0); r2 = karatsuba(p1,q1);
r1 = karatsuba(add(p0,p1),add(q0,q1));
r1 = add(r1,-r0); r1 = add(r1,-r2);
r = add( r0, [zeros(k,1);r1] );
r = add( r, [zeros(2*k,1);r2] );