function y = fft_gt(x,p)
n = length(x); q = n/p; 
y = zeros(n,1); m = zeros(p,q);
for i=0:n-1
    m(mod(i,p)+1, mod(i,q)+1) = x( mod(-i,n)+1 );    
end
m = fft2(m);
for s1 = 0:p-1
for s2 = 0:q-1
    y( mod(s1*q+s2*p,n)+1 ) = m(s1+1,s2+1);
end
end