function r = add(p,q)
m = length(p); n = length(q);
if( m>=n )
    r = p + [q;zeros(m-n,1)];
else
    r = q + [p;zeros(n-m,1)];
end