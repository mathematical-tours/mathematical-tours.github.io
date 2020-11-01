function y = reconstruct_ondelettes(x,psi)
p = length(psi); c = p*dot(psi,psi); y = zeros(p,1);
for( a=1:p-1 ) for( b=0:p-1 )
   ordre = mod(invmod(a,p)*((0:p-1)-b), p)+1; 
   y = y + x(a,b+1)*psi(ordre);
end; end;
y = y/c;