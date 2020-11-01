function y = transfo_ondelettes(x,psi)
p = length(psi); y = zeros(p-1,p);
for( a=1:p-1 ) for( b=0:p-1 )
   ordre = mod(invmod(a,p)*((0:p-1)-b), p)+1; 
   y(a,b+1) = dot(x,psi(ordre));
end; end;