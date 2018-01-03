function y = Walsh(x,nmax)

x = x(:);
if nargin<2
    nmax = 1;
end
if length(x)<=nmax
    y = x; return; 
end


M = @(a,b)(a+b)/sqrt(2);
D = @(a,b)(a-b)/sqrt(2);

u = Walsh(x(1:end/2),nmax);
v = Walsh(x(end/2+1:end),nmax);

y = [ M(u,v);  D(u,v) ];
    


end
