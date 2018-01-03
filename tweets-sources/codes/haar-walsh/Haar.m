function y = Haar(x,nmax)

x = x(:);
if nargin<2
    nmax = 1;
end
if length(x)<=nmax
    y = x; return; 
end

M = @(a,b)(a+b)/sqrt(2);
D = @(a,b)(a-b)/sqrt(2);

y = [ Haar( M(x(1:2:end),x(2:2:end)) ,nmax ) ; D(x(1:2:end),x(2:2:end))   ];


end
