function y = MLS(xi,yi,x,d,a)

n = length(x);
y = zeros(n,1);
for k=1:n
    w = 1./(1e-8 + abs(xi-x(k))).^a;
    P = polyfitweighted(xi,yi,d,w);    
    y(k) = polyval(P,x(k));
end

end
