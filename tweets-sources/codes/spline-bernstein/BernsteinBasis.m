function B = BernsteinBasis(x,n)

x = x(:);
m = length(x);
B = ones(m,1);
for i=1:n-1
    B = [zeros(m,1),B,zeros(m,1)];
    B = B(:,1:end-1).*(1-x) + B(:,2:end).*x;
end

end
