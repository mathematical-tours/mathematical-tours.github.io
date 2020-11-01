function y = fht2d(x)
n = length(x); y = zeros(n,n);
for(i=1:n) y(i,:) = fht(x(i,:)')'; end;
for(j=1:n) y(:,j) = fht(y(:,j)); end;