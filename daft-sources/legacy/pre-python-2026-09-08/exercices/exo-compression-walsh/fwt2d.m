function y = fwt2d(x)
y = x; n = length(x);
for i=1:n
    y(i,:) = fwt(y(i,:)')';
end
for j=1:n
    y(:,j) = fwt(y(:,j));
end