function y = frft2d(x,alpha)

n = length(x);

for i=1:n
    y(i,:) = frft( x(i,:), alpha );    
end
for j=1:n
    y(:,j) = frft( y(:,j), alpha );    
end