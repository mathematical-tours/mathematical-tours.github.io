function y = htgen(x,lambda)
% transforméee de Harley generalisee
n = length(x);
y = zeros(n,1);

for i=0:n-1
    z = cos( 2*pi/n*i*(0:n-1)' + lambda );
    y(i+1) = dot(x,z);
end