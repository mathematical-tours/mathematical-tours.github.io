% Transformée de Hartley rapide
function y = fht(x)

n = length(x);
y = zeros(n,1);

for i=0:n-1
    z = cos( 2*pi/n*i*(0:n-1)' ) + sin( 2*pi/n*i*(0:n-1)' );
    y(i+1) = dot(x,z);
end