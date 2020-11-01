function y = transformee_sinus_naif(x)
n = length(x)+1;
y = zeros(n-1,1);
for i=1:n-1
    for j=1:n-1
        y(i) = y(i) + x(j)*sin(i*j*pi/n);
    end
end