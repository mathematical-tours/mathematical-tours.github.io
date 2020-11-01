function y = decompose(x,n)

global N;

y = x.*0;
for i=1:n
    y = y + walsh2d(i)*sum(dot(walsh2d(i),x));
end

y = y/N^2;