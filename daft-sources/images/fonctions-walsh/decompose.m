function y = decompose(x,n)

y = x.*0;
for i=1:n
    y = y + walsh(i)*dot(walsh(i),x);
end

y = y/length(x);