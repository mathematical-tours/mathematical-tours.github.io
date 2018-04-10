function P = Lagrange(Q,x,i)

n = length(Q);
if nargin==2
    P = [];
    for i=1:n
        P(:,i) = Lagrange(Q,x,i);
    end
    return;
end

P = x*0+1;
for j = setdiff(1:n,i)
    P = P .* (x-Q(j))/(Q(i)-Q(j));
end



end