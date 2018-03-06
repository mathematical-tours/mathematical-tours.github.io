function C = MinPlusMult(A,B)

n = size(A,1);
m = size(B,2);
C = zeros(n,m);
for i=1:n
    for j=1:m
        C(i,j) = min( A(i,:) + B(:,j)' );
    end
end

end