function y = decompose_tensoriel(x,A)
s = length(A); m = length(x); m0 = m/s;
if(m==1) y = x; return; end;
B = zeros(m0,s); % résulats temporaires
for j=1:s
    sel = ((j-1)*m0+1):(j*m0);
    B(:,j) = decompose_tensoriel(x(sel),A);
end
y = zeros(m,1);
for i=1:s
    sel = ((i-1)*m0+1):(i*m0);
    for j=1:s
        y(sel) = y(sel) + A(i,j)*B(:,j);    
    end
end