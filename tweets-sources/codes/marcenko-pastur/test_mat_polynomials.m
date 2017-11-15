n = 1000;
S = @(x)(x+x')/2;
A = S(randn(n));
B = S(randn(n)); 
C = S(randn(n)); 

hist(eig(A*B+B*A), 100)
hist(eig(A*B*B*A + A*C*C*A + 2*B), 100);