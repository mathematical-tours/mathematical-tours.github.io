function gamma0 = ComputeLB(CX,CY)

n = size(CX,1);
m = size(CY,1);
a = ones(n,1)/n;
b = ones(m,1)/m;

CX1 = sort(CX,2);
CY1 = sort(CY,2);
D = reshape(CX1, [n 1 n]) - reshape(CY1, [1 m m]);
D = sum(abs(D),3);
[cost,gamma0] = mexEMD(a,b,D);

end
