function [I,J,w] = myot(X,Y,epsilon)

tol = 1e-10;

nX = length(X);
nY = length(Y);
c = abs(X-transpose(Y)).^2;
a = ones(nX,1)/nX;
b = ones(nY,1)/nY;
if epsilon>0
    options.niter = 50000;
    options.tol = 1e-8;
    [P,f,g,Err] = sinkhorn(c,a,b,epsilon,options);
    P = P/max(P(:));
    P = P .* (P>tol);
    [I,J,w] = find( P );
else
    [cost,gamma] = mexEMD(a,b,c);
    [I,J,w] = find(gamma);
end
w = w/max(w);

end