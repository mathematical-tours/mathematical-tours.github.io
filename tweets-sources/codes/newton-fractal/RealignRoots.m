function P1 = RealignRoots(P,Q)

n = length(P);
W = abs( repmat(P(:), [1 n]) - repmat(transpose(Q(:)), [n,1]) ).^2;

[I,COST] = Hungarian(W);
P1 = I'*P;

end