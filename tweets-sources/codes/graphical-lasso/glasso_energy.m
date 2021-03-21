function [E,gradE] = glasso_energy(P,C,lambda)

dotp = @(x,y)sum(x(:).*y(:));
E = dotp(P,C)-mylog(det(P)) + lambda*norm(P-diag(diag(P)),1);
if  min(eig(P))<0
    P = Inf;
end

if nargin>1
    gradE = C-inv(P);
end

end