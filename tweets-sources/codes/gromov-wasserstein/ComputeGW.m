function [gamma,GWlist] = ComputeGW(CX,CY,gamma)


n = size(CX,1);
m = size(CY,1);
a = ones(n,1)/n;
b = ones(m,1)/m;
niter = 15;


dotp = @(x,y)sum(x(:).*y(:));
GW = @(gamma)dotp(CX.^2*a,a) + dotp(CY.^2*b,b) - 2*dotp(gamma*CY,CX*gamma);

GWlist = GW(gamma);
for it=1:niter
    [~,gamma1] = mexEMD(a,b,-CX*gamma*CY);
    if norm(gamma1-gamma,'fro')<1e-6
        break;
    end
    gamma=gamma1;
    GWlist(end+1) = GW(gamma);   
end
    

end
