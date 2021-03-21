function [C,E] = glasso(C_hat,lambda)

% projection on |u|<r
projinf = @(u,r)min(max(u,-r),r);
dotp = @(x,y)sum(x(:).*y(:));

p = size(C_hat,1);

C = C_hat + lambda*eye(p);
niter = 1000;
niter_pg = 10;
E = [];
for it=1:niter
    % number of the row
    i = mod(it-1,p)+1; I = 1:p; I(i)=[];
    Q = inv(C(I,I));
    c = C(I,i);
    c_hat = C_hat(I,i);
    % solve
    %     min_{|c-c_hat|_inf<=lambda} 1/2 * <c,Q*c>
    % use projected gradient descent
    tau = 1.5 / norm(Q);
    Epg = [];
    for jt=1:niter_pg
        c = c - tau * Q*c;
        c = projinf(c-c_hat,lambda) + c_hat;
        Epg(end+1) = 1/2*dotp(c,Q*c);
    end
    C(i, I) = c;
    C(I, i) = c;
    % record global energy 
    E(end+1) = -mylog(det(C));
end

end
