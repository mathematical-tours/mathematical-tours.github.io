function X = ista(A,y,x,lambda,tau, niter)


Soft = @(x,tau)max(abs(x)-tau,0).*sign(x);

X = [x];
for i=1:niter-1
    x = Soft(x - tau*A'*(A*x-y), tau*lambda);
    X(:,end+1) = x;
end

end