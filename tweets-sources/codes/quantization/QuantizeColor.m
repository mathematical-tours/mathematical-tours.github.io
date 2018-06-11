function f1 = QuantizeColor(f0,K)

n = size(f0,1);
X = reshape(f0, [n*n 3])';
rand('state', 123);
C = X(:,round(rand*n*n));
for i=1:K-1
    D = distmat(X,C);
    D = min(D,[],2);
    [~,r] = max(D,[],1);
    C(:,end+1) = X(:,r);
end

niter = 30;
for i=1:niter
    D = distmat(X,C);
    [D,I] = min(D,[],2); 
    % clf; imagesc( reshape(I,[n n]) ); drawnow;
    for k=1:K
        C(:,k) = mean(X(:,I==k), 2);
    end
end

% generate quantized image
f1 = zeros(n,n,3);
for s=1:3
    f1(:,:,s) = reshape(C(s,I), [n n]);
end

end