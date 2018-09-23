function f = smooth3d(f,K)

n = size(f,1);
a = [1 1:n-1]; b = [2:n n];
for k=1:K
    f = f + ...
        f(a,:,:) + f(b,:,:) + ...
        f(:,a,:) + f(:,b,:) + ...
        f(:,:,a) + f(:,:,b);
    f = f/7;
end

end