%%
% Test for the Durand?Kerner method.

rep = 'results/';
[~,~] = mkdir(rep);
if not(exist('test'))
    test=1;
end

% n = 5;
% z0 = randn(n,1)+1i*randn(n,1);

% click and play
clf; hold on;
z0 = [];
while true
    axis([-1 1 -1 1]);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);   
    if button==3
        break;
    end
    z0(end+1) = a+1i*b;
end
z0 = z0(:);
n = length(z0);
% click and play
clf; hold on;
plot(z0, 'k.', 'MarkerSize', 25);
z1 = [];
for i=1:n
    axis([-1 1 -1 1]);
    [a,b,button] = ginput(1);
    plot(a,b, '.', 'MarkerSize', 15);   
    z1(end+1) = a+1i*b;
end
z1 = z1(:);


niter = 600;

tau = .02;

% z = randn(n,1)+1i*randn(n,1);
Z = [];
z = z1;
for i=1:niter
    Z(:,end+1) = z;
    q = z;
    for k=1:n
        J = setdiff(1:n,k);
        q(k) = z(k) - polyval(z(k),z0) ./ polyval(z(k),z(J));
    end
    z = (1-tau)*z + tau*q;
    
    D = abs( repmat(z,[1 length(z0)]) - repmat(permute(z0,[2 1]),[length(z) 1]) );
	d = mean(min(D));
    if d<=1e-2
        Z(:,end+1) = z;
        break;
    end
end

clf; hold on;
for i=1:size(Z,2)-1
    t = (i-1)/(size(Z,2)-2);
    plot( permute(Z(:,i:i+1), [2 1]), '-', 'Color', [t 0 1-t], 'LineWidth', 2 );
end
plot(z0, 'r.', 'MarkerSize', 25);
axis tight; axis equal; % axis([-1 1 -1 1]); 
axis off;

saveas(gcf, [rep 'roots-' num2str(test) '.eps'], 'epsc');
test = test+1;