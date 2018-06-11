%%
% Test for Farthest point sampling.


addpath('../toolbox/');
rep = MkResRep();

n = 256;
name = 'hibiscus';
f = rescale(sum(load_image(name, n),3));
[fS,I] = sort(f(:)); f(I) = linspace(0,1,length(f(:)));

% points
P = [1;1];

[Y,X] = meshgrid(1:n,1:n);
Q = [X(:)';Y(:)'];

fmin = .05;
U = rescale(-f,0.01,1).^2;


pmax = 1000;
for i=1:pmax
    if mod(i,3)==1
    clf; hold on;
    imagesc(1:n,1:n,f);
    plot(P(2,:), P(1,:), 'b.', 'MarkerSize', 20);
    axis image; axis off;
    drawnow;
    end
    %
    D = distmat(P,Q); 
    if 0
    w = U(P(1,:) + (n-1)*P(2,:))+fmin;
    D = D .* repmat(w(:), [1 n*n]);
    else
    D = D .* repmat(U(:)', [size(D,1) 1]);        
    end
    D = min(D,[],1);
    [~,j] = max(D(:));
    P(:,end+1) = [X(j);Y(j)];
end