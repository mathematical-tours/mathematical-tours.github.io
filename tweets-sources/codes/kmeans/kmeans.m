%%
% Test for kmeans.


addpath('../toolbox/');
rep = MkResRep();

n = 10000; % for compute
n1 = 4000; % for display
% mixture of gaussians
%
m = [.3+.2i, .2+.6i, .6+.8i, .8+.3i ];
s = [.2 .15 .4 .3]/3;
%
m = .5+.5i; s = .15;
%
p = round(n/length(m));
gauss = @(m,s)s*randn(p,1)+1i*s*randn(p,1) + m;
X = [];
for k=1:length(m)
    X = [X; gauss(m(k),s(k))];
end
X = rand(n,1)+1i*rand(n,1);
X = X(randperm(n));


clf;
plot(X, 'k.', 'MarkerSize', 15);
axis equal; axis([0 1 0 1]);

% centroids
K = 30; 

% random init
c = .5+.5i + .03*( randn(K,1)+1i*randn(K,1) );
c = .03*( rand(K,1)+1i*rand(K,1) );

c = rand+1i*rand; 
for k=2:K
    D = min(abs(X-transpose(c)), [], 2);
    P = D.^2/sum(D.^2); % probability to pick a point
    i = rand_discr(P, 1);
    c(end+1) = X(i); c=c(:);
end

% k-means++ init
X = [X;c];
Col = distinguishable_colors(K);

% kemeans
q = 50; 
tau = 1; % to slow down evolution
for it=1:q
    % NN assignement
    D = abs(X-transpose(c));
    [tmp,I] = min(D, [], 2);
    % display
    clf; hold on;
    for k=1:K
        L = find(I==k);
        plot( X(L(L<=n1)), '.', 'Color', Col(k,:), 'MarkerSize', 15 );
        if not(isempty(L))
            c(k) = mean(X(L));
        end
    end
    axis equal; axis([0 1 0 1]);
    set(gca, 'XTick', [], 'YTick', []); box on;
    drawnow;
    %
    saveas(gcf, [rep 'anim-' znum2str(it,2) '.png']);
end