%%
% Test for Floyd Warshall algorithm.

rep = '../results/floyd-warshall/';
[~,~] = mkdir(rep);
addpath('../toolbox/');

%% 
% Set D=inf for non-connected vertices.

% image graph
n = 10;
[Y,X] = meshgrid(1:n,1:n);
P = [X(:)'; Y(:)'];
D = distmat(P,P);

% 1D graph
n = 20;
D = abs( repmat((1:n)', [1 n]) - repmat((1:n), [n 1]) );

D(D>1)=Inf; %  | D<1
niter = size(D, 1);


i = n/2;
j = n/2;
s = sub2ind([n n],i,j);

G = D;
niter = size(D, 1);
for k=1:niter
    i2k = repmat(G(:,k), 1, niter);
    k2j = repmat(G(k,:), niter, 1);
    G = min(G, i2k+k2j);
    % clf; imagesc(reshape(G(s,:), [n n]));
    clf;  imagesc(G); drawnow; axis image; axis off;
    if mod(k,4)==1
        saveas(gcf, [rep 'fw-' num2str(k) '.png']);
    end
end
saveas(gcf, [rep 'fw-' num2str(k) '.png']);


G = D;
niter = size(D, 1);
for k=1:niter
    G = MinPlusMult(G,D);
    clf;  imagesc(G); drawnow; axis image; axis off;
    if mod(k,4)==1
        saveas(gcf, [rep 'mp-' num2str(k) '.png']);
    end
end
saveas(gcf, [rep 'mp-' num2str(k) '.png']);