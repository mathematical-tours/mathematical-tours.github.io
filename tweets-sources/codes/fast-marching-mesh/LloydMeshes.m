%%
% Check for Lloyd algorithm on meshes.

name = 'elephant';
addpath('../toolbox/');
rep = MkResRep('lloyd');

% Load a 3D mesh.
[V,F] = read_off([name '.off']);
n = size(V,2);

options.end_points = [];
options.nb_iter_max = Inf;
% Use a uniform, constant, metric for the propagation.
options.W = ones(n,1);

% #starting points
p = 50;
X = floor(rand(p,1)*n)+1;

% choose clustered points
i = floor(rand(1)*n)+1;
i = 16095;
D = sum( (V - repmat(V(:,i), [1 n]) ).^2 );
[~,I] = sort(D, 'ascend');
I = I(1:round(.05*end)); I = I(randperm(length(I)));
X = I(1:p);

options.start_points = X;
[D,S,Q] = perform_fast_marching_mesh(V, F, X, options);
%
options.voronoi_edges = [];
options.col = Q;
clf; hold on;
plot_fast_marching_mesh(V,F, D, [], options);
camlight;

% 
C = distinguishable_colors(p+1); C(4,:) = [];
q = 80;
for i=1:q
    options.start_points = X;
    [D,S,Q] = perform_fast_marching_mesh(V, F, X, options);
    %
    [Qexact,DQ, voronoi_edges] = compute_voronoi_mesh(V, F, X, options);
    options.voronoi_edges = voronoi_edges;
    % setup colors
    col = zeros(3,n);
    for k=1:p
        I = find(Q==X(k));
        col(:,I) = repmat(C(k,:)', [1 length(I)]);
    end
    options.col = col';
    %
    clf; hold on;
    plot_fast_marching_mesh(V,F, D, [], options); 
    camlight;
    saveas(gcf, [rep 'lloyd-' num2str(p) '-' znum2str(i,2) '.png'], 'png');
    % update point locations
    for k=1:p
        I = find(Q==X(k)); r = length(I);
        m = mean( V(:,I), 2 );
        D = sum( (V(:,I) - repmat(m, [1 r]) ).^2 );
        [~,l] = min(D);
        X(k) = I(l);
    end
end
% AutoCrop(rep, 'lloyd-');