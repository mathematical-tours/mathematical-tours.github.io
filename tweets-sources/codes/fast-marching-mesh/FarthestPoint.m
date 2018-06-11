%%
% Farthest point seeding

name = 'elephant';
addpath('../toolbox/');
rep = MkResRep(name);

% Load a 3D mesh.
[V,F] = read_off([name '.off']);
n = size(V,2);

options.end_points = [];
options.nb_iter_max = Inf;
% Use a uniform, constant, metric for the propagation.
options.W = ones(n,1);


%%
% Seeding

% #starting points
p = 50;
C = distinguishable_colors(p+1); C(4,:) = [];
X = floor(rand*n)+1;
for i=1:p
    progressbar(i,p);
    options.start_points = X;
    [D,S,Q] = perform_fast_marching_mesh(V, F, X, options);
    %
    [Qexact,DQ, voronoi_edges] = compute_voronoi_mesh(V, F, X, options);
    options.voronoi_edges = voronoi_edges;    
    % display distance
    clf;
    options.col = D;
    plot_fast_marching_mesh(V,F, D, [], options);
    colormap parula(256);
    camlight;
    saveas(gcf, [rep 'seeding-' znum2str(i,2) '.png'], 'png');    
    % display Voronoi
    col = zeros(3,n);
    for k=1:i
        I = find(Q==X(k));
        col(:,I) = repmat(C(k,:)', [1 length(I)]);
    end
    options.col = col';
    clf; hold on;
    plot_fast_marching_mesh(V,F, D, [], options); 
    camlight;
    saveas(gcf, [rep 'voronoi-' znum2str(i,2) '.png'], 'png');
    % next point
    D(D==Inf)=0; 
    [~,X(end+1)] = max(D);
end
