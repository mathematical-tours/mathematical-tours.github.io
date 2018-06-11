
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
% Starting points for the distance computation.


% #starting points
p = 10;
X = floor(rand*n)+1;
for i=2:p
    progressbar(i,p);
    options.start_points = X;
    [D,S,Q] = perform_fast_marching_mesh(V, F, X, options);
    D(D==Inf)=0; 
    [~,X(end+1)] = max(D);
    % display
    clf;
    plot_fast_marching_mesh(V,F, rescale(D), [], options);
    colormap parula(256);
    camlight;
    saveas(gcf, [rep 'seeding-' znum2str(i,2) '.png'], 'png');
end


%%
% Display propagation

q = 50; % #images
nblist = round( linspace(.01,1,q)*n );
for i = 1:length(nblist);
    options.nb_iter_max = nblist(i);
    [D,S,Q] = perform_fast_marching_mesh(V, F, X, options);
    col = D; col(col==Inf) = 0;
    clf;
    hold('on');
    options.face_vertex_color = rescale(col);
    plot_mesh(V,F,options);
    plot3( V(1,X), V(2,X), V(3,X), 'r.', 'MarkerSize', 25);
    colormap jet(256);
    camlight; drawnow;
    saveas(gcf, [rep 'fmm-' num2str(p) '-' znum2str(i,2) '.png'], 'png');
end

% AutoCrop(rep, 'fmm-');