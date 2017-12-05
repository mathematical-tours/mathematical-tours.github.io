function F = plot_isosurface(M,options)

% plot_isosurface - display an isosurface
%
% F = plot_isosurface(M,options)
%
%   M is a 3D matrix
%   display the iso-level at options.isolevel
%   F is the patch object of the surface.
%
%   Copyright (c) Gabriel Peyre 2015

clf;

n = max(size(M));
t = linspace(-1,1,n);

options.null = 0;
rho = getoptions(options, 'isolevel', (max(M(:))-min(M(:)))/2  );
t = getoptions(options, 't', []);
alpha_val = getoptions(options, 'alpha', 1);
color = getoptions(options, 'color', 'red');


F = isosurface( t,t,t,M,rho );
p = patch(F);
isonormals( t,t,t,M,p );
set(p, 'FaceColor', color, 'EdgeColor', 'none'); 

set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'ZTick', []);

lighting phong;
alpha(alpha_val);    
camproj('perspective');
view(3);

% axis equal;
% axis([1 size(W,1) 1 size(W,2) 1 size(W,3)]);
% daspect([1 1 1]);
% cameramenu;

axis equal; axis off;
camlight;