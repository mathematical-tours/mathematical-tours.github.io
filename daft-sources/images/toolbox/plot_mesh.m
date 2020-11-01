function plot_mesh(vertex,face,cf,ce,normal)

% plot_mesh - plot a 3D mesh.
%
%   plot_mesh(vertex,face,cf,ce);
%
%   'cf' is face color and 'ce' edge color.
%
%   Copyright (c) 2004 Gabriel Peyré


if nargin<2
    error('Not enough arguments.');
end
if nargin<3
    cf = 0.7;
end
if nargin<4
    ce = 0;
end

patch('vertices',vertex,'faces',face,'facecolor',[cf cf cf],'edgecolor',[ce ce ce]);
lighting phong;
camlight infinite; 
camproj('perspective');
axis square; 
axis off;

if nargin==5
    % plot the normals
    hold on;
    quiver3(vertex(:,1),vertex(:,2),vertex(:,3),normal(:,1),normal(:,2),normal(:,3),0.8);
    hold off;
end